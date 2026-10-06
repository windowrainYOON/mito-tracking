"""Cell ROI extraction from mitochondria images using weak autofluorescence.

Pushing the LUT on a mito image reveals dim cytoplasmic/nuclear autofluorescence
that outlines each cell, while gaps between cells stay near zero. Pipeline:
  1. mask burned-in annotations (scale bar), clip + log-compress so bright mito
     does not dominate, Gaussian smooth -> autofluorescence density map
  2. foreground = density > k * Li threshold, cleaned morphologically
  3. dark "valley" lines (Sato filter, black ridges) mark cell-cell borders
  4. seeds = maxima of distance-to-valley inside foreground -> watershed
  5. region merging: adjacent regions whose shared border is not darker than
     their interiors (nucleus/mito edges, not real cell borders) are merged,
     and undersized fragments are absorbed into their weakest-border neighbour
With --nuclei (same field, nuclear stain), steps 4-5 are replaced: each in-focus
nucleus (Otsu on the nuclear channel; dim out-of-focus nuclei fall below it) is
a watershed marker, so every cell ROI holds exactly one nucleus.
QC (nucleus mode): a cell is excluded as `edge` if its ROI touches the frame.
Neighbours with no dark border between them (continuous-looking cytoplasm) are
listed in `weak_border_with`; they are still split one-nucleus-one-cell unless
--exclude-binucleate is given, which drops them as `binucleate`.
Outputs ImageJ RoiSet.zip (all) + RoiSet_filtered.zip (QC pass only), label TIFF, CSV and overlay PNGs.
"""
import argparse, csv, os
import numpy as np, tifffile, roifile
from scipy import ndimage as ndi
from skimage import exposure, feature, filters, io, measure, morphology, segmentation


REF_PX_UM = 0.0990       # pixel size of the reference data set the size parameters were tuned on
REF_GREEN_UNIT = 4.79    # 75th pct of the background-subtracted, sigma=4 smoothed green image in that set
REF_GREEN_OFFSET = 0.077  # 1st pct of the same smoothed image (dark gaps between cells)


def normalize_green(g, px_scale=1.0):
    """Rescale green so its cytoplasmic level matches the reference set.

    Cell outlines come from dim autofluorescence/POI signal, and the segmentation constants
    (clip level, log offset) are in reference intensity units. Using the dark-gap level as the
    offset and the 75th percentile of the smoothed image as the unit makes the result independent
    of gain, exposure, bit depth and display stretch.
    """
    sm = ndi.gaussian_filter(g, 4 * px_scale)
    off = np.percentile(sm, 1)
    unit = np.percentile(sm, 75) - off
    if unit <= 0:
        raise ValueError('Green channel has no usable signal for cell segmentation')
    gain = REF_GREEN_UNIT / unit
    gn = np.clip((g - off) * gain + REF_GREEN_OFFSET, 0, None)
    return gn, dict(offset=float(off), gain=float(gain))


def load_green(path):
    a = tifffile.imread(path)
    if a.ndim == 3 and a.shape[-1] == 3:
        g = a[..., 1].astype(float)
        overlay = (a[..., 0] > 200) & (a[..., 2] > 200)  # white annotations (scale bar)
        g[ndi.binary_dilation(overlay, iterations=2)] = 0
    else:
        g = a.astype(float)
    with tifffile.TiffFile(path) as t:
        xr = t.pages[0].tags.get('XResolution')
        px_um = xr.value[1] / xr.value[0] if xr else 1.0
    return a, g, px_um


def segment_nuclei(b, min_area_px=3000, k=1.0):
    """Otsu nuclei; `k` = pixel-size scale factor (REF_PX_UM / pixel size)."""
    s = filters.gaussian(b, 2 * k, preserve_range=True)
    m = s > filters.threshold_otsu(s)
    m = ndi.binary_fill_holes(morphology.binary_opening(m, morphology.disk(max(1, round(3 * k)))))
    m = morphology.remove_small_objects(m, max_size=int(min_area_px * k * k))
    return measure.label(m)


def border_stats(lab, img):
    """Mean img along each shared border between adjacent labels, and border length."""
    s, n = {}, {}
    for a, b, ia, ib in ((lab[:, :-1], lab[:, 1:], img[:, :-1], img[:, 1:]),
                         (lab[:-1, :], lab[1:, :], img[:-1, :], img[1:, :])):
        m = (a != b) & (a > 0) & (b > 0)
        p, q = np.minimum(a[m], b[m]), np.maximum(a[m], b[m])
        v = (ia[m] + ib[m]) / 2
        for k, val in zip(zip(p.tolist(), q.tolist()), v.tolist()):
            s[k] = s.get(k, 0.0) + val; n[k] = n.get(k, 0) + 1
    return {k: (s[k] / n[k], n[k]) for k in s}


def merge_regions(lab, img, contrast_tau, min_area, min_border=30):
    lab = lab.copy()
    while True:
        props = {r.label: r for r in measure.regionprops(lab, intensity_image=img)}
        med = {l: np.median(img[lab == l]) for l in props}
        edges = border_stats(lab, img)
        best = None
        for (p, q), (bmean, blen) in edges.items():
            if blen < min_border:
                continue
            contrast = min(med[p], med[q]) - bmean  # how much darker the border is
            small = min(props[p].area, props[q].area) < min_area
            if contrast < contrast_tau or small:
                key = (0 if small else 1, contrast)
                if best is None or key < best[0]:
                    best = (key, p, q)
        if best is None:
            return segmentation.relabel_sequential(lab)[0]
        _, p, q = best
        lab[lab == q] = p


def segment(g, sigma=10, fg_k=0.75, valley_pct=88, seed_dist=90,
            contrast_tau=0.08, min_area_px=25000, nuclei=None, k=1.0):
    """Lengths/areas are in reference pixels and scaled by `k` (REF_PX_UM / pixel size).
    `g` should already be intensity-normalized (normalize_green)."""
    disk = lambda r: morphology.disk(max(1, round(r * k)))
    L = np.log1p(np.clip(g, 0, 20))
    dens = filters.gaussian(L, sigma * k)
    fg = dens > filters.threshold_li(dens) * fg_k
    fg = morphology.binary_opening(fg, disk(15))
    fg = morphology.remove_small_objects(fg, max_size=int(8000 * k * k))
    fg = morphology.remove_small_holes(fg, max_size=int(20000 * k * k))

    d4 = filters.gaussian(L, 4 * k)
    v = filters.sato(d4, sigmas=[10 * k, 16 * k], black_ridges=True)
    vn = v / np.percentile(v, 99.5)
    valley = vn > np.percentile(vn[fg], valley_pct)

    if nuclei is not None:
        fg = fg | (nuclei > 0)
        dn = dens / np.percentile(dens[fg], 99)
        land = vn + 0.5 * (1 - np.clip(dn, 0, 1))
        lab = segmentation.watershed(land, nuclei, mask=fg)
        pk = np.array([r.centroid for r in measure.regionprops(nuclei)]).astype(int)
        return dict(dens=dens, fg=fg, valley=valley, seeds=pk, raw=lab, nuclei=nuclei, landscape=land,
                    labels=smooth_labels(lab, nuclei, k))

    dist = filters.gaussian(ndi.distance_transform_edt(fg & ~valley), 3)
    pk = feature.peak_local_max(dist, min_distance=seed_dist, threshold_abs=25,
                                labels=measure.label(fg), exclude_border=False)
    markers = np.zeros(fg.shape, int)
    markers[tuple(pk.T)] = np.arange(1, len(pk) + 1)
    lab = segmentation.watershed(vn - 0.5 * dist / dist.max(), markers, mask=fg)
    raw_lab = lab.copy()
    lab = merge_regions(lab, dens, contrast_tau, min_area_px)
    return dict(dens=dens, fg=fg, valley=valley, seeds=pk, raw=raw_lab, nuclei=None,
                labels=smooth_labels(lab))


def find_dim_nuclei(nuc_img, nuclei, px, k=1.0, rel=0.4, min_frac=0.4, min_solidity=0.8):
    """Out-of-focus nuclei: above `rel` x the in-focus Otsu level, compact, at least `min_frac` of the
    minimum in-focus nucleus area and clear of the in-focus nuclei. Returned as labels numbered after
    the in-focus ones (0 elsewhere). They give their cell a seed so its mitochondria are not handed to
    a neighbour."""
    s = filters.gaussian(nuc_img, 2 * k, preserve_range=True)
    t = filters.threshold_otsu(s)
    m = (s > rel * t) & ~ndi.binary_dilation(nuclei > 0, iterations=max(1, round(1.5 / px)))
    m = morphology.binary_opening(m, morphology.disk(max(1, round(3 * k))))
    out = np.zeros_like(nuclei)
    nxt = int(nuclei.max())
    for r in measure.regionprops(measure.label(m)):
        if r.area >= min_frac * 3000 * k * k and r.solidity >= min_solidity:
            nxt += 1
            out[tuple(r.coords.T)] = nxt
    return out


def segment_morph(red, nuclei, landscape, px, k=1.0, fg_level=0.12, compactness=0.003, mito_weight=1.0):
    """Cell ROIs from shape and intensity: one nucleus at the centre of each cell, its cytoplasm spread
    around it, and neighbouring cells separated by thin lines where the mitochondria stop.

      foreground : mitochondria density (red, sigma 3 um) above `fg_level` of its 99th percentile, + nuclei
      cost       : low where mitochondria are dense (sigma 1.5 um), high on the dark mito-free lines
                   (Sato black-ridge filter at 1.5-2.5 um), + the green-autofluorescence border map
      watershed  : from every nucleus (in focus and dim), compact (`compactness` per reference pixel of
                   distance from the nucleus), so each cell grows evenly around its nucleus
    """
    r = red.astype(float)
    md = ndi.gaussian_filter(r, 1.5 / px)
    md = np.clip(md / max(np.percentile(md, 99), 1e-9), 0, 1)
    ridge = filters.sato(md, sigmas=[1.5 / px, 2.5 / px], black_ridges=True)
    ridge = np.clip(ridge / max(np.percentile(ridge, 99.5), 1e-9), 0, 2)
    md3 = ndi.gaussian_filter(r, 3 / px)
    fg = (md3 > fg_level * max(np.percentile(md3, 99), 1e-9)) | (nuclei > 0)
    fg = morphology.binary_closing(fg, morphology.disk(max(1, round(2 / px))))
    fg = morphology.remove_small_holes(fg, max_size=int(200 / px / px))
    cost = mito_weight * (1 - md) + mito_weight * ridge + landscape
    lab = segmentation.watershed(cost, nuclei, mask=fg, compactness=compactness / k)
    return dict(labels=smooth_labels(lab, nuclei, k), fg=fg, mito_density=md, ridge=ridge)


def smooth_labels(lab, keep=None, k=1.0):
    """Open each ROI, fill holes, keep its largest piece; `keep` pixels (nuclei) are never lost."""
    out = np.zeros_like(lab)
    for l in range(1, lab.max() + 1):
        m = morphology.binary_opening(lab == l, morphology.disk(max(1, round(9 * k))))
        if keep is not None:
            m |= keep == l
        m = ndi.binary_fill_holes(m)
        if m.any():
            cc = measure.label(m); m = cc == np.argmax(np.bincount(cc.ravel())[1:]) + 1
            out[m & (out == 0)] = l
    return segmentation.relabel_sequential(out)[0]


def qc_status(lab, dens, touches, binuc_tau, exclude_binuc=False, min_border=30):
    """Per-label status: 'edge', 'binucleate' or 'ok', plus binucleate partners."""
    med = {l: np.median(dens[lab == l]) for l in range(1, lab.max() + 1)}
    partners = {l: [] for l in med}
    for (p, q), (bmean, blen) in border_stats(lab, dens).items():
        m = min(med[p], med[q])
        if blen >= min_border and (m - bmean) / m < binuc_tau:
            partners[p].append(q); partners[q].append(p)
    status = {l: 'edge' if touches[l] else ('binucleate' if exclude_binuc and partners[l] else 'ok') for l in med}
    return status, partners


def save_outputs(a, g, res, px_um, outdir, prefix, edge_margin=3, binuc_tau=0.0, exclude_binuc=False):
    os.makedirs(outdir, exist_ok=True)
    lab = res['labels']; H, W = lab.shape
    tifffile.imwrite(os.path.join(outdir, f'{prefix}_labels.tif'), lab.astype(np.uint16), imagej=True,
                     resolution=(1 / px_um, 1 / px_um), metadata={'unit': 'micron'})
    rois, rows = [], []
    for r in measure.regionprops(lab):
        c = max(measure.find_contours(np.pad(lab == r.label, 1).astype(float), 0.5), key=len) - 1
        c = measure.approximate_polygon(c, 1.0)
        name = f'cell{r.label:02d}'
        roi = roifile.ImagejRoi.frompoints(np.c_[c[:, 1], c[:, 0]])
        roi.roitype = roifile.ROI_TYPE.POLYGON; roi.name = name
        rois.append(roi)
        minr, minc, maxr, maxc = r.bbox
        touches = minr <= edge_margin or minc <= edge_margin or maxr >= H - edge_margin or maxc >= W - edge_margin
        rows.append(dict(roi=name, area_px=int(r.area), area_um2=round(r.area * px_um ** 2, 1),
                         centroid_x=round(r.centroid[1], 1), centroid_y=round(r.centroid[0], 1),
                         bbox_x=minc, bbox_y=minr, bbox_w=maxc - minc, bbox_h=maxr - minr,
                         touches_border=touches, solidity=round(r.solidity, 3)))
    nuc_mode = res.get('nuclei') is not None
    if nuc_mode:
        status, partners = qc_status(lab, res['dens'], {i + 1: r['touches_border'] for i, r in enumerate(rows)},
                                     binuc_tau, exclude_binuc)
        for i, r in enumerate(rows):
            r['weak_border_with'] = ';'.join(f'cell{q:02d}' for q in partners[i + 1])
            r['status'] = status[i + 1]
    if res.get('nuclei') is not None:
        for r in measure.regionprops(res['nuclei'] * (lab > 0)):
            c = max(measure.find_contours(np.pad(res['nuclei'] == r.label, 1).astype(float), 0.5), key=len) - 1
            c = measure.approximate_polygon(c, 1.0)
            roi = roifile.ImagejRoi.frompoints(np.c_[c[:, 1], c[:, 0]])
            roi.roitype = roifile.ROI_TYPE.POLYGON; roi.name = f'nuc{lab[res["nuclei"] == r.label].max():02d}'
            rois.append(roi)
    keep = {r['roi'] for r in rows if r.get('status', 'ok') == 'ok'}
    for suffix, sel in (('RoiSet', rois), ('RoiSet_filtered', [x for x in rois if x.name.replace('nuc', 'cell') in keep])):
        if suffix == 'RoiSet_filtered' and not nuc_mode:
            continue
        zp = os.path.join(outdir, f'{prefix}_{suffix}.zip')
        if os.path.exists(zp): os.remove(zp)
        if sel: roifile.roiwrite(zp, sel)
    with open(os.path.join(outdir, f'{prefix}_rois.csv'), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    # overlays: (1) saturated-LUT view, (2) original image
    sat = np.clip(g, 0, 6) / 6
    sat = np.dstack([sat * 0.35, sat, sat * 0.35])
    for tag, base in (('overlay_saturated', sat), ('overlay_original', a[..., :3] / 255.0)):
        ov = segmentation.mark_boundaries(base, lab, color=(1, 0.25, 0.25), mode='thick')
        ov = (ov * 255).astype(np.uint8)
        io.imsave(os.path.join(outdir, f'{prefix}_{tag}.png'), ov, check_contrast=False)
    # numbered overlay with matplotlib
    import matplotlib
    from matplotlib.figure import Figure  # no pyplot: safe to call from a worker thread
    fig = Figure(figsize=(18, 6.4)); ax = fig.subplots(1, 3)
    ax[0].imshow(a[..., :3]); ax[0].set_title('original')
    ax[1].imshow(res['dens'], cmap='gray'); ax[1].contour(res['fg'], [0.5], colors='y', linewidths=0.8)
    ax[1].plot(res['seeds'][:, 1], res['seeds'][:, 0], 'c.', ms=6); ax[1].set_title('autofluorescence density + foreground + seeds')
    ax[2].imshow(sat); ax[2].contour(lab > 0, [0.5], colors='none')
    for l, r in enumerate(rows, 1):
        st = r.get('status', 'ok')
        ok = st == 'ok'
        ax[2].contour(lab == l, [0.5], colors=[matplotlib.colormaps['tab20'](l % 20) if not nuc_mode else ('#ff4040' if ok else '#9a9a9a')],
                      linewidths=2.5 if (ok and nuc_mode) else 1.2, linestyles='-' if ok else ':')
        txt = r['roi'].replace('cell', '') + ('' if ok else {'edge': '\nedge', 'binucleate': '\n2 nuclei'}[st])
        ax[2].text(r['centroid_x'], r['centroid_y'], txt, color='w' if ok else '#cccccc', fontsize=11 if ok else 8,
                   ha='center', va='center', weight='bold', bbox=dict(fc='k', alpha=0.5, lw=0))
    if res.get('nuclei') is not None:
        ax[2].contour(res['nuclei'] > 0, [0.5], colors='w', linewidths=0.8, linestyles='--')
        ax[0].imshow(np.dstack([sat[..., 1] * 0.3, sat[..., 1], np.clip(res['nuc_img'] / max(np.percentile(res['nuc_img'], 99.5), 1e-9), 0, 1)]))
        ax[0].set_title('mito (G, saturated) + nuclei (B)')
    n_ok = sum(r.get('status', 'ok') == 'ok' for r in rows)
    ax[2].set_title(f'{len(rows)} cell ROIs, {n_ok} pass QC (saturated LUT)' if nuc_mode else f'{len(rows)} cell ROIs (saturated LUT)')
    for x in ax: x.axis('off')
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f'{prefix}_summary.png'), dpi=110)
    return rows


def load_nuclei(path, k=1.0):
    """Nuclear-stain image (B channel of an RGB TIFF, or grayscale) -> (intensity, labelled nuclei)."""
    na = tifffile.imread(path)
    nuc_img = na[..., 2].astype(float) if na.ndim == 3 else na.astype(float)
    if na.ndim == 3:
        nuc_img[ndi.binary_dilation((na[..., 0] > 200) & (na[..., 1] > 200), iterations=2)] = 0
    return nuc_img, segment_nuclei(nuc_img, k=k)


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('image'); ap.add_argument('-o', '--outdir', default='cell_roi_out')
    ap.add_argument('--nuclei', help='nuclear-stain image of the same field (enables nucleus-seeded mode)')
    ap.add_argument('--sigma', type=float, default=10); ap.add_argument('--fg-k', type=float, default=0.75)
    ap.add_argument('--valley-pct', type=float, default=88); ap.add_argument('--seed-dist', type=int, default=90)
    ap.add_argument('--contrast-tau', type=float, default=0.08); ap.add_argument('--min-area-um2', type=float, default=250)
    ap.add_argument('--binuc-tau', type=float, default=0.0,
                    help='relative border contrast below which two neighbouring ROIs are flagged in weak_border_with')
    ap.add_argument('--exclude-binucleate', action='store_true',
                    help='drop weak-border neighbours as one binucleate cell instead of splitting them')
    A = ap.parse_args()
    a, g, px = load_green(A.image)
    nuc = nuc_img = None
    if A.nuclei:
        nuc_img, nuc = load_nuclei(A.nuclei)
    res = segment(g, A.sigma, A.fg_k, A.valley_pct, A.seed_dist, A.contrast_tau, int(A.min_area_um2 / px ** 2), nuc)
    res['nuc_img'] = nuc_img
    prefix = os.path.splitext(os.path.basename(A.image))[0]
    rows = save_outputs(a, g, res, px, A.outdir, prefix, binuc_tau=A.binuc_tau, exclude_binuc=A.exclude_binucleate)
    print(f'pixel size {px:.4f} um; seeds {len(res["seeds"])}; raw regions {res["raw"].max()}; final ROIs {len(rows)}')
    for r in rows: print(r)
