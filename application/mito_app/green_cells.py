"""Green-guided cell ROI refinement and green-positive / green-negative cell calls.

A protein of interest is usually expressed by a whole cell or not at all. So a bright green
patch that sits just outside a cell ROI, or at the edge of a non-expressing neighbour, most
likely belongs to the nearby expressing cell, and its ROI is extended to cover it.

  1. bright green mask: green (auto-levelled, sigma 1.5 px-ref smoothed) above a threshold found in
     two Otsu steps (cytoplasm vs. signal, then dim vs. bright signal), never below 3x the
     cytoplasm level, so an image with no expressing cell gives (almost) no mask
  2. patches = bright pixels grouped within 1 um of each other
  3. a cell's green score = fraction of its cytoplasm (cell minus nuclei) covered by the mask;
     score >= min_fraction -> green-positive
  4. a patch outside all cells, or in a cell that is green-negative without it, is moved to the
     green-positive cell within `near_um` (the highest-scoring one): its envelope is joined to that
     cell with a closing, never taking nuclei or crossing a third cell, and limited to
     `max_gain` of the receiving cell's area
  5. cells are re-scored on the refined ROIs
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation

from .cell_roi import REF_GREEN_UNIT


def bright_mask(gn, k=1.0):
    """Bright (expressing) green pixels of the auto-levelled green image. Returns (mask, threshold)."""
    sm = ndi.gaussian_filter(gn, 1.5 * k)
    v = np.log1p(np.clip(sm, 0, None))
    floor = 3 * REF_GREEN_UNIT
    try:
        t_sig = filters.threshold_multiotsu(v, 3)[1]
        hi = v[v > t_sig]
        t = float(np.expm1(filters.threshold_otsu(hi))) if hi.size > 100 and np.ptp(hi) > 0 else floor
    except ValueError:
        t = floor
    t = max(t, floor)
    m = morphology.remove_small_objects(sm > t, max_size=max(1, round(4 * k * k)))
    return m, t


def scores(lab, G, nuc_mask):
    """Fraction of each cell's cytoplasm covered by the bright mask, and the cytoplasm areas."""
    n = int(lab.max())
    cyto = (lab > 0) & ~nuc_mask
    area = np.bincount(lab[cyto], minlength=n + 1).astype(float)
    hit = np.bincount(lab[cyto & G], minlength=n + 1).astype(float)
    return np.divide(hit, area, out=np.zeros_like(hit), where=area > 0), area


def refine(gn, lab, nuc, px, k=1.0, min_fraction=0.02, near_um=2.0, max_gain=0.25, do_refine=True, log=print):
    """Returns (refined labels, {label: 'positive'|'negative'}, info dict)."""
    G, thr = bright_mask(gn, k)
    nuc_mask = nuc > 0
    lab = lab.copy()
    r_grp = max(1, round(1.0 / px))
    r_near = max(1, round(near_um / px))
    grp = measure.label(morphology.binary_dilation(G, morphology.disk(r_grp)))
    moved = []
    if do_refine and G.any():
        env_all = ndi.binary_fill_holes(morphology.binary_closing(G, morphology.disk(r_grp))) | G
        # biggest patches first, so a cell that is 'positive' only through a neighbour's patch loses it first
        for p in sorted(measure.regionprops(grp), key=lambda r: -G[r.slice][r.image].sum()):
            sl = tuple(slice(max(0, s.start - r_near - 2), s.stop + r_near + 2) for s in p.slice)
            comp = grp[sl] == p.label
            g_p = G[sl] & comp
            if g_p.sum() * px * px < 0.25:
                continue
            frac, area = scores(lab, G, nuc_mask)
            L = lab[sl]
            env = env_all[sl] & comp & ~nuc_mask[sl]
            ov = np.bincount(L[g_p], minlength=lab.max() + 1)
            owner = int(ov[1:].argmax() + 1) if ov[1:].any() and ov[1:].max() >= ov[0] else 0
            # owner's score without this patch
            if owner:
                own_frac = (frac[owner] * area[owner] - ((L == owner) & g_p & ~nuc_mask[sl]).sum()) / max(area[owner], 1)
                if own_frac >= min_fraction:
                    continue
            near = ndi.binary_dilation(g_p, morphology.disk(r_near))
            cands = [c for c in np.unique(L[near]) if c and c != owner and frac[c] >= min_fraction]
            if not cands:
                continue
            b = max(cands, key=lambda c: frac[c])
            take = env & np.isin(L, [0, owner])
            if take.sum() > max_gain * (lab == b).sum():
                continue
            allowed = np.isin(L, [0, owner, b]) & (~nuc_mask[sl] | (L == b))
            joined = morphology.binary_closing((L == b) | take, morphology.disk(r_near)) & allowed
            cc = measure.label(joined | (L == b))
            keep = np.isin(cc, np.unique(cc[L == b]))
            add = keep & (L != b) & allowed
            if add.sum() > max_gain * (lab == b).sum():
                continue
            gain = int(add.sum())
            L[add] = b
            if owner:  # the owner keeps its largest piece; cut-off bits go to the receiver if touching it
                cc_a = measure.label(L == owner)
                if cc_a.max() > 1:
                    main = np.bincount(cc_a.ravel())[1:].argmax() + 1
                    for piece in range(1, cc_a.max() + 1):
                        if piece == main:
                            continue
                        pm = cc_a == piece
                        L[pm] = b if (ndi.binary_dilation(pm) & (L == b)).any() else 0
            lab[sl] = L
            moved.append(dict(patch_um2=float(g_p.sum() * px * px), from_cell=f'cell{owner:02d}' if owner else 'outside',
                              to_cell=f'cell{b:02d}', added_um2=float(gain * px * px)))
            log(f'green patch {g_p.sum() * px * px:.1f} um2 '
                f'{"outside the cells" if not owner else f"in cell{owner:02d}"} -> cell{b:02d} (+{gain * px * px:.1f} um2)')
    frac, _ = scores(lab, G, nuc_mask)
    status = {c: 'positive' if frac[c] >= min_fraction else 'negative' for c in range(1, lab.max() + 1)}
    return lab, status, dict(threshold=thr, mask=G, fraction={c: float(frac[c]) for c in status}, moved=moved)


def _touching(lab, ids):
    """Labels touching any of `ids` (one-pixel dilation)."""
    m = ndi.binary_dilation(np.isin(lab, list(ids)), iterations=2)
    return set(np.unique(lab[m]).tolist()) - {0}


def _keep_largest(lab, ids):
    """Each label in `ids` keeps its largest piece; other pieces go to the neighbour with the longest contact."""
    for l in ids:
        cc = measure.label(lab == l)
        if cc.max() <= 1:
            continue
        main = np.bincount(cc.ravel())[1:].argmax() + 1
        for piece in range(1, cc.max() + 1):
            if piece == main:
                continue
            pm = cc == piece
            ring = ndi.binary_dilation(pm) & ~pm
            nb = lab[ring]; nb = nb[(nb != l) & (nb > 0)]
            lab[pm] = np.bincount(nb).argmax() if nb.size else 0
    return lab


def split_by_intensity(gn, lab, nuc, nuc_img, landscape, G, px, k, min_fraction, min_area_px, pos=None, log=print):
    """Redraw the ROIs of green-positive cells and their neighbours with green intensity as the main weight.

    Markers: the in-focus nuclei (cells keep their ids), dim (out-of-focus) nuclei, and bright green
    territories that hold no nucleus at all. Watershed over (1.5 x inverted log green + the nucleus-mode
    landscape), so borders between cells follow the dark gaps of the green image. A new cell (from a dim
    nucleus or a green territory) is kept only if it is at least `min_area_px` and itself green-positive;
    otherwise its pixels go back to the cells they came from.
    Returns (labels, {new id: 'dim nucleus' | 'green territory'})."""
    nuc_mask = nuc > 0
    if pos is None:
        frac, _ = scores(lab, G, nuc_mask)
        pos = [c for c in range(1, lab.max() + 1) if frac[c] >= min_fraction]
    if not pos:
        return lab, {}
    U = np.isin(lab, list(_touching(lab, pos) | set(pos)))
    markers = np.where(nuc_mask & U, lab, 0)
    nxt = int(lab.max()) + 1
    new = {}
    # dim nuclei: above 0.4 x the in-focus Otsu level, away from the in-focus nuclei
    s = ndi.gaussian_filter(nuc_img.astype(float), 2 * k)
    t = filters.threshold_otsu(s)
    dim = (s > 0.4 * t) & U & ~ndi.binary_dilation(nuc_mask, iterations=max(1, round(1.5 / px)))
    dim = morphology.binary_opening(dim, morphology.disk(max(1, round(3 * k))))
    dim = morphology.remove_small_objects(dim, max_size=int(800 * k * k))
    for r in measure.regionprops(measure.label(dim)):
        markers[tuple(r.coords.T)] = nxt; new[nxt] = 'dim nucleus'; nxt += 1
    # bright green territories with no nucleus (in focus or dim) inside
    dens = ndi.gaussian_filter(G.astype(float), 2.0 / px)
    # (pieces within 6 um of each other, with no nucleus between them, are one cell)
    terr = morphology.remove_small_objects((dens > 0.08) & U, max_size=max(1, int(0.3 * min_area_px)))
    tl = measure.label(terr)
    groups = measure.label(ndi.binary_dilation(terr, iterations=max(1, round(3.0 / px))) & U)
    seeded = set(np.unique(groups[(markers > 0) & terr]).tolist())
    gid = {}
    for r in measure.regionprops(tl, intensity_image=dens):
        if markers[tuple(r.coords.T)].any():
            continue
        g = int(groups[tuple(r.coords[0])])
        if g in seeded:
            continue
        if g not in gid:
            gid[g] = nxt; new[nxt] = 'green territory'; nxt += 1
        core = r.coords[dens[tuple(r.coords.T)] >= 0.5 * r.intensity_max]
        markers[tuple(core.T)] = gid[g]
    sm = ndi.gaussian_filter(gn, 1.0 / px)
    ref = np.log1p(np.percentile(sm[U], 99)) or 1.0
    inten = np.clip(np.log1p(np.clip(sm, 0, None)) / ref, 0, 1)
    out = lab.copy()
    ws = segmentation.watershed(landscape + 1.5 * (1 - inten), markers, mask=U)
    out[U & (ws > 0)] = ws[U & (ws > 0)]
    # only borders that involve a green-positive (or new) cell move; two green-negative neighbours keep theirs
    keep_old = ~np.isin(out, pos + list(new)) & ~np.isin(lab, pos)
    out[keep_old] = lab[keep_old]
    out = _keep_largest(out, list(_touching(lab, pos) | set(pos)) + list(new))
    # validate the new cells
    kept = {}
    for nid, why in new.items():
        m = out == nid
        if not m.any():
            continue
        f, _ = scores(np.where(m, 1, 0), G, nuc_mask)
        if m.sum() >= min_area_px and f[1] >= min_fraction:
            kept[nid] = why
            continue
        out[m] = lab[m]  # rejected: back to the cells it was taken from
    out = _keep_largest(out, list(set(np.unique(out[U]).tolist()) - {0}))
    # renumber kept new cells after the existing ones
    for i, nid in enumerate(sorted(kept), int(lab.max()) + 1):
        out[out == nid] = i
    kept = {i: kept[nid] for i, nid in enumerate(sorted(kept), int(lab.max()) + 1)}
    for nid, why in kept.items():
        log(f'new green-positive cell{nid:02d} split off by green intensity (seed: {why}, '
            f'{(out == nid).sum() * px * px:.0f} um2)')
    return out, kept


def process(gn, lab, nuc, nuc_img, landscape, px, k=1.0, min_fraction=0.02, min_area_px=0,
            do_split=True, do_refine=True, log=print):
    """Intensity split of green-positive cells, then the green-pattern patch refinement."""
    new = {}
    if do_split:
        G, _ = bright_mask(gn, k)
        # which cells really express: call them after a trial patch refinement, so a neighbour's patch
        # at a green-negative cell's edge does not make its borders count as green-positive ones
        _, st0, _ = refine(gn, lab, nuc, px, k, min_fraction, do_refine=do_refine, log=lambda s: None)
        pos = [c for c, v in st0.items() if v == 'positive']
        lab, new = split_by_intensity(gn, lab, nuc, nuc_img, landscape, G, px, k, min_fraction, min_area_px, pos, log)
    lab, status, info = refine(gn, lab, nuc, px, k, min_fraction, do_refine=do_refine, log=log)
    info['new_cells'] = new
    return lab, status, info
