"""Individual mitochondria segmentation (adaptive threshold + conservative watershed split).

Alternative to MiNA's single per-cell Otsu mask, for fragmentation analysis:
  1. background: rolling-ball subtraction (radius 1.5 um), Gaussian denoise
  2. foreground: local (adaptive) mean threshold AND per-cell Otsu x fg_floor, like
     Mitochondria Analyzer; small holes filled (sub-resolution holes would become fake MiNA
     loops), specks removed
  3. over-split: intensity watershed from h-maxima markers (one per bright peak)
  4. conservative re-merge: two touching fragments stay apart only if BOTH
       - the intensity dips at the contact: saddle / dimmer peak < saddle_ratio (0.75), and
       - the contact is a neck: contact length < neck_ratio x width of the thinner side
  5. separation: 1-px empty border wherever two objects touch (also diagonally), so MiNA's
     skeleton/graph analysis sees them as separate networks
Size constants are in reference pixels (0.099 um) and scale with the pixel size; the peak
height h scales with the red intensity range so other gains/bit depths behave the same.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, restoration, segmentation

REF_RED_P99 = 141.832  # 99th pct of the preprocessed red inside analysed cells in the reference set

K8 = np.ones((3, 3), bool)
P = dict(rb_radius_um=1.5, sigma=0.7, block_um=1.1, local_offset=0.0, fg_floor=0.5,
         hole_px=6, min_px=8, h=6.0, saddle_ratio=0.75, neck_ratio=1.0, manual_thr=0.0)


def scaled_params(sm, cells, k=1.0, p=None):
    """Scale the reference size (k = REF_PX_UM / px) and peak-height parameters to this image.
    `sm` is the preprocessed red; preprocess with sigma = P['sigma'] * k."""
    q = dict(P if p is None else p)
    q['hole_px'] = max(1, round(q['hole_px'] * k * k))
    q['min_px'] = max(2, round(q['min_px'] * k * k))
    p99 = np.percentile(sm[cells], 99) if cells.any() else REF_RED_P99
    q['h'] *= max(p99, 1e-9) / REF_RED_P99
    return q


def preprocess(red, px_um, p):
    img = red.astype(float)
    bg = restoration.rolling_ball(img, radius=p['rb_radius_um'] / px_um)
    return filters.gaussian(img - bg, p['sigma'], preserve_range=True)


def foreground(sm, cell, px_um, p):
    """Mito pixels of one cell: above the local mean AND above fg_floor x the cell's Otsu level, or, with a
    manual threshold (p['manual_thr'] > 0), simply above that value of the preprocessed red."""
    if p.get('manual_thr', 0) > 0:
        fg = cell & (sm > p['manual_thr'])
    else:
        blk = int(round(p['block_um'] / px_um)) | 1
        local = filters.threshold_local(sm, blk, method='mean', offset=p['local_offset'])
        floor = p['fg_floor'] * filters.threshold_otsu(sm[cell])
        fg = cell & (sm > local) & (sm > floor)
    holes = ndi.binary_fill_holes(fg) & ~fg
    hl = measure.label(holes, connectivity=1)
    small = np.bincount(hl.ravel()) <= p['hole_px']; small[0] = False
    fg |= small[hl]
    fg &= cell
    return morphology.remove_small_objects(fg, max_size=p['min_px'] - 1, connectivity=2)


def split_objects(sm, fg, px_um, p):
    peaks = morphology.h_maxima(np.where(fg, sm, 0), p['h']) & fg
    markers = measure.label(peaks, connectivity=2)
    lab = segmentation.watershed(-sm, markers, mask=fg, connectivity=2)
    # unmarked components (no peak above h) become their own objects
    rest = measure.label(fg & (lab == 0), connectivity=2)
    lab[rest > 0] = rest[rest > 0] + lab.max()
    dist = ndi.distance_transform_edt(fg)
    n = lab.max()
    peak = ndi.maximum(sm, lab, np.arange(n + 1)); peak[0] = 0
    width = 2 * np.asarray(ndi.maximum(dist, lab, np.arange(n + 1)))
    # contacts: per pair, list of (intensity at contact, contact pixels)
    cont = {}
    H, W = lab.shape
    LP = np.pad(lab, 1); SP = np.pad(sm, 1)
    for dy, dx in ((0, 1), (1, 0), (1, 1), (1, -1)):
        a, sa = lab, sm
        b = LP[1 + dy:1 + dy + H, 1 + dx:1 + dx + W]; sb = SP[1 + dy:1 + dy + H, 1 + dx:1 + dx + W]
        m = (a != b) & (a > 0) & (b > 0)
        for i, j, v in zip(a[m].tolist(), b[m].tolist(), np.minimum(sa[m], sb[m]).tolist()):
            k = (min(i, j), max(i, j))
            c = cont.setdefault(k, [0.0, 0])
            c[0] = max(c[0], v); c[1] += 1         # saddle = brightest contact point
    # width of region near contact: use region width (2 x max inscribed radius)
    parent = np.arange(n + 1)

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x

    def keep_cut(sad, clen, a, b):
        dim = min(peak[a], peak[b])
        narrow = clen < p['neck_ratio'] * min(width[a], width[b]) + 1e-9
        return (sad / dim < p['saddle_ratio']) and narrow

    # greedy merge, strongest (highest saddle ratio) first; recompute on merged regions
    changed = True
    while changed:
        changed = False
        merged = {}
        for (i, j), (sad, clen) in cont.items():
            a, b = find(i), find(j)
            if a == b:
                continue
            k = (min(a, b), max(a, b))
            m = merged.setdefault(k, [0.0, 0]); m[0] = max(m[0], sad); m[1] += clen
        order = sorted(merged.items(), key=lambda kv: -kv[1][0] / min(peak[kv[0][0]], peak[kv[0][1]]))
        for (a, b), (sad, clen) in order:
            a, b = find(a), find(b)
            if a == b:
                continue
            if not keep_cut(sad, clen, a, b):
                parent[b] = a
                peak[a] = max(peak[a], peak[b]); width[a] = max(width[a], width[b])
                changed = True
                break  # contacts of the merged region changed; recompute
    roots = np.array([find(x) for x in range(n + 1)])
    lab = roots[lab]
    lab, _, _ = segmentation.relabel_sequential(lab)
    return lab


def separate(lab):
    """Clear pixels of the higher label wherever two labels touch (8-conn): 1-px gaps."""
    mx = ndi.maximum_filter(lab, footprint=K8, mode='constant')
    big = np.where(lab > 0, lab, lab.max() + 1)
    mn = ndi.minimum_filter(big, footprint=K8, mode='constant', cval=lab.max() + 1)
    touch = (lab > 0) & (mn < lab)  # a lower nonzero label is an 8-neighbour
    out = lab.copy(); out[touch] = 0
    del mx
    return out


def object_metrics(lab, px):
    rows = []
    for r in measure.regionprops(lab):
        a = r.area * px * px
        per = r.perimeter * px
        sk = morphology.skeletonize(r.image, method='lee')
        rows.append(dict(label=r.label, area_um2=a, length_um=max(sk.sum(), 1) * px,
                         aspect_ratio=r.axis_major_length / max(r.axis_minor_length, 1),
                         form_factor=per * per / (4 * np.pi * a) if a > 0 else np.nan))
    return rows


def segment_cell(sm, cell, px, p):
    """Return (connected foreground, split labels, separated labels) for one cell mask."""
    fg = foreground(sm, cell, px, p)
    lab = split_objects(sm, fg, px, p)
    return fg, lab, separate(lab)


def fragmentation_metrics(fg, sep, cell, px):
    """Per-cell fragmentation summary of the separated objects."""
    a_px = px * px
    rp = measure.regionprops(sep)
    area = np.array([r.area * a_px for r in rp]) if rp else np.zeros(0)
    ar = np.array([r.axis_major_length / max(r.axis_minor_length, 1) for r in rp]) if rp else np.zeros(0)
    ff = np.array([(r.perimeter * px) ** 2 / (4 * np.pi * r.area * a_px) for r in rp]) if rp else np.zeros(0)
    fpa = fg.sum() * a_px
    small_round = (area < 0.5) & (ar < 2)
    nan = float('nan')
    return dict(
        n_objects_connected=int(measure.label(fg, connectivity=2).max()), n_objects_split=len(rp),
        objects_per_100um2_footprint=100 * len(rp) / fpa if fpa else nan,
        object_area_mean_um2=float(area.mean()) if area.size else nan,
        object_area_median_um2=float(np.median(area)) if area.size else nan,
        area_weighted_mean_object_um2=float((area ** 2).sum() / area.sum()) if area.sum() else nan,
        form_factor_mean=float(np.nanmean(ff)) if ff.size else nan,
        frac_footprint_small_round=float(area[small_round].sum() / area.sum()) if area.sum() else nan,
        frac_objects_small_round=float(small_round.mean()) if area.size else nan)
