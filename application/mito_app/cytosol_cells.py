"""Cell ROIs the way Cytosol Analysis (Cytosol Viewer) detects them.

Ported from cytosol-analysis cytosol/cellroi.py (main, b367f01): RoiParams, the nucleus
segmentation and segment_cells are copied unchanged, so the same image and settings give the same ROIs
as Cytosol Viewer's automatic detection. Only the result differs: plain label images instead of the
CellRois editing class, with the frame-edge handling (CellRois.apply_edge_mode) as edge_mode().

Pipeline (segment_cells):
  1. background channel (here the mito channel by default): smooth and saturate, so the dim signal fills
     the whole cytoplasm while the gaps between cells stay dark
  2. nuclei: Otsu on the nucleus channel, or a manual threshold in nucleus-channel units; touching nuclei
     split by a distance-transform watershed
  3. merge: max(saturated background, nuclei) -> each cell is one solid blob; foreground = Otsu on that
     merged image, or a manual threshold given in background-channel units
  4. one nucleus = one cell: watershed from the nuclei over the dark cell-cell valleys (Sato filter)
  5. cells cut by the image frame: 'trim' replaces the frame side by an inner border edge_margin_um inside
     the frame; 'exclude' drops them, 'keep' leaves them
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage as ndi
from skimage import feature, filters, measure, morphology, segmentation

EDGE_MODES = ('trim', 'exclude', 'keep')


@dataclass
class RoiParams:
    saturation: float | None = None  # background value shown as full white; None = auto
    sat_percentile: float = 50.0     # auto saturation = this percentile of the smoothed channel
    smooth_um: float = 0.8           # Gaussian smoothing of the merged image
    fg_k: float = 1.0                # foreground threshold = Otsu * fg_k
    min_cell_um2: float = 80.0
    min_nucleus_um2: float = 15.0
    nucleus_split_um: float = 6.0    # min distance between nucleus centers when splitting
    nucleus_neck_um: float = 1.0     # how much narrower the neck between two touching nuclei must be
    edge_mode: str = 'trim'
    edge_margin_um: float = 1.0
    work_px_um: float = 0.2          # computations run at about this pixel size
    # manual thresholds in raw channel values (None = Otsu): nucleus channel
    # (after light smoothing) and background channel (cell vs. gap)
    nuc_threshold: float | None = None
    cell_threshold: float | None = None


def _norm(img):
    lo, hi = np.percentile(img, [0.5, 99.5])
    return np.clip((img - lo) / max(hi - lo, 1e-9), 0, 1)


def auto_saturation(bg: np.ndarray, px_um: float, p: RoiParams = RoiParams()) -> float:
    s = filters.gaussian(bg.astype(np.float32), max(0.5 / px_um, 0.5), preserve_range=True)
    v = float(np.percentile(s, p.sat_percentile))
    return v if v > 0 else float(s.max() or 1.0)


def brighten(bg: np.ndarray, saturation: float) -> np.ndarray:
    """Background channel with a bright LUT: 0..saturation -> 0..1, clipped."""
    lo = float(np.percentile(bg, 0.5))
    return np.clip((bg.astype(np.float32) - lo) / max(saturation - lo, 1e-9), 0, 1)


def nucleus_signal(nuc: np.ndarray, px_um: float) -> np.ndarray:
    """The nucleus channel as thresholded: lightly smoothed, raw units."""
    return filters.gaussian(nuc.astype(np.float32), max(0.3 / px_um, 0.5), preserve_range=True)


def nucleus_threshold(s: np.ndarray, p: RoiParams) -> float:
    return float(p.nuc_threshold) if p.nuc_threshold is not None else float(filters.threshold_otsu(s))


def segment_nuclei(nuc: np.ndarray, px_um: float, p: RoiParams) -> np.ndarray:
    s = nucleus_signal(nuc, px_um)
    m = s > nucleus_threshold(s, p)
    m = ndi.binary_opening(m, morphology.disk(max(1, round(0.4 / px_um))))
    m = ndi.binary_fill_holes(m)
    m = morphology.remove_small_objects(m, max_size=int(p.min_nucleus_um2 / px_um ** 2))
    # split touching nuclei
    dist = ndi.distance_transform_edt(m)
    dist = filters.gaussian(dist, max(0.5 / px_um, 1))
    lab0 = measure.label(m)
    # a seed per distance maximum that stands out by `nucleus_neck_um` from the
    # neck joining it to its neighbour, so notched nuclei are not cut
    peaks = morphology.h_maxima(dist, p.nucleus_neck_um / px_um) & m
    pk = np.array([r.centroid for r in measure.regionprops(measure.label(peaks))], int).reshape(-1, 2)
    # seeds of one nucleus closer than nucleus_split_um: keep the deepest
    keep = []
    for y, x in sorted(pk.tolist(), key=lambda q: -dist[q[0], q[1]]):
        if all(lab0[y, x] != lab0[a, b] or np.hypot(y - a, x - b) * px_um >= p.nucleus_split_um
               for a, b in keep):
            keep.append((y, x))
    markers = np.zeros(m.shape, np.int32)
    for i, (y, x) in enumerate(keep, 1):
        markers[y, x] = i
    for r in measure.regionprops(lab0):  # every nucleus gets at least one seed
        if not markers[r.slice][r.image].any():
            cy, cx = map(int, r.coords[np.argmax(dist[tuple(r.coords.T)])])
            markers[cy, cx] = markers.max() + 1
    lab = segmentation.watershed(-dist, markers, mask=m)
    # drop fragments the split left too small
    small = [r.label for r in measure.regionprops(lab) if r.area < p.min_nucleus_um2 / px_um ** 2]
    if small:
        lab[np.isin(lab, small)] = 0
    return segmentation.relabel_sequential(lab)[0]


def _border_stats(lab, img):
    """{(a, b): (mean img along the shared border of labels a < b, border length)}."""
    keys, vals = [], []
    for a, b, ia, ib in ((lab[:, :-1], lab[:, 1:], img[:, :-1], img[:, 1:]),
                         (lab[:-1, :], lab[1:, :], img[:-1, :], img[1:, :])):
        m = (a != b) & (a > 0) & (b > 0)
        pa, pb = np.minimum(a[m], b[m]).astype(np.int64), np.maximum(a[m], b[m]).astype(np.int64)
        keys.append(pa * (int(lab.max()) + 1) + pb)
        vals.append((ia[m] + ib[m]) / 2)
    keys, vals = np.concatenate(keys), np.concatenate(vals)
    if not len(keys):
        return {}
    u, inv = np.unique(keys, return_inverse=True)
    n = np.bincount(inv)
    mean = np.bincount(inv, weights=vals) / n
    k = int(lab.max()) + 1
    return {(int(x // k), int(x % k)): (float(v), int(c)) for x, v, c in zip(u, mean, n)}


def _merge_regions(lab, img, contrast_tau, min_area, min_border):
    lab = lab.copy()
    # per-label area and median, updated only for the merged label
    fgl, fgv = lab[lab > 0], img[lab > 0]
    order = np.lexsort((fgv, fgl))
    fgl, fgv = fgl[order], fgv[order]
    ids, start, area = np.unique(fgl, return_index=True, return_counts=True)
    med = {int(i): float(np.median(fgv[s0:s0 + n])) for i, s0, n in zip(ids, start, area)}
    area = dict(zip(ids.tolist(), area.tolist()))
    while len(area) > 1:
        best = None
        for (a, b), (bmean, blen) in _border_stats(lab, img).items():
            if blen < min_border:
                continue
            contrast = min(med[a], med[b]) - bmean
            small = min(area[a], area[b]) < min_area
            if contrast < contrast_tau or small:
                key = (0 if small else 1, contrast)
                if best is None or key < best[0]:
                    best = (key, a, b)
        if best is None:
            break
        _, a, b = best
        lab[lab == b] = a
        area[a] += area.pop(b)
        med.pop(b)
        med[a] = float(np.median(img[lab == a]))
    return segmentation.relabel_sequential(lab)[0]


def _largest(mask):
    cc = measure.label(mask)
    if cc.max() <= 1:
        return mask
    return cc == np.argmax(np.bincount(cc.ravel())[1:]) + 1


def _smooth_labels(lab, radius, keep=None):
    """Open each ROI, fill holes, keep its largest piece; `keep` pixels are never lost."""
    out = np.zeros_like(lab)
    se = morphology.disk(max(1, radius))
    for r in measure.regionprops(lab):
        sl = tuple(slice(max(s.start - radius - 1, 0), s.stop + radius + 1) for s in r.slice)
        m = ndi.binary_opening(lab[sl] == r.label, se)
        if keep is not None:
            m |= keep[sl] == r.label
        m = ndi.binary_fill_holes(m)
        if m.any():
            m = _largest(m)
            o = out[sl]
            o[m & (o == 0)] = r.label
    return out


def _bg_floor(bg) -> float:
    return float(np.percentile(bg, 0.5))  # the 0 of brighten()


def _cell_threshold(dens, bg, sat, p: RoiParams) -> tuple[float, float]:
    """(threshold on the merged 0..1 image, the same in background-channel units)."""
    lo = _bg_floor(bg)
    span = max(sat - lo, 1e-9)
    if p.cell_threshold is not None:
        return (float(p.cell_threshold) - lo) / span, float(p.cell_threshold)
    t01 = float(filters.threshold_otsu(dens)) * p.fg_k
    return t01, lo + t01 * span


def segment_cells(bg: np.ndarray, nuc: np.ndarray | None, px_um: float | None,
                  p: RoiParams | None = None) -> dict:
    """Cell ROIs from a background-signal plane and an optional nucleus plane."""
    p = p or RoiParams()
    px_um = px_um or 0.1
    H, W = bg.shape
    k = max(1, int(round(p.work_px_um / px_um)))
    wpx = px_um * k

    def down(a):
        a = a.astype(np.float32)
        if k == 1:
            return a
        h, w = (H // k) * k, (W // k) * k
        a = a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))
        return a

    sat = p.saturation if p.saturation else auto_saturation(bg, px_um, p)
    b = down(brighten(bg, sat))
    nuclei = None
    nuc_thr = None
    if nuc is not None:
        dn_ = down(nuc)
        nuc_thr = nucleus_threshold(nucleus_signal(dn_, wpx), p)
        nuclei = segment_nuclei(dn_, wpx, p)
        merged = np.maximum(b, (nuclei > 0).astype(np.float32))
    else:
        merged = b
    sig = max(p.smooth_um / wpx, 0.5)
    dens = filters.gaussian(merged, sig)
    t01, cell_thr = _cell_threshold(dens, bg, sat, p)
    fg = dens > t01
    r_open = max(1, round(1.0 / wpx))
    fg = ndi.binary_opening(fg, morphology.disk(r_open))
    min_cell = p.min_cell_um2 / wpx ** 2
    fg = morphology.remove_small_objects(fg, max_size=int(min_cell * 0.3))
    fg = morphology.remove_small_holes(fg, max_size=int(min_cell * 0.3))

    d4 = filters.gaussian(merged, max(0.4 / wpx, 0.5))
    sig_v = [max(0.8 / wpx, 1), max(1.4 / wpx, 1.5)]
    v = filters.sato(d4, sigmas=sig_v, black_ridges=True)
    vn = v / max(np.percentile(v, 99.5), 1e-9)

    if nuclei is not None and nuclei.max() > 0:
        fg |= nuclei > 0
        dn = dens / max(np.percentile(dens[fg], 99), 1e-9)
        lab = segmentation.watershed(vn + 0.5 * (1 - np.clip(dn, 0, 1)), nuclei, mask=fg)
        lab = _smooth_labels(lab, max(1, round(0.6 / wpx)), keep=nuclei)
    else:
        valley = vn > np.percentile(vn[fg], 88) if fg.any() else np.zeros_like(fg)
        dist = filters.gaussian(ndi.distance_transform_edt(fg & ~valley), 1.5)
        pk = feature.peak_local_max(dist, min_distance=max(2, int(4 / wpx)),
                                    threshold_abs=max(1.0, 1.0 / wpx),
                                    labels=measure.label(fg), exclude_border=False)
        markers = np.zeros(fg.shape, np.int32)
        markers[tuple(pk.T)] = np.arange(1, len(pk) + 1)
        lab = segmentation.watershed(vn - 0.5 * dist / max(dist.max(), 1e-9), markers, mask=fg)
        lab = _merge_regions(lab, dens, 0.08, min_cell, max(3, int(1.5 / wpx)))
        lab = _smooth_labels(lab, max(1, round(0.6 / wpx)))
        nuclei = None

    # drop small cells (and their nuclei)
    for r in measure.regionprops(lab):
        if r.area < min_cell:
            lab[lab == r.label] = 0
    lab, fw, _ = segmentation.relabel_sequential(lab)
    if nuclei is not None:
        nuclei = np.asarray(fw[nuclei]) * (lab == np.asarray(fw[nuclei]))

    def up(a):
        if a is None:
            return None
        if k > 1:
            a = np.repeat(np.repeat(a, k, 0), k, 1)
        out = np.zeros((H, W), np.int32)
        h, w = min(H, a.shape[0]), min(W, a.shape[1])
        out[:h, :w] = a[:h, :w]
        if h < H:
            out[h:, :w] = out[h - 1:h, :w]
        if w < W:
            out[:, w:] = out[:, w - 1:w]
        return out

    lab, nuclei = up(lab), up(nuclei)
    if k > 1:  # smooth the blocky upsampled borders
        lab = _smooth_labels(lab, max(1, k // 2 + 1), keep=nuclei)

    lab = lab.astype(np.int32)
    before = lab.copy()
    status, cut = edge_mode(lab, nuclei, px_um, p.edge_mode, p.edge_margin_um, p.min_cell_um2)
    trimmed = np.isin(before, list(cut)) & (lab == 0)
    fw = renumber(lab)
    lab = fw[lab]
    if nuclei is not None:
        nuclei = np.where(lab > 0, fw[nuclei], 0).astype(np.int32)
    status = {int(fw[l]): s for l, s in status.items() if fw[l]}
    cut = {int(fw[l]): n for l, n in cut.items() if fw[l]}
    return dict(labels=lab.astype(np.int32), nuclei=nuclei, status=status, trimmed=trimmed, trimmed_px=cut,
                saturation=sat, nuc_threshold=nuc_thr, cell_threshold=cell_thr)


def renumber(lab) -> np.ndarray:
    """Lookup table old -> new label, 1..n top-to-bottom, left-to-right (CellRois.renumber), so the cells get
    the same numbers as in Cytosol Viewer."""
    ids = [int(i) for i in np.unique(lab) if i > 0]
    cy = {l: np.argwhere(lab == l).mean(0) for l in ids}
    order = sorted(ids, key=lambda l: (round(cy[l][0] / 50), cy[l][1]))
    lut = np.zeros(int(lab.max()) + 1, np.int32)
    for new, old in enumerate(order, 1):
        lut[old] = new
    return lut


def edge_mode(lab, nuclei, px_um, mode='trim', margin_um=1.0, min_cell_um2=0.0) -> dict:
    """Cells cut by the image frame, in place (CellRois.apply_edge_mode).
    Returns ({label: 'ok' | 'edge' | 'edge-trimmed'}, {trimmed label: pixels cut off})."""
    if mode not in EDGE_MODES:
        raise ValueError(f'edge mode must be one of {EDGE_MODES}')
    H, W = lab.shape
    px = px_um or 0.1
    m_px = max(2, int(round(margin_um / px)))
    yy, xx = np.ogrid[:H, :W]
    inner = (yy >= m_px) & (yy < H - m_px) & (xx >= m_px) & (xx < W - m_px)
    frame = ~((yy >= 1) & (yy < H - 1) & (xx >= 1) & (xx < W - 1))
    on_frame = set(np.unique(lab[frame]).tolist()) - {0}
    status, cut = {}, {}

    def erase(l):
        lab[lab == l] = 0
        if nuclei is not None:
            nuclei[nuclei == l] = 0

    for l in [int(i) for i in np.unique(lab) if i > 0]:
        if l not in on_frame:
            status[l] = 'ok'
        elif mode == 'keep':
            status[l] = 'edge'
        elif mode == 'exclude':
            erase(l)
        else:
            cell = (lab == l) & inner
            nuc = nuclei == l if nuclei is not None else None
            cc = measure.label(cell)
            if cc.max() == 0:
                erase(l)
                continue
            # keep the piece holding most of the nucleus, else the largest
            if nuc is not None and (nuc & cell).any():
                keep = np.argmax(np.bincount(cc[nuc & cell], minlength=cc.max() + 1)[1:]) + 1
            else:
                keep = np.argmax(np.bincount(cc.ravel())[1:]) + 1
            new = cc == keep
            if new.sum() * px ** 2 < min_cell_um2:
                erase(l)
                continue
            cut[l] = int(((lab == l) & ~new).sum())
            lab[(lab == l) & ~new] = 0
            if nuclei is not None:
                nuclei[(nuclei == l) & ~new] = 0
            status[l] = 'edge-trimmed'
    return status, cut


def _work_scale(shape, px_um, p):
    k = max(1, int(round(p.work_px_um / px_um)))
    H, W = shape

    def down(a):
        a = a.astype(np.float32)
        if k == 1:
            return a
        h, w = (H // k) * k, (W // k) * k
        return a[:h, :w].reshape(h // k, k, w // k, k).mean((1, 3))

    def up(a):
        if k > 1:
            a = np.repeat(np.repeat(a, k, 0), k, 1)
        out = np.zeros((H, W), a.dtype)
        h, w = min(H, a.shape[0]), min(W, a.shape[1])
        out[:h, :w] = a[:h, :w]
        out[h:, :w] = out[h - 1:h, :w]
        out[:, w:] = out[:, w - 1:w]
        return out
    return k, down, up


def nucleus_image(nuc, px_um, p: RoiParams | None = None):
    """(image the nucleus threshold applies to, its Otsu value), full size: `image > t` is the nucleus mask
    segment_cells starts from (as in Cytosol Viewer's threshold preview)."""
    p = p or RoiParams()
    k, down, up = _work_scale(nuc.shape, px_um, p)
    s = nucleus_signal(down(nuc), px_um * k)
    return up(s), float(filters.threshold_otsu(s))


def cell_image(bg, nuc, px_um, p: RoiParams | None = None):
    """(image the cell threshold applies to, in background-channel units, its auto value), full size:
    `image > t` is the cell foreground of segment_cells before opening and hole filling."""
    p = p or RoiParams()
    k, down, up = _work_scale(bg.shape, px_um, p)
    wpx = px_um * k
    sat = p.saturation if p.saturation else auto_saturation(bg, px_um, p)
    merged = down(brighten(bg, sat))
    if nuc is not None:
        s = nucleus_signal(down(nuc), wpx)
        merged = np.maximum(merged, (s > nucleus_threshold(s, p)).astype(np.float32))
    dens = filters.gaussian(merged, max(p.smooth_um / wpx, 0.5))
    lo = _bg_floor(bg)
    span = max(sat - lo, 1e-9)
    t01 = float(filters.threshold_otsu(dens)) * p.fg_k
    return up(lo + dens * span), lo + t01 * span
