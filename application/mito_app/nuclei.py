"""Nucleus detection on a messy nuclear stain (diffuse cytoplasmic blue, shot noise, chromatin holes).

Real nuclei have chromatin texture; the diffuse blue in cell bodies is flat. Steps:
  1. smooth (sigma 2 px), subtract a coarse background (grey opening, 25 um), threshold at 0.5 x Otsu
  2. split the mask by a distance-transform watershed
  3. texture filter per piece: noise-corrected coefficient of variation of the sigma-1.5 image;
     keep a piece if tex >= 0.22, or tex >= 0.17 and it is as bright as the image's textured nuclei
     (median >= 0.7 x ref, ref = median of pieces with tex >= 0.25)
  4. merge kept neighbouring pieces unless their contact is both dark (dip) and narrow (neck),
     so curved nuclei stay whole and truly touching nuclei stay apart
  5. refine each nucleus with its own threshold (halfway between core and surrounding ring),
     close rim bays inside the convex hull, drop non-compact regions (solidity < 0.8)
  6. re-apply the texture rule; in-focus = median >= 0.5 x ref and area >= the minimum nucleus area,
     dim (out of focus) = the rest, if at least max(4000 px, half the median in-focus area) at reference scale
`k` is the pixel-size scale factor (reference pixel size 0.099 um / pixel size); all pixel constants scale with it.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation
from skimage.transform import resize

REF_PX_UM = 0.0990  # reference pixel size the constants are tuned on (cell_roi.REF_PX_UM)
MIN_AREA_PX = 3000   # smallest in-focus nucleus at the reference pixel size
TEX_HI, TEX_LO, TEX_REF = 0.22, 0.17, 0.25
BRIGHT_REL = 0.7     # tex_lo pieces are kept only if this bright relative to the textured nuclei
DIM_REL = 0.5        # below this x ref a nucleus counts as dim
DIM_MIN_PX = 4000    # smallest dim nucleus at the reference pixel size
DIP_REL, NECK_REL = 0.8, 0.6


def flatten(s, k=1.0, um=25.0, px=0.099):
    """Subtract a coarse background (grey opening on a 4x downsampled image)."""
    f = 4
    r = max(2, int(round(um / 2 / px * k / f)))
    bg = ndi.grey_opening(s[::f, ::f], footprint=morphology.disk(r))
    bg = resize(ndi.gaussian_filter(bg, r / 2), s.shape, order=1, preserve_range=True)
    return np.clip(s - bg, 0, None)


def texture(b, lab, k=1.0):
    """{label: (texture, median)}: chromatin texture = noise-corrected SD / median of the sigma-1.5 image
    inside the eroded region; the noise level comes from the high-pass image outside all regions."""
    b = b.astype(float)
    g = ndi.gaussian_filter(b, 1.5 * k)
    hp = b - ndi.gaussian_filter(b, 1.0 * k)
    bgv = hp[lab == 0]
    noise = 1.4826 * np.median(np.abs(bgv - np.median(bgv))) if bgv.size else 0
    noise_s = noise / (2 * np.sqrt(np.pi) * 1.5 * k)
    out = {}
    for r in measure.regionprops(lab):
        m = lab == r.label
        er = ndi.binary_erosion(m, iterations=max(1, int(10 * k)))
        if er.sum() < 100:
            er = m
        v = g[er]; med = np.median(v)
        out[r.label] = (float(np.sqrt(max(v.var() - noise_s ** 2, 0)) / max(med, 1e-9)), float(med))
    return out


def keep_rule(info, ref):
    return [l for l, (t, med) in info.items() if t >= TEX_HI or (t >= TEX_LO and med >= BRIGHT_REL * ref)]


def _split(m, k, px):
    d = filters.gaussian(ndi.distance_transform_edt(m), 2 * k)
    mk = measure.label(morphology.h_maxima(d, max(1.0, 1.0 / px * k * 0.3)))
    return segmentation.watershed(-d, mk, mask=m)


def _merge(s, lab):
    """Merge adjacent pieces unless their contact is dark (boundary median < DIP_REL x the smaller piece's
    median) and narrow (contact length < NECK_REL x the smaller piece's minor axis)."""
    lab = lab.copy()
    st = np.ones((3, 3), bool)
    apart = set()
    changed = True
    while changed:
        changed = False
        props = {r.label: r for r in measure.regionprops(lab)}
        for a, ra in props.items():
            y0, x0, y1, x1 = ra.bbox
            sl = (slice(max(0, y0 - 2), y1 + 2), slice(max(0, x0 - 2), x1 + 2))
            ma = lab[sl] == a
            for b_ in sorted(set(np.unique(lab[sl][ndi.binary_dilation(ma, st) & ~ma])) - {0, a}):
                if b_ < a or (a, b_) in apart:
                    continue
                rb = props[b_]
                sl2 = (slice(max(0, min(ra.bbox[0], rb.bbox[0]) - 2), max(ra.bbox[2], rb.bbox[2]) + 2),
                       slice(max(0, min(ra.bbox[1], rb.bbox[1]) - 2), max(ra.bbox[3], rb.bbox[3]) + 2))
                A, B = lab[sl2] == a, lab[sl2] == b_
                band = (A & ndi.binary_dilation(B, st)) | (B & ndi.binary_dilation(A, st))
                small, rs = (A, ra) if ra.area <= rb.area else (B, rb)
                dip = np.median(s[sl2][band]) < DIP_REL * np.median(s[sl2][small])
                neck = band.sum() / 2.0 < NECK_REL * rs.minor_axis_length
                if dip and neck:
                    apart.add((a, b_))
                    continue
                lab[lab == b_] = a
                changed = True
                break
            if changed:
                break
    return segmentation.relabel_sequential(lab)[0]


def _refine(s, lab, k, band_px=15):
    """Per-nucleus threshold halfway between its core and the surrounding ring; close rim bays in the hull."""
    band = max(2, int(round(band_px * k)))
    out = np.zeros_like(lab)
    for r in measure.regionprops(lab):
        y0, x0, y1, x1 = r.bbox
        sl = (slice(max(0, y0 - 2 * band), y1 + 2 * band), slice(max(0, x0 - 2 * band), x1 + 2 * band))
        m0 = lab[sl] == r.label
        near = ndi.binary_dilation(m0, iterations=band) & ((lab[sl] == 0) | m0)
        ring = (ndi.binary_dilation(m0, iterations=2 * band) & ~ndi.binary_dilation(m0, iterations=band // 2)
                & (lab[sl] == 0))
        core = ndi.binary_erosion(m0, iterations=band // 2)
        m = m0
        if ring.sum() > 20 and core.sum() > 20:
            t = 0.5 * (np.median(s[sl][core]) + np.median(s[sl][ring]))
            mm = ndi.binary_fill_holes((s[sl] > t) & near)
            mm = morphology.opening(mm, morphology.disk(max(1, round(3 * k))))
            cc = measure.label(mm)
            hit = np.unique(cc[m0 & (cc > 0)])
            mm = np.isin(cc, hit[hit > 0])
            if 0.6 * m0.sum() < mm.sum() < 1.5 * m0.sum():
                m = mm
        hull = morphology.convex_hull_image(m) if m.any() else m
        m = ndi.binary_fill_holes(morphology.closing(m, morphology.disk(max(1, round(10 * k)))) & hull)
        out[sl][m & (out[sl] == 0)] = r.label
    return out


def _compact(lab, min_area, min_solidity=0.8):
    out = np.zeros_like(lab); n = 0
    for r in measure.regionprops(lab):
        if r.area >= min_area and r.solidity >= min_solidity:
            n += 1; out[tuple(r.coords.T)] = n
    return out


METHODS = (('texture_merge', 'Texture + merge (recommended)'), ('texture', 'Texture, no merge'),
           ('otsu', 'Global Otsu (old)'))


def detect(b, k=1.0, method='texture_merge'):
    """Nuclear-stain image -> (in-focus labels, dim labels numbered after the in-focus ones, info).

    'texture_merge': steps 1-6 above. 'texture': without the merge step (4) and without the dim size floor,
    so curved nuclei can come out in pieces. 'otsu': global Otsu + compact faint regions as dim nuclei."""
    b = b.astype(float)
    px = REF_PX_UM / k
    if method == 'otsu':
        from . import cell_roi
        nuc = cell_roi.segment_nuclei(b, k=k)
        z = np.zeros(b.shape, int)
        return nuc, cell_roi.find_dim_nuclei(b, nuc, px, k), dict(method=method, ref=0.0, regions=z, texture={},
                                                                     rejected=z, n_rejected=0)
    min_area = MIN_AREA_PX * k * k
    s = flatten(filters.gaussian(b, 2 * k, preserve_range=True), k, px=px)
    if not s.any():
        z = np.zeros(b.shape, int)
        return z, z, dict(method=method, ref=0.0, regions=z, texture={}, rejected=z, n_rejected=0)
    m = s > 0.5 * filters.threshold_otsu(s)
    m = ndi.binary_fill_holes(morphology.opening(m, morphology.disk(max(1, round(4 * k)))))
    m = morphology.remove_small_objects(m, max_size=int(0.4 * min_area))
    lab = _split(m, k, px)
    info = texture(b, lab, k)
    hi = [med for t, med in info.values() if t >= TEX_REF]
    ref = float(np.median(hi)) if hi else float(np.median([med for _, med in info.values()] or [0.0]))
    if method == 'texture_merge':
        lab = segmentation.relabel_sequential(np.where(np.isin(lab, keep_rule(info, ref)), lab, 0))[0]
        lab = _merge(s, lab)
    lab = _compact(_refine(s, lab, k), 0.4 * min_area)
    info = texture(b, lab, k)
    keep = keep_rule(info, ref)
    area = {r.label: r.area for r in measure.regionprops(lab)}
    foc = [l for l in keep if info[l][1] >= DIM_REL * ref and area[l] >= min_area]
    dmin = 0
    if method == 'texture_merge':
        fa = np.median([area[l] for l in foc]) if foc else 0
        dmin = max(DIM_MIN_PX * k * k, 0.5 * fa)
    dim = [l for l in keep if l not in foc and area[l] >= dmin]
    nuc = segmentation.relabel_sequential(np.where(np.isin(lab, foc), lab, 0))[0].astype(int)
    dl = segmentation.relabel_sequential(np.where(np.isin(lab, dim), lab, 0))[0].astype(int)
    dl = np.where(dl > 0, dl + nuc.max(), 0)
    rejected = np.where(np.isin(lab, [l for l in info if l not in foc and l not in dim]), lab, 0)
    return nuc, dl, dict(method=method, ref=ref, regions=lab, texture=info, rejected=rejected,
                         n_rejected=len(info) - len(foc) - len(dim))
