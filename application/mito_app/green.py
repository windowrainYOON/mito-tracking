"""Protein-of-interest (green) signal on mitochondria.

Two levels:
  cell : green on vs off the MiNA mito mask, enrichment, colocalisation (Pearson, Manders),
         green puncta counts -> correlate with whole-cell MiNA metrics.
  mito : every connected piece of the MiNA mito mask is one mito object; its shape
         (skeleton length, branches, aspect ratio, ...) and the green on it, including the
         green puncta that overlap it -> correlate green load with mito morphology.

Green puncta ("green ROIs"): white top-hat (removes the diffuse cytoplasmic level), then one
threshold for the whole image = max(Otsu, median + 6 MAD) of top-hat values inside the analysed
cells, times a user sensitivity factor. One threshold per image keeps dim and bright cells
comparable within an image.
Intensities are background-subtracted (dark gap level between cells).
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, measure, morphology

from .mina import analyze_skeleton

NAN = float('nan')


def detect_puncta(green, cell_mask, k=1.0, sensitivity=1.0):
    """Label green puncta inside `cell_mask`. Returns (labels, threshold, tophat image)."""
    sm = ndi.gaussian_filter(green, 1.0 * k)
    th = morphology.white_tophat(sm, morphology.disk(max(2, round(8 * k))))
    v = th[cell_mask]
    if v.size == 0 or v.max() <= v.min():
        return np.zeros(green.shape, int), NAN, th
    med = np.median(v)
    t = max(filters.threshold_otsu(v), med + 6 * 1.4826 * np.median(np.abs(v - med))) * sensitivity
    m = morphology.remove_small_objects((th > t) & cell_mask, max_size=max(1, round(3 * k * k)))
    return measure.label(m, connectivity=2), float(t), th


def _safe_div(a, b):
    return a / b if b else NAN


def quantify(green, red, lab, mito, rows, px, bg=0.0, k=1.0, sensitivity=1.0, min_mito_area_um2=0.05, log=print):
    """Add green-on-mito columns to `rows` (in place) and return (mito_rows, puncta_rows, mito_labels, puncta_labels)."""
    g = np.clip(green.astype(float) - bg, 0, None)
    r = red.astype(float)
    cell_ids = [int(row['cell'][4:]) for row in rows]
    cells = np.isin(lab, cell_ids)
    plab, thr, _ = detect_puncta(g, cells, k, sensitivity)
    log(f'green puncta threshold {thr:.1f} (top-hat units), {plab.max()} puncta in analysed cells')
    a_px = px * px

    # mito objects: connected pieces of the mito mask, never crossing a cell boundary;
    # specks below min_mito_area_um2 stay in MiNA but are not mito objects here
    mlab = np.zeros(lab.shape, int)
    nxt = 0
    min_px = max(1, round(min_mito_area_um2 / a_px))
    for cid in cell_ids:
        l = measure.label(morphology.remove_small_objects(mito & (lab == cid), max_size=min_px - 1, connectivity=2),
                          connectivity=2)
        mlab[l > 0] = l[l > 0] + nxt
        nxt += l.max()

    # punctum -> mito object with the largest overlap
    pm_overlap = {}
    both = (plab > 0) & (mlab > 0)
    for p, m in zip(plab[both], mlab[both]):
        pm_overlap.setdefault(p, {}).setdefault(m, 0)
        pm_overlap[p][m] += 1
    dist_to_mito = ndi.distance_transform_edt(~mito) * px

    puncta_rows = []
    p_props = {pr.label: pr for pr in measure.regionprops(plab, intensity_image=g)}
    for p, pr in p_props.items():
        cid = int(np.bincount(lab[tuple(pr.coords.T)]).argmax())
        ov = pm_overlap.get(p, {})
        best = max(ov, key=ov.get) if ov else 0
        vals = g[tuple(pr.coords.T)]
        puncta_rows.append(dict(
            punctum=f'p{p:04d}', cell=f'cell{cid:02d}', area_um2=pr.area * a_px,
            green_mean=float(vals.mean()), green_max=float(vals.max()), green_integrated=float(vals.sum() * a_px),
            centroid_x=pr.centroid[1], centroid_y=pr.centroid[0],
            on_mito=bool(ov), mito=f'm{best:04d}' if best else '',
            overlap_fraction=sum(ov.values()) / pr.area,
            distance_to_mito_um=float(dist_to_mito[tuple(pr.coords.T)].min())))

    # per mito object
    mito_rows = []
    for mr in measure.regionprops(mlab, intensity_image=g):
        sl = mr.slice
        obj = mr.image
        cid = int(lab[sl][obj][0])
        sk = morphology.skeletonize(obj, method='lee').astype(bool)
        edges, ng, br, summed, donuts, end, junc = analyze_skeleton(np.pad(sk, 1), px)
        vals = g[sl][obj]
        on_p = plab[sl][obj]
        p_full = plab[sl] > 0
        p_ids = [p for p, ov in pm_overlap.items() if mr.label in ov]
        major, minor = max(mr.axis_major_length, 1.0), max(mr.axis_minor_length, 1.0)
        mito_rows.append(dict(
            mito=f'm{mr.label:04d}', cell=f'cell{cid:02d}',
            area_um2=mr.area * a_px, length_um=float(sum(e[3] for e in edges)),
            n_branches=len(edges), n_junctions=int(ndi.label(junc, structure=np.ones((3, 3)))[1]),
            n_endpoints=int(end.sum()), donut=bool(donuts),
            major_axis_um=major * px, minor_axis_um=minor * px,
            aspect_ratio=major / minor, solidity=mr.solidity,
            red_mean=float(r[sl][obj].mean()),
            green_mean=float(vals.mean()), green_integrated=float(vals.sum() * a_px),
            green_max=float(vals.max()),
            n_green_puncta=len(p_ids),
            puncta_overlap_area_um2=float((on_p > 0).sum() * a_px),
            puncta_coverage=float((on_p > 0).mean()),
            puncta_green_integrated=float(g[sl][obj & p_full].sum() * a_px),
            puncta_mean_area_um2=float(np.mean([p_props[p].area for p in p_ids]) * a_px) if p_ids else NAN,
            centroid_x=mr.centroid[1], centroid_y=mr.centroid[0]))

    # per cell
    by_cell_m = {}
    for m in mito_rows:
        by_cell_m.setdefault(m['cell'], []).append(m)
    for row in rows:
        cid = int(row['cell'][4:])
        cm = lab == cid
        on, off = cm & mito, cm & ~mito
        gc, rc = g[cm], r[cm]
        p_in = [p for p in puncta_rows if p['cell'] == row['cell']]
        pm = plab[cm] > 0
        ms = by_cell_m.get(row['cell'], [])
        g_on, g_off = (float(g[on].mean()) if on.any() else NAN), (float(g[off].mean()) if off.any() else NAN)
        row.update(
            green_mean_on_mito=g_on, green_mean_off_mito=g_off,
            green_enrichment_on_mito=_safe_div(g_on, g_off),
            green_integrated_on_mito=float(g[on].sum() * a_px),
            green_fraction_on_mito=_safe_div(float(g[on].sum()), float(gc.sum())),
            pearson_red_green=float(np.corrcoef(rc, gc)[0, 1]) if gc.std() > 0 and rc.std() > 0 else NAN,
            manders_m_green=_safe_div(float(g[cm & mito & (plab > 0)].sum()), float(g[cm & (plab > 0)].sum())),
            manders_m_mito=_safe_div(float(r[on & (plab > 0)].sum()), float(r[on].sum())),
            n_green_puncta=len(p_in), n_puncta_on_mito=sum(p['on_mito'] for p in p_in),
            fraction_puncta_on_mito=_safe_div(sum(p['on_mito'] for p in p_in), len(p_in)),
            green_puncta_area_um2=float(pm.sum() * a_px),
            green_puncta_density_per_100um2=_safe_div(len(p_in) * 100, float(cm.sum() * a_px)),
            n_mito_objects=len(ms),
            mito_length_mean_um=float(np.mean([m['length_um'] for m in ms])) if ms else NAN,
            mito_aspect_ratio_mean=float(np.mean([m['aspect_ratio'] for m in ms])) if ms else NAN,
            fraction_mito_with_puncta=_safe_div(sum(m['n_green_puncta'] > 0 for m in ms), len(ms)))
    return mito_rows, puncta_rows, mlab, plab, thr
