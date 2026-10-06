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
from skimage import filters, measure, morphology

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
