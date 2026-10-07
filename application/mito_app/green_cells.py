"""Green-positive / green-negative cell calls. The cell ROIs are never changed here: they come from
the nucleus + mitochondria segmentation (cell_roi.segment_morph).

  1. bright green mask: green (auto-levelled, sigma 1.5 px-ref smoothed) above a threshold found in
     two Otsu steps (cytoplasm vs. signal, then dim vs. bright signal), never below 3x the
     cytoplasm level, so an image with no expressing cell gives (almost) no mask
  2. a cell's green score = fraction of its cytoplasm (cell minus nuclei) covered by the mask;
     score >= min_fraction -> green-positive
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import filters, morphology

from .cell_roi import REF_GREEN_UNIT


def bright_image(gn, k=1.0):
    """Smoothed auto-levelled green, the image the bright-green threshold is applied to."""
    return ndi.gaussian_filter(gn, 1.5 * k)


def auto_threshold(sm):
    """Automatic bright-green threshold (two Otsu steps on log intensity, never below 3x the cytoplasm level)."""
    v = np.log1p(np.clip(sm, 0, None))
    floor = 3 * REF_GREEN_UNIT
    try:
        t_sig = filters.threshold_multiotsu(v, 3)[1]
        hi = v[v > t_sig]
        t = float(np.expm1(filters.threshold_otsu(hi))) if hi.size > 100 and np.ptp(hi) > 0 else floor
    except ValueError:
        t = floor
    return max(t, floor)


def bright_mask(gn, k=1.0, thr=0.0):
    """Bright (expressing) green pixels of the auto-levelled green image (`thr` > 0: manual threshold).
    Returns (mask, threshold)."""
    sm = bright_image(gn, k)
    t = thr if thr > 0 else auto_threshold(sm)
    m = morphology.remove_small_objects(sm > t, max_size=max(1, round(4 * k * k)))
    return m, t


def scores(lab, G, nuc_mask):
    """Fraction of each cell's cytoplasm covered by the bright mask, and the cytoplasm areas."""
    n = int(lab.max())
    cyto = (lab > 0) & ~nuc_mask
    area = np.bincount(lab[cyto], minlength=n + 1).astype(float)
    hit = np.bincount(lab[cyto & G], minlength=n + 1).astype(float)
    return np.divide(hit, area, out=np.zeros_like(hit), where=area > 0), area


def classify(gn, lab, nuc, k=1.0, min_fraction=0.02, thr=0.0):
    """Returns ({label: 'positive'|'negative'}, info dict with threshold, mask and per-cell fraction).
    `thr` > 0: manual bright-green threshold (auto-levelled units) instead of the automatic one."""
    G, thr = bright_mask(gn, k, thr)
    frac, _ = scores(lab, G, nuc > 0)
    status = {c: 'positive' if frac[c] >= min_fraction else 'negative' for c in range(1, lab.max() + 1)}
    return status, dict(threshold=thr, mask=G, fraction={c: float(frac[c]) for c in status})
