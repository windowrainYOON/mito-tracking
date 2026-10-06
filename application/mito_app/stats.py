"""Metric catalogue and rank correlations between mito morphology and green load."""
import math

import numpy as np
from scipy import stats

# (column, label) pairs offered as X (mito) and Y (green) in the correlation view and table
CELL_X = [('footprint_fraction', 'Mito footprint fraction'), ('mitochondrial_footprint_um2', 'Mito footprint (µm²)'),
          ('n_networks', 'Networks'), ('n_branches', 'Branches'), ('total_branch_length_um', 'Total branch length (µm)'),
          ('branch_length_mean_um', 'Branch length mean (µm)'),
          ('summed_branch_lengths_mean_um', 'Network length mean (µm)'),
          ('network_branches_mean', 'Branches per network'), ('donuts', 'Donuts'),
          ('n_mito_objects', 'Mito objects'), ('mito_length_mean_um', 'Mito object length mean (µm)'),
          ('mito_aspect_ratio_mean', 'Mito aspect ratio mean'),
          ('objects_per_100um2_footprint', 'Objects per 100 µm² footprint'),
          ('area_weighted_mean_object_um2', 'Area-weighted object size (µm²)'),
          ('form_factor_mean', 'Form factor mean'), ('frac_footprint_small_round', 'Footprint in small round objects')]
CELL_Y = [('green_poi_mean_cell', 'Green mean, whole cell'), ('green_mean_on_mito', 'Green mean on mito'),
          ('green_mean_off_mito', 'Green mean off mito'), ('green_enrichment_on_mito', 'Green enrichment on/off mito'),
          ('green_fraction_on_mito', 'Fraction of green on mito'), ('green_integrated_on_mito', 'Green integrated on mito'),
          ('pearson_red_green', 'Pearson red–green'), ('manders_m_green', 'Manders: puncta green on mito'),
          ('manders_m_mito', 'Manders: mito red under puncta'), ('n_green_puncta', 'Green puncta'),
          ('n_puncta_on_mito', 'Green puncta on mito'), ('fraction_puncta_on_mito', 'Fraction of puncta on mito'),
          ('green_puncta_density_per_100um2', 'Puncta per 100 µm²'),
          ('fraction_mito_with_puncta', 'Fraction of mito objects with puncta'),
          ('green_bright_area_percent', 'Bright green area (% cytoplasm)')]
MITO_X = [('length_um', 'Length (µm)'), ('aspect_ratio', 'Aspect ratio'), ('area_um2', 'Area (µm²)'),
          ('n_branches', 'Branches'), ('n_junctions', 'Junctions'), ('n_endpoints', 'End points'),
          ('solidity', 'Solidity'), ('form_factor', 'Form factor'), ('major_axis_um', 'Major axis (µm)'), ('red_mean', 'Red mean')]
MITO_Y = [('green_mean', 'Green mean'), ('green_integrated', 'Green integrated'), ('green_max', 'Green max'),
          ('n_green_puncta', 'Green puncta'), ('puncta_coverage', 'Puncta coverage'),
          ('puncta_overlap_area_um2', 'Puncta overlap area (µm²)'),
          ('puncta_green_integrated', 'Puncta green integrated'), ('puncta_mean_area_um2', 'Puncta mean size (µm²)')]


def available(metrics, rows):
    """Metrics that are present in the data (fragmentation metrics exist only for split objects)."""
    keys = set(rows[0]) if rows else set()
    return [m for m in metrics if m[0] in keys]


def numeric(rows, key):
    out = []
    for r in rows:
        try:
            v = float(r.get(key, 'nan'))
        except (TypeError, ValueError):
            v = math.nan
        out.append(v)
    return np.array(out)


def spearman(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return math.nan, math.nan, n
    res = stats.spearmanr(x[ok], y[ok])
    return float(res.statistic), float(res.pvalue), n


def regression(x, y):
    """Least-squares line y = slope * x + intercept with Pearson r, R² and the two-sided p of the slope."""
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    nan = dict(slope=math.nan, intercept=math.nan, r=math.nan, r2=math.nan, p=math.nan, n=n)
    if n < 3 or np.ptp(x[ok]) == 0 or np.ptp(y[ok]) == 0:
        return nan
    res = stats.linregress(x[ok], y[ok])
    return dict(slope=float(res.slope), intercept=float(res.intercept), r=float(res.rvalue),
                r2=float(res.rvalue ** 2), p=float(res.pvalue), n=n)


def pair_stats(x, y):
    """Pearson regression and Spearman rank correlation of one metric pair."""
    reg = regression(x, y)
    rho, p_s, _ = spearman(x, y)
    return dict(pearson_r=reg['r'], r2=reg['r2'], slope=reg['slope'], intercept=reg['intercept'], p_value=reg['p'],
                spearman_rho=rho, spearman_p=p_s, n=reg['n'])


def correlation_table(cell_rows, mito_rows):
    """Pearson r / linear regression and Spearman rho for every mito metric x green metric pair, at cell and mito level.
    Mito-level rows are pooled over all analysed cells; objects within a cell are not independent."""
    out = []
    for level, rows, xs, ys in (('cell', cell_rows, CELL_X, CELL_Y), ('mito', mito_rows, MITO_X, MITO_Y)):
        for xk, xl in available(xs, rows):
            for yk, yl in available(ys, rows):
                out.append(dict(level=level, x=xk, y=yk, **pair_stats(numeric(rows, xk), numeric(rows, yk))))
    return out
