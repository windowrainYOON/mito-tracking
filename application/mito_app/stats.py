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
          ('mito_aspect_ratio_mean', 'Mito aspect ratio mean')]
CELL_Y = [('green_poi_mean_cell', 'Green mean, whole cell'), ('green_mean_on_mito', 'Green mean on mito'),
          ('green_mean_off_mito', 'Green mean off mito'), ('green_enrichment_on_mito', 'Green enrichment on/off mito'),
          ('green_fraction_on_mito', 'Fraction of green on mito'), ('green_integrated_on_mito', 'Green integrated on mito'),
          ('pearson_red_green', 'Pearson red–green'), ('manders_m_green', 'Manders: puncta green on mito'),
          ('manders_m_mito', 'Manders: mito red under puncta'), ('n_green_puncta', 'Green puncta'),
          ('n_puncta_on_mito', 'Green puncta on mito'), ('fraction_puncta_on_mito', 'Fraction of puncta on mito'),
          ('green_puncta_density_per_100um2', 'Puncta per 100 µm²'),
          ('fraction_mito_with_puncta', 'Fraction of mito objects with puncta')]
MITO_X = [('length_um', 'Length (µm)'), ('aspect_ratio', 'Aspect ratio'), ('area_um2', 'Area (µm²)'),
          ('n_branches', 'Branches'), ('n_junctions', 'Junctions'), ('n_endpoints', 'End points'),
          ('solidity', 'Solidity'), ('major_axis_um', 'Major axis (µm)'), ('red_mean', 'Red mean')]
MITO_Y = [('green_mean', 'Green mean'), ('green_integrated', 'Green integrated'), ('green_max', 'Green max'),
          ('n_green_puncta', 'Green puncta'), ('puncta_coverage', 'Puncta coverage'),
          ('puncta_overlap_area_um2', 'Puncta overlap area (µm²)'),
          ('puncta_green_integrated', 'Puncta green integrated'), ('puncta_mean_area_um2', 'Puncta mean size (µm²)')]


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


def correlation_table(cell_rows, mito_rows):
    """Spearman rho for every mito metric x green metric pair, at cell and mito level.
    Mito-level rows are pooled over all analysed cells; objects within a cell are not independent."""
    out = []
    for level, rows, xs, ys in (('cell', cell_rows, CELL_X, CELL_Y), ('mito', mito_rows, MITO_X, MITO_Y)):
        for xk, xl in xs:
            for yk, yl in ys:
                rho, p, n = spearman(numeric(rows, xk), numeric(rows, yk))
                out.append(dict(level=level, x=xk, y=yk, spearman_rho=rho, p_value=p, n=n))
    return out
