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


def why_nan(x, y, xl='X', yl='Y'):
    """Why a regression / correlation of x and y is undefined ('' if it is defined)."""
    ok = np.isfinite(x) & np.isfinite(y)
    n = int(ok.sum())
    if n < 3:
        return f'n = {n}: at least 3 values are needed'
    for v, l in ((x, xl), (y, yl)):
        if np.ptp(v[ok]) == 0:
            return f'{l} is the same ({v[ok][0]:.4g}) for all {n} values'
    return ''


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
                spearman_rho=rho, spearman_p=p_s, n=reg['n'], note=why_nan(x, y))


def correlation_table(cell_rows, mito_rows):
    """Pearson r / linear regression and Spearman rho for every mito metric x green metric pair, at cell and mito level.
    Mito-level rows are pooled over all analysed cells; objects within a cell are not independent."""
    out = []
    for level, rows, xs, ys in (('cell', cell_rows, CELL_X, CELL_Y), ('mito', mito_rows, MITO_X, MITO_Y)):
        for xk, xl in available(xs, rows):
            for yk, yl in available(ys, rows):
                out.append(dict(level=level, x=xk, y=yk, **pair_stats(numeric(rows, xk), numeric(rows, yk))))
    return out


# ---- group comparison ----
def corr_matrix(rows, xs, ys, method='r'):
    """(len(ys), len(xs)) matrix of Pearson r ('r') or Spearman rho ('rho') and the n behind each value."""
    vals = np.full((len(ys), len(xs)), np.nan); ns = np.zeros(vals.shape, int)
    Y = [numeric(rows, yk) for yk, _ in ys]   # convert each column once
    for j, (xk, _) in enumerate(xs):
        x = numeric(rows, xk)
        for i, (yk, _) in enumerate(ys):
            y = Y[i]
            if method == 'r':
                reg = regression(x, y); vals[i, j], ns[i, j] = reg['r'], reg['n']
            else:
                vals[i, j], _, ns[i, j] = spearman(x, y)
    return vals, ns


def corr_diff_p(r1, n1, r2, n2):
    """Two-sided p that two independent correlations differ (Fisher z test); nan where not computable."""
    r1, r2 = np.clip(np.asarray(r1, float), -0.999999, 0.999999), np.clip(np.asarray(r2, float), -0.999999, 0.999999)
    n1, n2 = np.asarray(n1, float), np.asarray(n2, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        se = np.sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
        z = (np.arctanh(r2) - np.arctanh(r1)) / se
    p = 2 * stats.norm.sf(np.abs(z))
    return np.where((n1 > 3) & (n2 > 3), p, np.nan)


def corr_heterogeneity_p(vals, ns):
    """Cochran's Q test that k independent correlations (one matrix per group) are equal; per cell of the
    matrices. Groups with n <= 3 or a missing value are left out of that cell. nan where < 2 groups remain."""
    V = np.clip(np.asarray(vals, float), -0.999999, 0.999999); N = np.asarray(ns, float)
    ok = np.isfinite(V) & (N > 3)
    z = np.where(ok, np.arctanh(np.where(ok, V, 0)), 0); w = np.where(ok, N - 3, 0)
    sw = w.sum(0)
    with np.errstate(divide='ignore', invalid='ignore'):
        zbar = (w * z).sum(0) / sw
        q = (w * (z - zbar) ** 2).sum(0)
    df = ok.sum(0) - 1
    p = stats.chi2.sf(q, np.maximum(df, 1))
    return np.where(df >= 1, p, np.nan)


def holm(pvals):
    """Holm-Bonferroni adjusted p values (nan entries are ignored and stay nan)."""
    p = np.asarray(pvals, float); out = np.full(p.shape, np.nan)
    idx = [i for i in np.argsort(p) if np.isfinite(p[i])]
    m = len(idx); run = 0.0
    for r, i in enumerate(idx):
        run = max(run, min(1.0, (m - r) * p[i])); out[i] = run
    return out


def pair_test(a, b):
    """Welch t and Mann-Whitney U p values for two samples (nan if not computable)."""
    a = np.asarray(a, float); b = np.asarray(b, float); a = a[np.isfinite(a)]; b = b[np.isfinite(b)]
    if len(a) < 2 or len(b) < 2 or np.ptp(np.concatenate([a, b])) == 0:
        return math.nan, math.nan
    return (float(stats.ttest_ind(a, b, equal_var=False).pvalue),
            float(stats.mannwhitneyu(a, b, alternative='two-sided').pvalue))


def group_summary_long(rows, metrics, groups, control=None, key='group'):
    """Long-format group comparison, one row per metric x group: n, mean, SD, SEM, 95 % CI, median, IQR,
    the across-group test (Welch / Mann-Whitney for 2 groups, ANOVA / Kruskal-Wallis for more) and, when a
    `control` group is given, every other group vs the control (Welch, Mann-Whitney), Holm-adjusted over the
    groups compared with it for that metric. Without a control and with 2 groups the pair is the test above."""
    out = []
    for mk, ml in metrics:
        arrays = {g: numeric([r for r in rows if str(r.get(key, '')) == g], mk) for g in groups}
        p1, p2, n1, n2 = group_test(list(arrays.values()))
        vs = {}
        if control in arrays:
            others = [g for g in groups if g != control]
            raw = [pair_test(arrays[g], arrays[control]) for g in others]
            adj_t = holm([r[0] for r in raw]); adj_u = holm([r[1] for r in raw])
            vs = {g: (raw[i][0], raw[i][1], adj_t[i], adj_u[i]) for i, g in enumerate(others)}
        for g in groups:
            d = dict(metric=ml, group=g, **describe(arrays[g]))
            d[f'p all groups ({n1})'] = p1; d[f'p all groups ({n2})'] = p2
            if control in arrays:
                t, u, ta, ua = vs.get(g, (math.nan,) * 4)
                d.update({'vs control': control if g != control else '(control)', 'p Welch vs control': t,
                          'p MWU vs control': u, 'p Welch vs control (Holm)': ta, 'p MWU vs control (Holm)': ua})
            d['column'] = mk
            out.append(d)
    return out


def describe(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    n = len(v)
    if not n:
        return dict(n=0, mean=math.nan, sd=math.nan, sem=math.nan, ci95=math.nan, median=math.nan, q1=math.nan,
                    q3=math.nan)
    sd = float(np.std(v, ddof=1)) if n > 1 else math.nan
    sem = sd / math.sqrt(n) if n > 1 else math.nan
    ci = float(stats.t.ppf(0.975, n - 1) * sem) if n > 1 else math.nan
    q1, med, q3 = np.percentile(v, [25, 50, 75])
    return dict(n=n, mean=float(v.mean()), sd=sd, sem=sem, ci95=ci, median=float(med), q1=float(q1), q3=float(q3))


def group_test(arrays):
    """Parametric and rank test across groups: Welch t + Mann-Whitney U for two groups,
    one-way ANOVA + Kruskal-Wallis for more. Returns (p_param, p_rank, name_param, name_rank)."""
    arrays = [np.asarray(a, float)[np.isfinite(a)] for a in arrays]
    arrays = [a for a in arrays if len(a)]
    two = len(arrays) == 2
    names = ('Welch t', 'Mann-Whitney U') if two else ('ANOVA', 'Kruskal-Wallis')
    if len(arrays) < 2 or any(len(a) < 2 for a in arrays):
        return math.nan, math.nan, *names
    allv = np.concatenate(arrays)
    if np.ptp(allv) == 0:
        return math.nan, math.nan, *names
    if two:
        p1 = stats.ttest_ind(*arrays, equal_var=False).pvalue
        p2 = stats.mannwhitneyu(*arrays, alternative='two-sided').pvalue
    else:
        p1 = stats.f_oneway(*arrays).pvalue
        p2 = stats.kruskal(*arrays).pvalue
    return float(p1), float(p2), *names


def sample_means(rows, keys, by=('group', 'dataset', 'preset', 'sample')):
    """One row per image (sample): the mean of each metric over its cells / mito objects."""
    out = {}
    for r in rows:
        out.setdefault(tuple(r.get(k, '') for k in by), []).append(r)
    res = []
    for kv, rs in out.items():
        d = dict(zip(by, kv), n_rows=len(rs))
        for k in keys:
            v = numeric(rs, k); v = v[np.isfinite(v)]
            d[k] = float(v.mean()) if len(v) else math.nan
        res.append(d)
    return res


def group_summary(rows, metrics, groups, key='group'):
    """Per metric: n / mean / SD / SEM / 95 % CI / median / IQR of every group and the across-group tests."""
    out = []
    for mk, ml in metrics:
        arrays = [numeric([r for r in rows if str(r.get(key, '')) == g], mk) for g in groups]
        p1, p2, n1, n2 = group_test(arrays)
        d = {'metric': ml, f'p ({n1})': p1, f'p ({n2})': p2}
        for g, a in zip(groups, arrays):
            for s, v in describe(a).items():
                d[f'{g}: {s}'] = v
        d['column'] = mk
        out.append(d)
    return out
