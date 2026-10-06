"""End-to-end pipeline: nucleus-seeded cell ROIs -> per-cell MiNA on the red channel.

Inputs are three images of the same field (ImageJ RGB TIFF exports, signal in the
matching channel, or single-channel images):
  red   = mitochondria (MiNA is run on it)
  green = protein of interest (mean intensity per cell; its weak autofluorescence
          also outlines the cells for the ROI step)
  blue  = nuclei (one in-focus nucleus = one cell)
"""
import csv, os, re
from dataclasses import asdict, dataclass

import numpy as np

from . import cell_roi, green_cells, mina, mito_objects, stats
from . import green as green_q


@dataclass
class Params:
    min_area_um2: float = 250.0       # smallest cell ROI kept (cell_roi --min-area-um2)
    binuc_tau: float = 0.0            # border contrast below which neighbours are "weak border" pairs
    exclude_binucleate: bool = False  # drop weak-border pairs instead of keeping them split per nucleus
    include_edge_cells: bool = False  # also run MiNA on cells touching the frame
    pixel_size_um: float = 0.0        # 0 = read from the TIFF (falls back to the reference size)
    puncta_sensitivity: float = 1.0   # multiplies the green puncta threshold (>1 = fewer, brighter puncta)
    min_mito_area_um2: float = 0.05   # smallest mito object in the per-mito table
    mito_method: str = 'split'        # 'split' = adaptive threshold + watershed objects, 'otsu' = MiNA classic
    green_refine: bool = True         # extend cell ROIs over nearby green patches of expressing cells
    green_pos_percent: float = 2.0    # cell is green-positive if bright green covers >= this % of its cytoplasm


def sample_name(red_path):
    """'Processed_New01_Ch1_Red.tif' -> 'Processed_New01'."""
    stem = os.path.splitext(os.path.basename(red_path))[0]
    return re.sub(r'[_\- ]*(ch\d+)?[_\- ]*(red|mito)$', '', stem, flags=re.I) or stem


def find_siblings(path):
    """Guess the other two channel files from one file name (…Ch1_Red / …Ch2_Green / …Ch3_Blue)."""
    d, name = os.path.split(path)
    m = re.search(r'(ch\d+[_\- ]*)?(red|green|blue)', name, re.I)
    if not m:
        return {}
    files = os.listdir(d or '.')
    out = {}
    for color in ('red', 'green', 'blue'):
        rx = re.compile(re.escape(name[:m.start()]) + r'(ch\d+[_\- ]*)?' + color + re.escape(name[m.end():]), re.I)
        hits = [f for f in files if rx.fullmatch(f)]
        if len(hits) == 1:
            out[color] = os.path.join(d, hits[0])
    return out


def find_sets(paths):
    """Group TIFF files (or the TIFFs inside folders) into red/green/blue image sets by file name.

    Returns (sets, leftover): sets is a list of {'red', 'green', 'blue'} dicts in name order,
    leftover the files that could not be placed in a complete set."""
    files = []
    for p in paths:
        if os.path.isdir(p):
            files += [os.path.join(p, f) for f in sorted(os.listdir(p))
                      if f.lower().endswith(('.tif', '.tiff')) and not f.startswith('.')]
        elif os.path.isfile(p):
            files.append(p)
    sets, used = {}, set()
    for f in files:
        s = find_siblings(f)
        if len(s) == 3:
            sets[s['red']] = s
            used.update(s.values())
    leftover = [f for f in files if f not in used]
    return [sets[k] for k in sorted(sets)], leftover


def safe_name(s):
    """Folder-safe version of a dataset / preset / sample name."""
    return re.sub(r'[\\/:*?"<>|]+', '_', str(s)).strip(' .') or 'unnamed'


def job_outdir(base, dataset, preset, name):
    """<output>/<dataset>/<preset>-<dataset>/<sample>/"""
    d, p = safe_name(dataset), safe_name(preset)
    return os.path.join(base, d, f'{p}-{d}', safe_name(name))


def write_group_tables(results):
    """Pool the per-cell / per-mito tables of all samples in each <preset>-<dataset> folder
    (`<preset>-<dataset>_all_cells.csv`, `…_all_mito.csv`, `…_all_results.xlsx`)."""
    groups = {}
    for r in results:
        groups.setdefault(os.path.dirname(r['outdir']), []).append(r)
    for d, rs in groups.items():
        tag = os.path.basename(d)
        cells = [dict(sample=r['name'], **row) for r in rs for row in r['rows']]
        mito = [dict(sample=r['name'], **row) for r in rs for row in r['mito_rows']]
        write_csv(os.path.join(d, f'{tag}_all_cells.csv'), cells)
        write_csv(os.path.join(d, f'{tag}_all_mito.csv'), mito)
        sheets = []
        for g_key, g_tag in GREEN_GROUPS:
            gc = [r for r in cells if r.get('green_status') == g_key]
            gm = [r for r in mito if r.get('green_status') == g_key]
            write_csv(os.path.join(d, f'{tag}_all_cells_{g_tag}.csv'), gc)
            write_csv(os.path.join(d, f'{tag}_all_mito_{g_tag}.csv'), gm)
            sheets += [(f'cells_{g_tag}', gc), (f'mito_{g_tag}', gm)]
        write_xlsx(os.path.join(d, f'{tag}_all_results.xlsx'), sheets + [('all_cells', cells), ('all_mito', mito)])


def _as_rgb(a, g):
    """cell_roi overlays expect an RGB uint8 array; wrap single-channel input as green."""
    if a.ndim == 3 and a.shape[-1] == 3:
        return a
    g8 = np.clip(g / max(np.percentile(g, 99.5), 1) * 255, 0, 255).astype(np.uint8)
    return np.dstack([np.zeros_like(g8), g8, np.zeros_like(g8)])


def run(red_path, green_path, blue_path, outdir, params=None, name=None, log=print):
    """Run the full pipeline and write all outputs to `outdir`. Returns a dict of results and file paths."""
    p = params or Params()
    name = name or sample_name(red_path)
    os.makedirs(outdir, exist_ok=True)

    log('1/5  Loading images')
    a, g, px = cell_roi.load_green(green_path)
    if p.pixel_size_um > 0:
        px = p.pixel_size_um
    elif px == 1.0:  # no calibration in the file
        px = cell_roi.REF_PX_UM
        log(f'     No pixel size in the TIFF; assuming {px} um/px (set it in Options if different)')
    k = float(np.clip(cell_roi.REF_PX_UM / px, 0.2, 5))
    a = _as_rgb(a, g)
    nuc_img, nuc = cell_roi.load_nuclei(blue_path, k)
    red, _ = mina.load_channel(red_path, 0)
    green, _ = mina.load_channel(green_path, 1)
    if not (red.shape == g.shape == nuc_img.shape):
        raise ValueError(f'Image sizes differ: red {red.shape}, green {g.shape}, blue {nuc_img.shape}')
    gn, norm = cell_roi.normalize_green(g, k)
    log(f'     {red.shape[1]}x{red.shape[0]} px, pixel size {px:.4f} um, {nuc.max()} in-focus nuclei; '
        f'green auto-levels: background {norm["offset"]:.2f}, gain x{norm["gain"]:.2f}')

    log('2/5  Segmenting cells (one nucleus = one cell)')
    res = cell_roi.segment(gn, min_area_px=int(p.min_area_um2 / px ** 2), nuclei=nuc, k=k)
    res['nuc_img'] = nuc_img
    lab0 = res['labels']
    log('     Green pattern: ' + ('refining cell ROIs and ' if p.green_refine else '') + 'calling green+/- cells')
    res['labels'], gstatus, ginfo = green_cells.refine(
        gn, lab0, nuc, px, k, min_fraction=p.green_pos_percent / 100, do_refine=p.green_refine,
        log=lambda s: log('     ' + s))
    rois = cell_roi.save_outputs(a, gn, res, px, outdir, name, binuc_tau=p.binuc_tau,
                                 exclude_binuc=p.exclude_binucleate)
    n_ok = sum(r['status'] == 'ok' for r in rois)
    for r in rois:
        r['green_status'] = gstatus[int(r['roi'][4:])]
    log(f'     {len(rois)} cell ROIs, {n_ok} pass QC (not touching the border'
        + (', not binucleate)' if p.exclude_binucleate else ')'))

    lab = res['labels']
    mask_fn = None
    if p.mito_method == 'split':
        log('3/5  Mito objects (adaptive threshold + watershed split), then MiNA inside each cell')
        p0 = dict(mito_objects.P, sigma=mito_objects.P['sigma'] * k)
        sm = mito_objects.preprocess(red, px, p0)
        used = [int(r['roi'][4:]) for r in rois if r['status'] == 'ok' or p.include_edge_cells]
        q = mito_objects.scaled_params(sm, np.isin(lab, used), k, p0)

        def mask_fn(cell):
            fg, _, sep = mito_objects.segment_cell(sm, cell, px, q)
            return sep > 0, mito_objects.fragmentation_metrics(fg, sep, cell, px)
    else:
        log('3/5  MiNA (Otsu per cell) on the red channel inside each cell')
    rows, b, sk, end, junc = mina.analyze_cells(red, lab, rois, px, green, p.include_edge_cells,
                                                log=lambda s: log('     ' + s), mask_fn=mask_fn)

    log('4/5  Green on mitochondria (per cell, per mito object, green puncta)')
    mito_rows, puncta_rows, mlab, plab, thr = green_q.quantify(
        green, red, lab, b, rows, px, bg=norm['offset'], k=k, sensitivity=p.puncta_sensitivity, min_mito_area_um2=p.min_mito_area_um2,
        log=lambda s: log('     ' + s))
    for row in rows:
        c = int(row['cell'][4:])
        row['green_status'] = gstatus[c]
        row['green_bright_area_percent'] = 100 * ginfo['fraction'][c]
    cstat = {row['cell']: row['green_status'] for row in rows}
    for t in (mito_rows, puncta_rows):
        for row in t:
            row['green_status'] = cstat.get(row['cell'], '')
    n_pos = sum(s == 'positive' for s in cstat.values())
    log(f'     {n_pos} green-positive, {len(cstat) - n_pos} green-negative analysed cells '
        f'(bright green >= {ginfo["threshold"]:.1f} auto-levelled units)')
    corr = stats.correlation_table(rows, mito_rows)

    log('5/5  Writing results')
    mina.write_outputs(red, lab, rois, rows, b, sk, end, junc, outdir, name)
    files = dict(per_cell_csv=os.path.join(outdir, f'{name}_per_cell.csv'),
                 per_mito_csv=os.path.join(outdir, f'{name}_per_mito.csv'),
                 puncta_csv=os.path.join(outdir, f'{name}_green_puncta.csv'),
                 corr_csv=os.path.join(outdir, f'{name}_correlations.csv'),
                 xlsx=os.path.join(outdir, f'{name}_results.xlsx'),
                 overlay=os.path.join(outdir, f'{name}_mina_overlay.png'),
                 green_overlay=os.path.join(outdir, f'{name}_green_on_mito.png'),
                 cells=os.path.join(outdir, f'{name}_mina_cells.png'),
                 summary=os.path.join(outdir, f'{name}_summary.png'),
                 roiset=os.path.join(outdir, f'{name}_RoiSet.zip'),
                 roiset_filtered=os.path.join(outdir, f'{name}_RoiSet_filtered.zip'),
                 mito_roiset=os.path.join(outdir, f'{name}_mito_RoiSet.zip'),
                 puncta_roiset=os.path.join(outdir, f'{name}_green_puncta_RoiSet.zip'),
                 green_cells=os.path.join(outdir, f'{name}_green_cells.png'))
    cell_rows = merge_cell_rows(rois, rows)
    settings = dict(sample=name, red=red_path, green=green_path, blue=blue_path, pixel_size_um=px,
                    green_background=norm['offset'], green_gain=norm['gain'], puncta_threshold=thr,
                    green_bright_threshold=ginfo['threshold'],
                    **{f'param_{k_}': v for k_, v in asdict(p).items()})
    tables = [('cells', cell_rows), ('mito', mito_rows), ('green_puncta', puncta_rows), ('correlations', corr)]
    for (_, t), f in zip(tables, (files['per_cell_csv'], files['per_mito_csv'], files['puncta_csv'], files['corr_csv'])):
        write_csv(f, t)
    # green-positive and green-negative cells, measured and reported separately
    groups = {}
    for g_key, g_tag in GREEN_GROUPS:
        gc = [r for r in cell_rows if r['green_status'] == g_key]
        gm = [r for r in mito_rows if r['green_status'] == g_key]
        gp = [r for r in puncta_rows if r['green_status'] == g_key]
        gr = stats.correlation_table(gc, gm)
        groups[g_key] = dict(rows=gc, mito_rows=gm, puncta_rows=gp, correlations=gr)
        for label, t in (('per_cell', gc), ('per_mito', gm), ('green_puncta', gp), ('correlations', gr)):
            write_csv(os.path.join(outdir, f'{name}_{label}_{g_tag}.csv'), t)
    group_sheets = [(f'{s}_{g_tag}', groups[g_key][k]) for g_key, g_tag in GREEN_GROUPS
                    for s, k in (('cells', 'rows'), ('mito', 'mito_rows'), ('puncta', 'puncta_rows'), ('corr', 'correlations'))]
    write_xlsx(files['xlsx'], group_sheets + [(f'all_{n_}', t) for n_, t in tables]
               + [('roi_refinement', ginfo['moved']), ('settings', [settings])])
    render_green_cells(gn, lab0, res['labels'], ginfo, gstatus, {r['roi'] for r in rois if r['status'] == 'ok'
                       or p.include_edge_cells}, files['green_cells'])
    render_green(red, green, lab, mito_rows, b, mlab, plab, rows, files['green_overlay'])
    write_label_rois(mlab, 'm', files['mito_roiset'])
    write_label_rois(plab, 'p', files['puncta_roiset'])
    log(f'     Done. Results in {outdir}')
    return dict(name=name, outdir=outdir, rows=cell_rows, mito_rows=mito_rows, puncta_rows=puncta_rows,
                correlations=corr, rois=rois, groups=groups, roi_moves=ginfo['moved'], **files)


GREEN_GROUPS = (('positive', 'green_pos'), ('negative', 'green_neg'))


PER_CELL_ROI_COLS = ('centroid_x', 'centroid_y', 'touches_border', 'weak_border_with')


def merge_cell_rows(rois, rows):
    """One row per analysed cell: MiNA + green metrics, plus ROI geometry."""
    by = {r['roi']: r for r in rois}
    out = []
    for row in rows:
        d = {'cell': row['cell'], 'green_status': row['green_status'], **row}
        for k in PER_CELL_ROI_COLS:
            d[k] = by[row['cell']].get(k, '')
        out.append(d)
    return out


def _fmt(v):
    return f'{v:.4f}' if isinstance(v, float) else v


def write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        if not rows:
            return
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        for d in rows:
            w.writerow({k: _fmt(v) for k, v in d.items()})


def write_xlsx(path, sheets):
    from openpyxl import Workbook
    wb = Workbook(); wb.remove(wb.active)
    for title, rows in sheets:
        ws = wb.create_sheet(title)
        if not rows:
            continue
        cols = list(rows[0])
        ws.append(cols)
        for r in rows:
            ws.append([(None if isinstance(v, float) and v != v else (bool(v) if isinstance(v, np.bool_) else v))
                       for v in (r.get(c) for c in cols)])
        ws.freeze_panes = 'B2'
        for i, c in enumerate(cols, 1):
            ws.column_dimensions[ws.cell(1, i).column_letter].width = max(9, min(32, len(c) + 2))
    wb.save(path)


def write_label_rois(lab, prefix, path):
    """ImageJ RoiSet with one polygon per labelled object (mito objects / green puncta)."""
    import roifile
    from skimage import measure
    rois = []
    for r in measure.regionprops(lab):
        c = max(measure.find_contours(np.pad(r.image, 1).astype(float), 0.5), key=len) - 1
        c = c + np.array(r.bbox[:2])
        roi = roifile.ImagejRoi.frompoints(np.c_[c[:, 1], c[:, 0]])
        roi.roitype = roifile.ROI_TYPE.POLYGON; roi.name = f'{prefix}{r.label:04d}'
        rois.append(roi)
    if os.path.exists(path):
        os.remove(path)
    if rois:
        roifile.roiwrite(path, rois)


def render_green(red, green, lab, mito_rows, mito, mlab, plab, rows, path):
    """Red mito (grey) + green POI, mito outlines magenta, puncta on mito yellow, off mito cyan."""
    from matplotlib.figure import Figure
    from scipy import ndimage as ndi
    from skimage import segmentation
    rh, gh = np.percentile(red, 99.5) or 1, np.percentile(green, 99.8) or 1
    rgb = np.dstack([np.clip(red / rh, 0, 1) * 0.55, np.clip(green / gh, 0, 1), np.clip(red / rh, 0, 1) * 0.55])
    rgb[segmentation.find_boundaries(mlab, mode='outer')] = (1.0, 0.2, 0.9)
    on = np.isin(plab, np.unique(plab[(plab > 0) & mito]))
    pb = segmentation.find_boundaries(plab, mode='outer')
    rgb[pb & ndi.binary_dilation(on, iterations=1)] = (1.0, 1.0, 0.0)
    rgb[pb & ~ndi.binary_dilation(on, iterations=1)] = (0.0, 0.9, 1.0)
    analysed = np.isin(lab, [int(r['cell'][4:]) for r in rows])
    rgb[segmentation.find_boundaries(np.where(analysed, lab, 0), mode='inner')] = (1, 1, 1)
    fig = Figure(figsize=(12, 12), dpi=150); ax = fig.subplots()
    ax.imshow(rgb); ax.set_axis_off()
    for r in rows:
        cy, cx = ndi.center_of_mass(lab == int(r['cell'][4:]))
        ax.text(cx, cy, r['cell'][4:], color='white', fontsize=13, weight='bold', ha='center', va='center',
                bbox=dict(facecolor='black', alpha=0.5, pad=1.5, edgecolor='none'))
    ax.set_title('Green POI on mito. grey/purple: red mito, green: POI; magenta: mito objects; '
                 'yellow: green puncta on mito, cyan: puncta off mito; white: analysed cells', fontsize=9)
    fig.tight_layout(); fig.savefig(path)


def render_green_cells(gn, lab0, lab, info, status, analysed, path):
    """Bright green mask, cell ROIs before (dashed) and after green refinement, green+ / green- labels."""
    from matplotlib.figure import Figure
    from scipy import ndimage as ndi
    v = np.clip(gn / max(np.percentile(gn, 99.5), 1e-9), 0, 1)
    rgb = np.dstack([v * 0.15, v, v * 0.15])
    rgb[info['mask']] = (1.0, 0.95, 0.2)
    fig = Figure(figsize=(12, 12), dpi=120); ax = fig.subplots()
    ax.imshow(rgb); ax.set_axis_off()
    if (lab0 != lab).any():
        ax.contour(lab0 > 0, [0.5], colors='#7fa7ff', linewidths=0.8, linestyles='--')
        for l in np.unique(lab0[lab0 != lab]):
            ax.contour(lab0 == l, [0.5], colors='#7fa7ff', linewidths=1.0, linestyles='--')
    for l in range(1, lab.max() + 1):
        m = lab == l
        if not m.any():
            continue
        pos = status[l] == 'positive'
        ax.contour(m, [0.5], colors='#ff5050' if pos else '#c8c8c8', linewidths=2.2 if pos else 1.0)
        cy, cx = ndi.center_of_mass(m)
        name = f'cell{l:02d}'
        ax.text(cx, cy, f"{l}{'+' if pos else '−'}" + ('' if name in analysed else '\n(not analysed)'),
                color='white', fontsize=(13 if pos else 10) if name in analysed else 8, weight='bold', ha='center', va='center',
                bbox=dict(facecolor='#b00000' if pos else 'black', alpha=0.6, pad=1.5, edgecolor='none'))
    n_mv = len(info['moved'])
    ax.set_title(f'Green+ cells (red outline) / green− (grey). Yellow: bright green >= {info["threshold"]:.1f}. '
                 + (f'{n_mv} green patch(es) moved into the expressing cell; blue dashed = ROI before refinement.'
                    if n_mv else 'No ROI changed by the green pattern.'), fontsize=9)
    fig.tight_layout(); fig.savefig(path)
