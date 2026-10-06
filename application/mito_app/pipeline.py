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

from . import cell_roi, mina, stats
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
    rois = cell_roi.save_outputs(a, gn, res, px, outdir, name, binuc_tau=p.binuc_tau,
                                 exclude_binuc=p.exclude_binucleate)
    n_ok = sum(r['status'] == 'ok' for r in rois)
    log(f'     {len(rois)} cell ROIs, {n_ok} pass QC (not touching the border'
        + (', not binucleate)' if p.exclude_binucleate else ')'))

    log('3/5  MiNA on the red channel inside each cell')
    lab = res['labels']
    rows, b, sk, end, junc = mina.analyze_cells(red, lab, rois, px, green, p.include_edge_cells,
                                                log=lambda s: log('     ' + s))

    log('4/5  Green on mitochondria (per cell, per mito object, green puncta)')
    mito_rows, puncta_rows, mlab, plab, thr = green_q.quantify(
        green, red, lab, b, rows, px, bg=norm['offset'], k=k, sensitivity=p.puncta_sensitivity, min_mito_area_um2=p.min_mito_area_um2,
        log=lambda s: log('     ' + s))
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
                 puncta_roiset=os.path.join(outdir, f'{name}_green_puncta_RoiSet.zip'))
    cell_rows = merge_cell_rows(rois, rows)
    settings = dict(sample=name, red=red_path, green=green_path, blue=blue_path, pixel_size_um=px,
                    green_background=norm['offset'], green_gain=norm['gain'], puncta_threshold=thr,
                    **{f'param_{k_}': v for k_, v in asdict(p).items()})
    tables = [('cells', cell_rows), ('mito', mito_rows), ('green_puncta', puncta_rows), ('correlations', corr)]
    for (_, t), f in zip(tables, (files['per_cell_csv'], files['per_mito_csv'], files['puncta_csv'], files['corr_csv'])):
        write_csv(f, t)
    write_xlsx(files['xlsx'], tables + [('settings', [settings])])
    render_green(red, green, lab, mito_rows, b, mlab, plab, rows, files['green_overlay'])
    write_label_rois(mlab, 'm', files['mito_roiset'])
    write_label_rois(plab, 'p', files['puncta_roiset'])
    log(f'     Done. Results in {outdir}')
    return dict(name=name, outdir=outdir, rows=cell_rows, mito_rows=mito_rows, puncta_rows=puncta_rows,
                correlations=corr, rois=rois, **files)


PER_CELL_ROI_COLS = ('centroid_x', 'centroid_y', 'touches_border', 'weak_border_with')


def merge_cell_rows(rois, rows):
    """One row per analysed cell: MiNA + green metrics, plus ROI geometry."""
    by = {r['roi']: r for r in rois}
    out = []
    for row in rows:
        d = dict(row)
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
