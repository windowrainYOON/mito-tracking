"""End-to-end pipeline: nucleus-seeded cell ROIs -> per-cell MiNA on the red channel.

Inputs are three images of the same field (ImageJ RGB TIFF exports, signal in the
matching channel, or single-channel images):
  red   = mitochondria (MiNA is run on it)
  green = protein of interest (mean intensity per cell; its weak autofluorescence
          also outlines the cells for the ROI step)
  blue  = nuclei (one in-focus nucleus = one cell)
"""
import csv, os, re
from dataclasses import dataclass

import numpy as np

from . import cell_roi, mina


@dataclass
class Params:
    min_area_um2: float = 250.0       # smallest cell ROI kept (cell_roi --min-area-um2)
    binuc_tau: float = 0.0            # border contrast below which neighbours are "weak border" pairs
    exclude_binucleate: bool = False  # drop weak-border pairs instead of keeping them split per nucleus
    include_edge_cells: bool = False  # also run MiNA on cells touching the frame


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

    log('1/4  Loading images')
    a, g, px = cell_roi.load_green(green_path)
    a = _as_rgb(a, g)
    nuc_img, nuc = cell_roi.load_nuclei(blue_path)
    red, px_red = mina.load_channel(red_path, 0)
    green, _ = mina.load_channel(green_path, 1)
    if not (red.shape == g.shape == nuc_img.shape):
        raise ValueError(f'Image sizes differ: red {red.shape}, green {g.shape}, blue {nuc_img.shape}')
    log(f'     {red.shape[1]}x{red.shape[0]} px, pixel size {px:.4f} um, {nuc.max()} in-focus nuclei')

    log('2/4  Segmenting cells (one nucleus = one cell)')
    res = cell_roi.segment(g, min_area_px=int(p.min_area_um2 / px ** 2), nuclei=nuc)
    res['nuc_img'] = nuc_img
    rois = cell_roi.save_outputs(a, g, res, px, outdir, name, binuc_tau=p.binuc_tau,
                                 exclude_binuc=p.exclude_binucleate)
    n_ok = sum(r['status'] == 'ok' for r in rois)
    log(f'     {len(rois)} cell ROIs, {n_ok} pass QC (not touching the border'
        + (', not binucleate)' if p.exclude_binucleate else ')'))

    log('3/4  MiNA on the red channel inside each cell')
    lab = res['labels']
    rows, b, sk, end, junc = mina.analyze_cells(red, lab, rois, px, green, p.include_edge_cells,
                                                log=lambda s: log('     ' + s))

    log('4/4  Writing results')
    mina.write_outputs(red, lab, rois, rows, b, sk, end, junc, outdir, name)
    per_cell = os.path.join(outdir, f'{name}_per_cell.csv')
    write_per_cell(per_cell, rois, rows)
    log(f'     Done. Results in {outdir}')
    return dict(name=name, outdir=outdir, rows=rows, rois=rois, per_cell_csv=per_cell,
                overlay=os.path.join(outdir, f'{name}_mina_overlay.png'),
                cells=os.path.join(outdir, f'{name}_mina_cells.png'),
                summary=os.path.join(outdir, f'{name}_summary.png'),
                roiset=os.path.join(outdir, f'{name}_RoiSet.zip'),
                roiset_filtered=os.path.join(outdir, f'{name}_RoiSet_filtered.zip'))


PER_CELL_ROI_COLS = ('centroid_x', 'centroid_y', 'touches_border', 'weak_border_with')


def write_per_cell(path, rois, rows):
    """One row per analysed cell: ROI geometry + MiNA metrics + green mean."""
    by = {r['roi']: r for r in rois}
    out = []
    for row in rows:
        r = by[row['cell']]
        d = dict(row)
        for k in PER_CELL_ROI_COLS:
            d[k] = r.get(k, '')
        out.append(d)
    with open(path, 'w', newline='') as f:
        if not out:
            f.write('cell\n'); return
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        for d in out:
            w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in d.items()})
