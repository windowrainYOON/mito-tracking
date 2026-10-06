# Mito Analyzer (desktop app)

Per-cell mitochondrial network analysis for three-channel images of the same field:

| Channel | Content | Used for |
|---|---|---|
| Red | mitochondria | MiNA analysis inside each cell |
| Green | protein of interest | mean intensity per cell; its weak autofluorescence outlines cells |
| Blue | nuclei | one in-focus nucleus = one cell (watershed seeds) |

Inputs are ImageJ RGB TIFF exports with the signal in the matching channel (e.g. `…_Ch1_Red.tif`, `…_Ch2_Green.tif`, `…_Ch3_Blue.tif`); single-channel 8-bit TIFFs also work. Picking one file auto-fills the other two when the names follow that pattern.

GUI toolkit: **PySide6 (Qt)**, cross-platform (macOS, Windows, Linux).

## Pipeline
1. **Cell ROIs** (`mito_app/cell_roi.py`): nuclei by Otsu on the blue channel; watershed from each nucleus over the green-channel autofluorescence and dark cell-cell "valleys". Cells touching the image border are excluded (`edge`). Neighbours with no dark border between them are listed in `weak_border_with`; they stay split one nucleus per cell unless *Exclude cells that share a weak border* is ticked.
2. **MiNA per cell** (`mito_app/mina.py`): Python reproduction of Fiji MiNA with its defaults (no preprocessing, Otsu, skeletonize, AnalyzeSkeleton without pruning, population SD). Unlike MiNA, the threshold and skeleton use only pixels inside the cell mask, so neighbouring cells do not leak in.

## Outputs (in the output folder, prefixed with the sample name)
- `_per_cell.csv` – one row per analysed cell: MiNA metrics, green mean, centroid, border flags
- `_mina_overlay.png`, `_mina_cells.png` – full-field and per-cell overlays (magenta footprint, green skeleton, yellow ends, blue junctions)
- `_RoiSet.zip` (all cells + nuclei), `_RoiSet_filtered.zip` (QC-passed only) – open in Fiji/ImageJ ROI Manager
- `_rois.csv`, `_labels.tif`, `_summary.png`, `_overlay_*.png`, `_mito_binary.tif`, `_mito_skeleton.tif`, `_mina_per_cell.csv`

## Build the macOS app
```bash
cd application
./build_mac.sh          # creates application/Mito Analyzer.app
open "Mito Analyzer.app"
```
Requires Python 3.11+. The app is unsigned and built locally; if macOS still blocks it, right-click → Open once.

## Run from source
```bash
cd application
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_app.py                       # GUI
.venv/bin/python run_app.py --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR   # headless
```
