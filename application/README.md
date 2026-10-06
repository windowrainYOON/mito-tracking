# Mito Analyzer (desktop app)

Per-cell mitochondrial network analysis for three-channel images of the same field:

| Channel | Content | Used for |
|---|---|---|
| Red | mitochondria | MiNA analysis inside each cell |
| Green | protein of interest | mean intensity per cell; its weak autofluorescence outlines cells |
| Blue | nuclei | one in-focus nucleus = one cell (watershed seeds) |

Inputs are ImageJ RGB TIFF exports with the signal in the matching channel (e.g. `…_Ch1_Red.tif`, `…_Ch2_Green.tif`, `…_Ch3_Blue.tif`); single-channel 8-bit TIFFs also work. Files are grouped into red/green/blue sets by that naming pattern.

GUI toolkit: **PySide6 (Qt)**, cross-platform (macOS, Windows, Linux).

## Pipeline
0. **Auto-levels**: the green image is rescaled per image set (dark gap level -> 0, 75th percentile of the smoothed image -> reference level), and all size constants scale with the pixel size, so other image sets with different gain, exposure, bit depth or magnification segment the same way.
1. **Cell ROIs** (`mito_app/cell_roi.py`): nuclei by Otsu on the blue channel; watershed from each nucleus over the green-channel autofluorescence and dark cell-cell "valleys". Cells touching the image border are excluded (`edge`). Neighbours with no dark border between them are listed in `weak_border_with`; they stay split one nucleus per cell unless *Exclude cells that share a weak border* is ticked.
1b. **Green intensity split** (`mito_app/green_cells.py`, option *Split touching green+ cells by green intensity*, default on): green-positive cells and their neighbours are redrawn by a watershed whose main weight is the (log) green intensity, so borders between expressing cells follow the dark gaps of the green image. Markers are the in-focus nuclei, dim (out-of-focus) nuclei (> 0.4 × the nuclear Otsu level) and bright green territories with no nucleus (pieces within 6 µm are one cell). A new cell is kept only if it is at least the minimum cell area and green-positive; a border between two green-negative cells never moves. New cells are listed with `seed` = `dim nucleus` / `green territory`.
1c. **Green pattern fine-tuning**: the protein of interest is usually expressed by a whole cell or not at all. Bright green (two-step Otsu on the auto-levelled green, never below 3× the cytoplasm level) is grouped into patches; a cell is **green-positive** when bright green covers ≥ 2 % of its cytoplasm (option). With *Refine cell ROIs with the green pattern* (default on), a patch outside all cells, or at the edge of a cell that is green-negative without it, is moved into the green-positive cell within 2 µm (never taking a nucleus, at most +25 % of the receiving cell's area). Cells are re-scored on the refined ROIs.
2. **Mito segmentation** (option, default *Split objects*, `mito_app/mito_objects.py`): adaptive (local mean) threshold with a per-cell Otsu floor, rolling-ball background, small holes filled; intensity watershed from h-maxima, re-merged unless the contact is both dark (saddle < 0.75 × dimmer peak) and narrow (contact length < thinner width); 1-px gaps between objects. Adds fragmentation metrics per cell (objects per 100 µm² footprint, area-weighted object size, form factor, small round fraction). *MiNA classic* uses one Otsu threshold per cell instead.
3. **MiNA per cell** (`mito_app/mina.py`): Python reproduction of Fiji MiNA with its defaults (no preprocessing, Otsu, skeletonize, AnalyzeSkeleton without pruning, population SD). Unlike MiNA, the threshold and skeleton use only pixels inside the cell mask, so neighbouring cells do not leak in.

4. **Green on mito** (`mito_app/green.py`):
   - *mito objects* = connected pieces of the MiNA mito mask inside a cell (pieces under 0.05 µm² are left out); each gets skeleton length, branches, junctions, aspect ratio, solidity.
   - *green puncta* = white top-hat of green, one threshold per image (max of Otsu and median + 6 MAD inside analysed cells, × the sensitivity option).
   - per cell: green mean on / off mito, enrichment, fraction of green on mito, Pearson, Manders, puncta counts and density.
   - per mito object: green mean / integrated / max on it, number of overlapping puncta, their area, coverage and intensity.
   - Pearson r with the regression line (slope, intercept, R², p) and Spearman ρ of every mito metric × green metric, per cell and per mito object.
5. **Green-positive and green-negative cells are reported separately**: every table gets a `green_status` column and is also written as a `_green_pos` and a `_green_neg` version, with the correlations computed within each group only. The puncta threshold stays one per image.

## In the app
**Image sets table** (batch input): add files or whole folders (*Add files… / Add folder…*, or drop them on the table; folders are searched recursively, so `PRD 24h/Processed/New-01_results/…` is found from `PRD 24h`); every complete red/green/blue set becomes one row with editable *Dataset* (default: the added folder's name, or the file's folder), *Preset* (default: the mito method, `split` / `otsu`; type e.g. a condition name instead) and *Sample* (from the file name). Double-click a channel cell to swap its file. *Run all* analyses the rows one after another (*Stop* finishes the current set and stops); a failed set is marked and the batch goes on. Click a finished row to show its results. Regression / correlation needs at least 3 cells (or objects) with varying values in the current selection; otherwise the app says why instead of showing a number (heatmap squares show –, tables have a `note` column).

The **Analysis (all samples)** tab loads every `<sample>_per_cell.csv` / `_per_mito.csv` under a results folder (filled in automatically after a batch) and pools them with `dataset`, `preset` and `sample` columns: filter by green+ / green−, dataset, preset or sample, colour by any of them, and export the pooled tables plus a regression table (`pooled_*`).

Correlation views (per sample and pooled): heatmap of Pearson r (or Spearman ρ) for every mito × green pair; click a square for its scatter with the least-squares line, 95 % band, r, R², p and n, split into quadrants at the mean (or median) of X and Y with the share of points in each, so negative relations (quadrants II / IV) stand out. With log X / log Y the fit uses log10 values.

Tabs: **Cells** (all per-cell metrics; column headers show `(+)` / `(−)`), **Mito objects** and **Green puncta** (sortable, filter by green+ / green− and by cell), **Correlation** (Spearman heatmap of all pairs within all / green+ / green− cells; click a square for the scatter plot, coloured by cell), **Green+ / − cells** (ROI refinement and the green calls), and overlay images.

## Outputs
Folder layout: `<output>/<dataset>/<preset>-<dataset>/<sample>/`. Each `<preset>-<dataset>` folder also gets the pooled tables of its samples: `<preset>-<dataset>_all_cells.csv`, `…_all_mito.csv`, their `_green_pos` / `_green_neg` versions, and `…_all_results.xlsx` (with a `sample` column).

Per sample (in its own folder, prefixed with the sample name):
- `_results.xlsx` – sheets `cells_green_pos`, `mito_green_pos`, `puncta_green_pos`, `corr_green_pos`, the same four for `green_neg`, then `all_*` (every analysed cell), `roi_refinement` (patches moved between ROIs) and `settings`
- `_per_cell_green_pos.csv`, `_per_mito_green_pos.csv`, `_green_puncta_green_pos.csv`, `_correlations_green_pos.csv`, the same for `_green_neg`, and the all-cell `_per_cell.csv`, `_per_mito.csv`, `_green_puncta.csv`, `_correlations.csv`
- `_green_cells.png` – green+ (red) / green− (grey) cells, bright green (yellow), ROIs before refinement (blue dashed)
- `_green_on_mito.png` – mito objects (magenta), puncta on mito (yellow) / off mito (cyan)
- `_mito_RoiSet.zip`, `_green_puncta_RoiSet.zip` – mito objects and puncta as ImageJ ROIs
- `_mina_overlay.png`, `_mina_cells.png` – full-field and per-cell overlays (magenta footprint, green skeleton, yellow ends, blue junctions)
- `_RoiSet.zip` (all cells + nuclei), `_RoiSet_filtered.zip` (QC-passed only) – open in Fiji/ImageJ ROI Manager
- `_rois.csv`, `_labels.tif`, `_summary.png`, `_overlay_*.png`, `_mito_binary.tif`, `_mito_skeleton.tif`, `_mina_per_cell.csv`

## Build the macOS app
```bash
cd application
./build_mac.sh                 # creates application/Mito Analyzer.app
./build_mac.sh --add-to-dock   # same, and pins it to the Dock once (rebuilds keep the same path)
open "Mito Analyzer.app"
```
Requires Python 3.11+. The app is unsigned and built locally; if macOS still blocks it, right-click → Open once.

## Run from source
```bash
cd application
python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
.venv/bin/python run_app.py                       # GUI
.venv/bin/python run_app.py --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR   # headless, one set
.venv/bin/python run_app.py --cli --batch FOLDER [FOLDER…] -o OUTDIR [--dataset NAME] [--preset NAME]  # headless batch
```

## Icon
`resources/icon.{png,icns,ico}` are generated by `python tools/make_icon.py`.

## Branches
Development happens on `develop`; the Mac build is rebuilt from it after each change.
