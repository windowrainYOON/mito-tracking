# Mito Analyzer (desktop app)

Per-cell mitochondrial network analysis for three-channel images of the same field:

| Channel | Content | Used for |
|---|---|---|
| Red | mitochondria | MiNA analysis inside each cell |
| Green | protein of interest | green on mito, green puncta, green+ / green− call; its weak autofluorescence helps place cell borders |
| Blue | nuclei | one nucleus (in focus or dim) = one cell (watershed seeds) |

Inputs are ImageJ RGB TIFF exports with the signal in the matching channel (e.g. `…_Ch1_Red.tif`, `…_Ch2_Green.tif`, `…_Ch3_Blue.tif`); single-channel 8-bit TIFFs also work. Files are grouped into red/green/blue sets by that naming pattern.

GUI toolkit: **PySide6 (Qt)**, cross-platform (macOS, Windows, Linux).

## Pipeline
Every metric in the tables and heatmaps is explained (in Korean) in [`docs/metrics.md`](docs/metrics.md); the full algorithm, with the cell-ROI step in detail, in [`docs/algorithm.md`](docs/algorithm.md).

A step-by-step explanation with figures from one example image is in [`docs/Mito_Analyzer_algorithm.pdf`](docs/Mito_Analyzer_algorithm.pdf) (regenerate with `tools/make_algorithm_pdf.py RED GREEN BLUE -o docs/Mito_Analyzer_algorithm.pdf`).

0. **Auto-levels**: the green image is rescaled per image set (dark gap level -> 0, 75th percentile of the smoothed image -> reference level), and all size constants scale with the pixel size, so other image sets with different gain, exposure, bit depth or magnification segment the same way.
1. **Cell ROIs** (`mito_app/cell_roi.py`, shape + intensity): each cell has one nucleus at its centre with the cytoplasm spread around it, and neighbouring cells are separated by thin lines where the mitochondria stop.
   - seeds (`mito_app/nuclei.py`, option *Nucleus detection*, default *Texture + merge*): real nuclei are told from diffuse cytoplasmic blue by their chromatin texture. Blue is smoothed (σ 2 px), a 25 µm background is subtracted, candidates = 0.5 × Otsu, split by a distance watershed; each piece is kept if its noise-corrected texture (SD / median) is ≥ 0.22, or ≥ 0.17 and it is as bright as the image's textured nuclei (≥ 0.7 × their median). Kept neighbouring pieces are merged unless their contact is both dark (< 0.8 × the smaller piece) and narrow (< 0.6 × its minor axis), so curved nuclei stay whole and touching nuclei stay apart. Each nucleus then gets its own threshold (halfway between core and surrounding ring), rim bays are closed inside the hull, non-compact regions (solidity < 0.8) are dropped. Nuclei dimmer than 0.5 × the reference (or small) are dim (out-of-focus) nuclei, kept if ≥ max(40 µm², half the median in-focus area); option *Dim nuclei also get their own cell*, default on, so their mitochondria are not given to a neighbour. ROIs grown from them have `seed` = `dim nucleus`. *Texture, no merge* skips the merge step (curved nuclei may come out in pieces); *Global Otsu (old)* is the earlier one-threshold method.
   - foreground: mitochondria density (red smoothed at 3 µm) above 12 % of its 99th percentile, plus the nuclei.
   - borders: compact watershed from the nuclei over a cost that is low where mitochondria are dense, high on the dark mito-free lines (black-ridge filter at 1.5–2.5 µm), plus the green-autofluorescence valleys. The compactness keeps each cell growing evenly around its nucleus.
   - Edge tips: a cell whose nucleus is well inside the frame (≥ 3 µm) but whose ROI reaches the border is split again by a two-seed watershed on the same cost (nucleus + 2 µm ring vs. the border pixels). The border-side part is cut off when the cut follows a mito-free line (median mito density on it ≤ 0.8 × the kept cytoplasm, or the part is a sliver < 10 % of the ROI), the kept cell holds ≥ 55 % of the ROI, is compact (solidity ≥ 0.70), its nucleus is ≥ 3 µm from the cut, and it no longer runs along the border. Such cells are analysed; `edge_trimmed_um2` gives the area cut off, hatched in the summary (option *Keep cells whose tip touches the border*, default on; CLI `--no-edge-trim`).
   - Cells still touching the image border are excluded (`edge`). Neighbours with no dark border between them are listed in `weak_border_with`; they stay split one nucleus per cell unless *Exclude cells that share a weak border* is ticked.
   - **Green+ / green− call** (`mito_app/green_cells.py`, does not change the ROIs): bright green (two-step Otsu on the auto-levelled green, never below 3× the cytoplasm level); a cell is green-positive when bright green covers ≥ 2 % of its cytoplasm (option).
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

**Groups.** Each image set has a *Group* (default: its dataset); select rows and press *Set group…* to put several sets in one group. Groups are saved in `<output>/groups.csv` and can be changed later with *Groups…* in the Analysis tab (CLI: `--group NAME`). The Analysis tab compares them:
- *Groups: correlation*: the mito × green correlation heatmaps of group A and group B side by side, and a third heatmap of the difference (B − A) with `*` p < 0.05 / `**` p < 0.01 (Fisher z test for two independent correlations); the largest differences are listed below.
- *Groups: mean / median*: per metric, the individual values with mean ± SD / SEM / 95 % CI or median ± IQR per group, and a table of n, mean, SD, SEM, 95 % CI, median, Q1, Q3 per group with Welch t + Mann-Whitney U (two groups) or one-way ANOVA + Kruskal-Wallis (more). The unit can be cells / mito objects or images (the mean of each image, so groups with many cells from few images are not over-weighted).

**Multi-group analysis.** *Groups: mean / median* takes a *Control* group (each group vs the control, Welch t and Mann-Whitney U, Holm-adjusted, stars on the plot), a *Groups shown* picker, and a long-format table (one row per metric × group). *Groups: correlation* has an *All groups* mode (one heatmap per group and a Cochran's Q map of where the correlation differs between groups) besides *A vs B* (◀ ▶ step B through the groups). Per-group correlation matrices are cached, so switching is fast with many groups.

**Group-wise export.** After every batch (and with *Export by group…* in the Analysis tab) `<output>/groups/` gets one folder per group (`<group>_cells.csv`, `_mito.csv`, `_cells_green_pos/neg.csv`, `_per_image.csv`, `_correlations.csv`, `_results.xlsx`), `group_comparison.xlsx` (per metric × group: n, mean, SD, SEM, 95 % CI, median, IQR, the across-group tests and, with a control, each group vs control with Holm adjustment; all / green+ / green− cells, per cell, per image and per mito object) and `groups_overview.csv`. CLI: `--control NAME`, `--no-group-export`.

**Pixel size.** Read from OME `PhysicalSizeX`, the ImageJ `unit=` (micron, nm, mm, cm, inch) with `XResolution`, or `ResolutionUnit` (cm / inch); values outside 0.005–10 µm/px (e.g. 72 dpi) count as missing. The table's *µm/px* column shows it (red `?` if missing); *Run all* asks for the missing values once, and each cell can be edited. Recorded as `pixel_size_um` / `pixel_size_source` in the settings sheet. CLI: `--pixel-size-um` (else 0.099 with a warning).

**Manual thresholds.** Options → *Thresholds*: nuclei, cell area (mito density level, default 0.12), mitochondria, bright green (green+) and green puncta are automatic per image, or one manual value for every image of the batch. *Preview / adjust thresholds…* shows the mask live on any image set as you move a slider (with that image's automatic value), then *Apply to all images* or *This image only* (stored in the table's *Thresholds* column; it wins over the batch value). Values used and the automatic ones go to the settings sheet. CLI: `--thr-nuclei`, `--cell-fg-level`, `--thr-mito`, `--thr-green-bright`, `--thr-puncta`.

Correlation views (per sample and pooled): heatmap of Pearson r (or Spearman ρ) for every mito × green pair; click a square for its scatter with the least-squares line, 95 % band, r, R², p and n, split into quadrants at the mean (or median) of X and Y with the share of points in each, so negative relations (quadrants II / IV) stand out. With log X / log Y the fit uses log10 values.

Tabs: **Cells** (all per-cell metrics; column headers show `(+)` / `(−)`), **Mito objects** and **Green puncta** (sortable, filter by green+ / green− and by cell), **Correlation** (Spearman heatmap of all pairs within all / green+ / green− cells; click a square for the scatter plot, coloured by cell), **Green+ / − cells** (the green calls), and overlay images.

## Outputs
Folder layout: `<output>/<dataset>/<preset>-<dataset>/<sample>/`; with no output folder given, `<output>` is `dataset/` inside the common parent folder of the input folders. Each `<preset>-<dataset>` folder also gets the pooled tables of its samples: `<preset>-<dataset>_all_cells.csv`, `…_all_mito.csv`, their `_green_pos` / `_green_neg` versions, and `…_all_results.xlsx` (with a `sample` column).

Per sample (in its own folder, prefixed with the sample name):
- `_results.xlsx` – sheets `cells_green_pos`, `mito_green_pos`, `puncta_green_pos`, `corr_green_pos`, the same four for `green_neg`, then `all_*` (every analysed cell) and `settings`
- `_per_cell_green_pos.csv`, `_per_mito_green_pos.csv`, `_green_puncta_green_pos.csv`, `_correlations_green_pos.csv`, the same for `_green_neg`, and the all-cell `_per_cell.csv`, `_per_mito.csv`, `_green_puncta.csv`, `_correlations.csv`
- `_green_cells.png` – green+ (red) / green− (grey) cells, bright green (yellow)
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
./build_mac.sh --locked        # exact, validated versions (requirements-lock.txt) instead of the newest releases
./build_mac.sh --clean         # recreate the build environment (.venv)
./build_mac.sh --dmg           # also make "Mito Analyzer-<version>-<arch>.dmg" to give to others (Apple Silicon only, unsigned)
open "Mito Analyzer.app"
```
Requires Python 3.11+ (installed with Homebrew if missing and `brew` is available). `requirements.txt` has lower bounds only, so the newest releases are installed; `requirements-lock.txt` keeps the last validated versions. The app is unsigned and built locally; if macOS still blocks it, right-click → Open once.

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
