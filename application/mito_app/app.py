"""Entry point. With no arguments it opens the GUI; with `--cli` it runs headless:

    python -m mito_app --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR [--exclude-binucleate] [--include-edge-cells]
    python -m mito_app --cli --batch FOLDER_OR_TIFF [...] -o OUTDIR [--dataset NAME] [--preset NAME] [--group NAME]
    python -m mito_app --cli IMAGE.czi -o OUTDIR [--czi-channels mito=0,protein=1,nucleus=2]
    ... --rois ROISET.zip | FOLDER     cell ROIs from Cytosol Viewer instead of the detected ones

Batch mode groups the TIFFs (folders are searched recursively) into red/green/blue sets by file name and writes each set to
OUTDIR/<dataset>/<preset>-<dataset>/<sample>/ (dataset defaults to the name of the folder given,
preset to the mito method). The group of each set (default: its dataset) goes to OUTDIR/groups.csv, and the
group-wise export (one folder per group + group_comparison.xlsx) to OUTDIR/groups/.
A Zeiss .czi file is one set by itself; --czi-channels gives the channel (0-based) of each role for every CZI
file of the run (default: guessed from the channel names, e.g. DAPI = nucleus, longest emission = mito).

--rois takes one RoiSet .zip / .roi (single set) or files / folders of them (batch); each is matched to its image set
by name ('New-01_RoiSet.zip' -> New-01). Sets without a matching file use the detected cell ROIs.

Thresholds are automatic per image unless a manual value is given (one value for every image):
    --thr-nuclei, --cell-fg-level, --thr-mito, --thr-green-bright, --thr-puncta
Images without a pixel size in the TIFF use --pixel-size-um (or, if 0, the 0.099 um reference, with a warning).
"""
import argparse, os, sys, warnings

# scikit-image deprecation notices (binary_opening etc.) are noise for end users
warnings.filterwarnings('ignore', category=FutureWarning, module='mito_app')


def cli(argv):
    from . import pipeline
    ap = argparse.ArgumentParser(prog='mito_app --cli')
    ap.add_argument('inputs', nargs='+', help='RED GREEN BLUE, one .czi file, or with --batch any TIFFs / CZIs / folders')
    ap.add_argument('--czi-channels', default='',
                    help='CZI channel of each role for all CZI files, e.g. mito=0,protein=1,nucleus=2 '
                         '(default: guessed from the channel names)')
    ap.add_argument('--rois', nargs='+', default=[],
                    help='cell ROIs from Cytosol Viewer: RoiSet .zip / .roi files or folders, matched to the image sets '
                         'by name (one file and one set: used as is)')
    ap.add_argument('--batch', action='store_true')
    ap.add_argument('--dataset'); ap.add_argument('--preset')
    ap.add_argument('--group', help='group of all sets in this batch (for the group comparison; default: the dataset)')
    ap.add_argument('-o', '--outdir', required=True)
    ap.add_argument('--name')
    ap.add_argument('--min-area-um2', type=float, default=250)
    ap.add_argument('--binuc-tau', type=float, default=0.0)
    ap.add_argument('--exclude-binucleate', action='store_true')
    ap.add_argument('--include-edge-cells', action='store_true')
    ap.add_argument('--pixel-size-um', type=float, default=0.0, help='0 = read from the TIFF')
    ap.add_argument('--puncta-sensitivity', type=float, default=1.0)
    ap.add_argument('--min-mito-area-um2', type=float, default=0.05)
    ap.add_argument('--mito-method', choices=('split', 'otsu'), default='split')
    ap.add_argument('--no-dim-nuclei', action='store_true', help='only in-focus nuclei seed a cell')
    ap.add_argument('--no-edge-trim', action='store_true',
                    help='do not keep frame-touching cells by trimming them (cytosol: to 1 um inside the frame; '
                         'mito: cutting their tip off along a mito-free line)')
    ap.add_argument('--cell-method', choices=('cytosol', 'mito'), default='cytosol',
                    help='cell ROI detection: cytosol = Cytosol Analysis / Cytosol Viewer algorithm '
                         '(mito_app/cytosol_cells.py), mito = the earlier detection of this app')
    ap.add_argument('--cell-channel', choices=('red', 'green'), default='red',
                    help='cytosol detection: channel whose dim signal outlines the cells')
    ap.add_argument('--nuclei-method', choices=('texture_merge', 'texture', 'otsu'), default='texture_merge',
                    help='nucleus detection of --cell-method mito (see mito_app/nuclei.py)')
    ap.add_argument('--green-pos-percent', type=float, default=2.0,
                    help='green-positive cell: bright green covers at least this %% of the cytoplasm')
    ap.add_argument('--thr-nuclei', type=float, default=0.0,
                    help='manual nucleus threshold on the smoothed blue (cytosol, otsu) / flattened blue (texture); '
                         '0 = auto')
    ap.add_argument('--thr-cell', type=float, default=0.0,
                    help='cytosol detection: manual cell threshold in cell-channel units; 0 = Otsu')
    ap.add_argument('--cell-fg-level', type=float, default=0.12,
                    help='--cell-method mito: cell area = mito density (red at 3 um / its 99th percentile) above this')
    ap.add_argument('--thr-mito', type=float, default=0.0,
                    help='manual mito threshold on the preprocessed red (rolling ball 1.5 um, sigma 0.7 px); 0 = auto')
    ap.add_argument('--thr-green-bright', type=float, default=0.0,
                    help='manual bright-green threshold for green+ cells (auto-levelled green); 0 = auto')
    ap.add_argument('--thr-puncta', type=float, default=0.0,
                    help='manual green puncta threshold (top-hat units); 0 = auto')
    ap.add_argument('--control', help='control group for the group comparison export (each group vs control)')
    ap.add_argument('--no-group-export', action='store_true', help='do not write OUTDIR/groups/')
    A = ap.parse_args(argv)
    p = pipeline.Params(min_area_um2=A.min_area_um2, binuc_tau=A.binuc_tau, exclude_binucleate=A.exclude_binucleate,
                        include_edge_cells=A.include_edge_cells, pixel_size_um=A.pixel_size_um,
                        puncta_sensitivity=A.puncta_sensitivity, min_mito_area_um2=A.min_mito_area_um2,
                        mito_method=A.mito_method, green_pos_percent=A.green_pos_percent,
                        dim_nuclei=not A.no_dim_nuclei, trim_edge_cells=not A.no_edge_trim,
                        nuclei_method=A.nuclei_method, thr_nuclei=A.thr_nuclei, cell_fg_level=A.cell_fg_level,
                        cell_method=A.cell_method, cell_channel=A.cell_channel, thr_cell=A.thr_cell,
                        thr_mito=A.thr_mito, thr_green_bright=A.thr_green_bright, thr_puncta=A.thr_puncta)
    try:
        roles = pipeline.imgio.parse_roles(A.czi_channels) if A.czi_channels else None
    except ValueError as e:
        ap.error(str(e))
    roi_files = pipeline.roi_io.collect_roi_files(A.rois)
    if A.rois and not roi_files:
        ap.error('--rois: no .zip / .roi files found')

    def one_roi(s):
        if not roi_files:
            return None
        if len(roi_files) == 1:
            return roi_files[0]
        m, _ = pipeline.roi_io.match_roi_files(roi_files, [s])
        if 0 not in m:
            ap.error(f"--rois: no file matches {s['name']}")
        return m[0]

    if not A.batch and len(A.inputs) == 1 and pipeline.imgio.is_multi(A.inputs[0]):
        sets, _ = pipeline.find_sets(A.inputs, roles)
        if not sets:
            ap.error(f'{A.inputs[0]}: cannot tell the channel roles; give --czi-channels')
        s = sets[0]
        print('channels: ' + ', '.join(f'{lab} = {pipeline.imgio.display_name(s[k])}' for k, lab in pipeline.imgio.ROLES))
        pipeline.run(s['red'], s['green'], s['blue'], A.outdir, p, A.name or s['name'], roi_file=one_roi(s))
        return 0
    if not A.batch:
        if len(A.inputs) != 3:
            ap.error('give RED GREEN BLUE, one .czi file, or use --batch')
        s = dict(red=A.inputs[0], green=A.inputs[1], blue=A.inputs[2], name=A.name or pipeline.sample_name(A.inputs[0]))
        pipeline.run(*A.inputs, A.outdir, p, A.name, roi_file=one_roi(s))
        return 0
    sets, leftover = pipeline.find_sets(A.inputs, roles)
    for f in leftover:
        print(f'skipped (no complete red/green/blue set{"; give --czi-channels" if pipeline.imgio.is_czi(f) else ""}): {f}')
    for s in sets:
        if 'czi' in s:
            print(f"{s['name']}: " + ', '.join(f'{lab} = {pipeline.imgio.display_name(s[k])}' for k, lab in pipeline.imgio.ROLES))
    if not A.pixel_size_um:
        nopx = [s['name'] for s in sets if not pipeline.cell_roi.read_pixel_size(s['green'])[0]]
        if nopx:
            print(f'WARNING: {len(nopx)} image set(s) have no pixel size in the TIFF ({", ".join(nopx[:5])}'
                  f'{" …" if len(nopx) > 5 else ""}); they use {pipeline.cell_roi.REF_PX_UM} um/px. '
                  'Give --pixel-size-um if that is wrong.')
    rmap, unmatched = pipeline.roi_io.match_roi_files(roi_files, sets)
    if roi_files:
        print(f'cell ROIs from Cytosol Viewer for {len(rmap)} of {len(sets)} image set(s)'
              + (f'; not matched: {", ".join(os.path.basename(f) for f in unmatched)}' if unmatched else ''))
    results, failed = [], 0
    entries = [dict(outdir=pipeline.job_outdir(A.outdir, A.dataset or s['dataset'], A.preset or A.mito_method, s['name']),
                    dataset=A.dataset or s['dataset'], sample=s['name'], group=A.group or A.dataset or s['dataset'])
               for s in sets]
    if entries:
        pipeline.write_groups(A.outdir, entries)
    for n, (s, e) in enumerate(zip(sets, entries), 1):
        name = s['name']
        out = e['outdir']
        print(f'[{n}/{len(sets)}] {name} -> {out}')
        try:
            results.append(pipeline.run(s['red'], s['green'], s['blue'], out, p, name, roi_file=rmap.get(n - 1)))
        except Exception as e:
            failed += 1
            print(f'     FAILED: {type(e).__name__}: {e}')
    pipeline.write_group_tables(results)
    if results and not A.no_group_export:
        pipeline.write_group_exports(A.outdir, control=A.control)
    print(f'{len(results)} done, {failed} failed')
    return 1 if failed or not sets else 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['--cli']:
        return cli(argv[1:])
    from .gui import main as gui_main
    return gui_main()
