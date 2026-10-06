"""Entry point. With no arguments it opens the GUI; with `--cli` it runs headless:

    python -m mito_app --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR [--exclude-binucleate] [--include-edge-cells]
    python -m mito_app --cli --batch FOLDER_OR_TIFF [...] -o OUTDIR [--dataset NAME] [--preset NAME]

Batch mode groups the TIFFs (folders are searched recursively) into red/green/blue sets by file name and writes each set to
OUTDIR/<dataset>/<preset>-<dataset>/<sample>/ (dataset defaults to the name of the folder given,
preset to the mito method).
"""
import argparse, os, sys, warnings

# scikit-image deprecation notices (binary_opening etc.) are noise for end users
warnings.filterwarnings('ignore', category=FutureWarning, module='mito_app')


def cli(argv):
    from . import pipeline
    ap = argparse.ArgumentParser(prog='mito_app --cli')
    ap.add_argument('inputs', nargs='+', help='RED GREEN BLUE, or with --batch any TIFFs / folders')
    ap.add_argument('--batch', action='store_true')
    ap.add_argument('--dataset'); ap.add_argument('--preset')
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
    ap.add_argument('--no-green-split', action='store_true', help='do not redraw green+ cells by green intensity')
    ap.add_argument('--no-green-refine', action='store_true', help='do not move green patches between ROIs')
    ap.add_argument('--green-pos-percent', type=float, default=2.0,
                    help='green-positive cell: bright green covers at least this %% of the cytoplasm')
    A = ap.parse_args(argv)
    p = pipeline.Params(A.min_area_um2, A.binuc_tau, A.exclude_binucleate, A.include_edge_cells,
                        A.pixel_size_um, A.puncta_sensitivity, A.min_mito_area_um2, A.mito_method,
                        not A.no_green_refine, A.green_pos_percent, not A.no_green_split)
    if not A.batch:
        if len(A.inputs) != 3:
            ap.error('give RED GREEN BLUE, or use --batch')
        pipeline.run(*A.inputs, A.outdir, p, A.name)
        return 0
    sets, leftover = pipeline.find_sets(A.inputs)
    for f in leftover:
        print(f'skipped (no complete red/green/blue set): {f}')
    results, failed = [], 0
    for n, s in enumerate(sets, 1):
        name = s['name']
        dataset = A.dataset or s['dataset']
        out = pipeline.job_outdir(A.outdir, dataset, A.preset or A.mito_method, name)
        print(f'[{n}/{len(sets)}] {name} -> {out}')
        try:
            results.append(pipeline.run(s['red'], s['green'], s['blue'], out, p, name))
        except Exception as e:
            failed += 1
            print(f'     FAILED: {type(e).__name__}: {e}')
    pipeline.write_group_tables(results)
    print(f'{len(results)} done, {failed} failed')
    return 1 if failed or not sets else 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['--cli']:
        return cli(argv[1:])
    from .gui import main as gui_main
    return gui_main()
