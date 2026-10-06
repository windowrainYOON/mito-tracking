"""Entry point. With no arguments it opens the GUI; with `--cli` it runs headless:

    python -m mito_app --cli RED.tif GREEN.tif BLUE.tif -o OUTDIR [--exclude-binucleate] [--include-edge-cells]
"""
import argparse, sys, warnings

# scikit-image deprecation notices (binary_opening etc.) are noise for end users
warnings.filterwarnings('ignore', category=FutureWarning, module='mito_app')


def cli(argv):
    from . import pipeline
    ap = argparse.ArgumentParser(prog='mito_app --cli')
    ap.add_argument('red'); ap.add_argument('green'); ap.add_argument('blue')
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
    A = ap.parse_args(argv)
    p = pipeline.Params(A.min_area_um2, A.binuc_tau, A.exclude_binucleate, A.include_edge_cells,
                        A.pixel_size_um, A.puncta_sensitivity, A.min_mito_area_um2, A.mito_method)
    pipeline.run(A.red, A.green, A.blue, A.outdir, p, A.name)
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ['--cli']:
        return cli(argv[1:])
    from .gui import main as gui_main
    return gui_main()
