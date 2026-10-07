"""Per-cell MiNA (Mitochondrial Network Analysis) on the red mitochondria channel.

Python reproduction of MiNA_Analyze_Morphology.py with MiNA's GUI defaults
(no Median / Unsharp / CLAHE, threshold op = otsu, no ridge detection,
AnalyzeSkeleton without pruning). One deliberate difference from MiNA:
MiNA crops the ROI's bounding box and thresholds/skeletonizes everything in it,
so neighbouring cells leak in. Here the Otsu threshold is computed only from
pixels inside each cell mask, and only mask pixels can be foreground.

Steps per cell (mask from cell_roi.py --nuclei labels):
  1. Otsu threshold on red intensities inside the mask; foreground = value > t
  2. mitochondrial footprint = foreground area (um^2)
  3. Skeletonize (Lee thinning == ImageJ "Skeletonize (2D/3D)")
  4. AnalyzeSkeleton-like graph: end points (1 neighbour), slab (2),
     junctions (>=3; 8-connected junction pixels form one vertex);
     a vertex-free closed ring gets one artificial vertex and one self-loop edge
  5. branch lengths, summed branch length per network (incl. 0 for single-pixel
     networks), branches per network, donuts; population SD as in MiNA.

Usage: python3 mina_per_cell.py RED.tif LABELS.tif ROIS.csv -o OUTDIR [--green GREEN.tif]
"""
import argparse, csv, math, os
import numpy as np, tifffile
from scipy import ndimage as ndi
from skimage import filters, measure, morphology, segmentation
try:
    from . import imgio
except ImportError:  # run as a script
    import imgio

K8 = np.ones((3, 3), int)


def load_channel(path, ch):
    a = imgio.imread(path)
    if imgio.is_czi(imgio.split_ref(path)[0]):
        px_um = imgio.pixel_size(path)[0] or 1.0
    else:
        with tifffile.TiffFile(path) as t:
            xr = t.pages[0].tags.get('XResolution')
            px_um = xr.value[1] / xr.value[0] if xr else 1.0
    if a.ndim == 3 and a.shape[-1] == 3:
        img = a[..., ch].astype(np.uint8).copy()
        others = [c for c in range(3) if c != ch]
        bar = (a[..., others[0]] > 200) & (a[..., others[1]] > 200)  # white scale bar / text
        img[ndi.binary_dilation(bar, iterations=2)] = 0
    else:
        img = a
    return img, px_um


def analyze_skeleton(sk, px):
    """Return (edges, n_graphs, branches_per_graph, donuts, end_px, junc_px) like AnalyzeSkeleton_."""
    sk = sk.astype(bool)
    nb = ndi.convolve(sk.astype(int), K8, mode='constant') - 1
    nb[~sk] = 0
    graphs, n_graphs = ndi.label(sk, structure=K8)
    junc = sk & (nb >= 3)
    end = sk & (nb <= 1)
    slab = sk & (nb == 2)
    vlab = np.zeros(sk.shape, int)
    jl, nj = ndi.label(junc, structure=K8)
    vlab[junc] = jl[junc]
    ey, ex = np.nonzero(end)
    vlab[ey, ex] = nj + 1 + np.arange(len(ey))
    nv = nj + len(ey)
    H, W = sk.shape

    def nbrs(y, x):
        for dy in (-1, 0, 1):
            for dx in (-1, 0, 1):
                if (dy or dx) and 0 <= y + dy < H and 0 <= x + dx < W and sk[y + dy, x + dx]:
                    yield y + dy, x + dx

    edges = []  # (graph_id, v1, v2, length_um)
    # slab segments
    sl, ns = ndi.label(slab, structure=K8)
    objs = ndi.find_objects(sl)
    for i in range(ns):
        sy, sx = objs[i]
        pts = [(y + sy.start, x + sx.start) for y, x in zip(*np.nonzero(sl[sy, sx] == i + 1))]
        g = graphs[pts[0]]
        ptset = set(pts)
        # find a start pixel touching a vertex
        start, v_start = None, None
        for p in pts:
            vs = [q for q in nbrs(*p) if vlab[q]]
            if vs:
                start, v_start = p, vs[0]
                break
        if start is None:  # closed ring without vertices -> artificial vertex, self-loop
            nv += 1
            length = 0.0
            prev, cur = None, pts[0]
            for _ in range(len(pts)):
                nxt = [q for q in nbrs(*cur) if q in ptset and q != prev]
                if not nxt:
                    break
                length += math.dist(cur, nxt[0]); prev, cur = cur, nxt[0]
                if cur == pts[0]:
                    break
            edges.append((g, nv, nv, length * px))
            continue
        length = math.dist(v_start, start)
        prev, cur, seen = v_start, start, {start}
        while True:
            nxt = [q for q in nbrs(*cur) if q != prev and q not in seen]
            slab_n = [q for q in nxt if q in ptset]
            if slab_n:
                length += math.dist(cur, slab_n[0]); prev, cur = cur, slab_n[0]; seen.add(cur)
                continue
            vend = [q for q in nbrs(*cur) if vlab[q] and q != prev]
            if vend:
                length += math.dist(cur, vend[0])
            edges.append((g, vlab[v_start], vlab[vend[0]] if vend else vlab[v_start], length * px))
            break
    # vertex-vertex edges with no slab in between
    vy, vx = np.nonzero(vlab)
    done = set()
    for y, x in zip(vy, vx):
        for q in nbrs(y, x):
            a, b = vlab[y, x], vlab[q]
            if b and a != b:
                key = (min(a, b), max(a, b))
                if key not in done:
                    done.add(key)
                    edges.append((graphs[y, x], a, b, math.dist((y, x), q) * px))
    # per-graph stats
    branches = np.zeros(n_graphs + 1, int)
    summed = np.zeros(n_graphs + 1)
    vcount = [dict() for _ in range(n_graphs + 1)]
    for g, a, b, L in edges:
        branches[g] += 1; summed[g] += L
        for v in (a, b):
            vcount[g][v] = vcount[g].get(v, 0) + 1
    donuts = sum(1 for g in range(1, n_graphs + 1)
                 if branches[g] >= 1 and all(c >= 2 for c in vcount[g].values()))
    return edges, n_graphs, branches[1:], summed[1:], donuts, end, junc


def pstats(v):
    v = np.asarray(v, float)
    if v.size == 0:
        return (float('nan'),) * 3
    return float(v.mean()), float(np.median(v)), float(v.std(ddof=0))


def analyze_cells(red, lab, rois, px, green=None, all_cells=False, log=print, mask_fn=None):
    """MiNA per cell. `rois` are cell_roi rows (roi, status, ...). Returns (rows, binary, skeleton, ends, junctions).
    By default the mito mask is MiNA's Otsu threshold inside the cell; `mask_fn(cell_mask) -> (mask, extra_cols)`
    replaces it (e.g. separated mito objects)."""
    full_bin = np.zeros(red.shape, bool); full_sk = np.zeros(red.shape, bool)
    full_end = np.zeros(red.shape, bool); full_junc = np.zeros(red.shape, bool)
    rows = []
    for r in rois:
        if r.get('status', 'ok') != 'ok' and not all_cells:
            continue
        lid = int(r['roi'].replace('cell', ''))
        m = lab == lid
        vals = red[m]
        if vals.size == 0 or vals.min() == vals.max():
            log(f"{r['roi']}: skipped (no intensity variation in red channel)")
            continue
        extra = {}
        if mask_fn is None:
            t = filters.threshold_otsu(vals)
            b = m & (red > t)
        else:
            t = float('nan')
            b, extra = mask_fn(m)
        if not b.any():
            log(f"{r['roi']}: skipped (no mitochondria found)")
            continue
        sk = morphology.skeletonize(b, method='lee').astype(bool)
        edges, ng, br, summed, donuts, end, junc = analyze_skeleton(sk, px)
        full_bin |= b; full_sk |= sk; full_end |= end; full_junc |= junc
        bl = pstats([e[3] for e in edges]); sl = pstats(summed); nb = pstats(br)
        cell_area = m.sum() * px * px; fp = b.sum() * px * px
        row = dict(cell=r['roi'], status=r.get('status', 'ok'), cell_area_um2=cell_area, otsu_threshold=float(t),
                   mitochondrial_footprint_um2=fp, footprint_fraction=fp / cell_area,
                   red_mean_in_mito=float(red[b].mean()),
                   n_networks=ng, n_branches=len(edges),
                   total_branch_length_um=float(sum(e[3] for e in edges)),
                   branch_length_mean_um=bl[0], branch_length_median_um=bl[1], branch_length_stdev_um=bl[2],
                   summed_branch_lengths_mean_um=sl[0], summed_branch_lengths_median_um=sl[1],
                   summed_branch_lengths_stdev_um=sl[2],
                   network_branches_mean=nb[0], network_branches_median=nb[1], network_branches_stdev=nb[2],
                   donuts=donuts, **extra)
        if green is not None:
            row['green_poi_mean_cell'] = float(green[m].mean())
        rows.append(row)
        log(f"{r['roi']}: footprint {fp:.1f} um2, {ng} networks, {len(edges)} branches, {donuts} donuts")
    return rows, full_bin, full_sk, full_end, full_junc


def write_outputs(red, lab, rois, rows, b, sk, end, junc, outdir, prefix):
    if rows:
        with open(os.path.join(outdir, f'{prefix}_mina_per_cell.csv'), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            for row in rows:
                w.writerow({k: (f'{v:.4f}' if isinstance(v, float) else v) for k, v in row.items()})
    tifffile.imwrite(os.path.join(outdir, f'{prefix}_mito_binary.tif'), b.astype(np.uint8) * 255)
    tifffile.imwrite(os.path.join(outdir, f'{prefix}_mito_skeleton.tif'), sk.astype(np.uint8) * 255)
    render(red, lab, rois, rows, b, sk, end, junc, outdir, prefix)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('red'); ap.add_argument('labels'); ap.add_argument('rois_csv')
    ap.add_argument('-o', '--outdir', default='mina_out')
    ap.add_argument('--green', help='protein-of-interest channel; adds mean intensity per cell')
    ap.add_argument('--all-cells', action='store_true', help='also analyse cells with status != ok')
    A = ap.parse_args()
    os.makedirs(A.outdir, exist_ok=True)
    red, px = load_channel(A.red, 0)
    green = load_channel(A.green, 1)[0] if A.green else None
    lab = tifffile.imread(A.labels)
    rois = list(csv.DictReader(open(A.rois_csv)))
    res = analyze_cells(red, lab, rois, px, green, A.all_cells)
    prefix = os.path.splitext(os.path.basename(A.red))[0]
    write_outputs(red, lab, rois, *res, A.outdir, prefix)


def render(red, lab, rois, rows, b, sk, end, junc, outdir, prefix):
    from matplotlib.figure import Figure  # no pyplot: safe to call from a worker thread
    ok = {r['cell'] for r in rows}
    hi = np.percentile(red, 99.5)
    base = np.clip(red / hi, 0, 1)
    rgb = np.dstack([base] * 3) * 0.85
    rgb[b] = rgb[b] * 0.6 + np.array([0.9, 0.1, 0.9]) * 0.4                 # footprint: magenta tint (MiNA)
    for r in rois:                                                         # cell outlines
        lid = int(r['roi'].replace('cell', ''))
        ed = segmentation.find_boundaries(lab == lid, mode='inner')
        ed = ndi.binary_dilation(ed, iterations=1)
        rgb[ed] = (0.0, 0.9, 1.0) if r['roi'] in ok else (0.45, 0.45, 0.45)
    rgb[sk] = (0.1, 1.0, 0.1)                                              # skeleton: green (MiNA)
    for mask, col in ((end, (1.0, 1.0, 0.0)), (junc, (0.2, 0.4, 1.0))):    # end: yellow, junction: blue
        rgb[ndi.binary_dilation(mask, iterations=1)] = col
    fig = Figure(figsize=(12, 12), dpi=150); ax = fig.subplots()
    ax.imshow(rgb); ax.set_axis_off()
    for r in rois:
        ax.text(float(r['centroid_x']), float(r['centroid_y']),
                r['roi'].replace('cell', '') + ('' if r['roi'] in ok else '\n' + r.get('status', 'excluded')),
                color='white' if r['roi'] in ok else '#aaaaaa', fontsize=13 if r['roi'] in ok else 9,
                ha='center', va='center', weight='bold',
                bbox=dict(facecolor='black', alpha=0.5, pad=1.5, edgecolor='none'))
    ax.set_title('MiNA per cell (red = mito). cyan: analysed cell ROI, grey: excluded; '
                 'magenta: footprint, green: skeleton, yellow: ends, blue: junctions', fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f'{prefix}_mina_overlay.png'))
    # per-cell zoom panels (up to 6 cells per row; red on top, overlay below)
    n = len(rows)
    if not n:
        return
    nc = min(n, 6); nr = -(-n // nc)
    fig = Figure(figsize=(4.5 * nc, 9 * nr), dpi=150); axes = fig.subplots(2 * nr, nc, squeeze=False)
    for a in axes.ravel():
        a.set_axis_off()
    for j, row in enumerate(rows):
        r = next(x for x in rois if x['roi'] == row['cell'])
        x0, y0 = int(r['bbox_x']), int(r['bbox_y']); w, h = int(r['bbox_w']), int(r['bbox_h'])
        sl = (slice(y0, y0 + h), slice(x0, x0 + w))
        top, bot = axes[2 * (j // nc), j % nc], axes[2 * (j // nc) + 1, j % nc]
        top.imshow(np.clip(red[sl] / hi, 0, 1), cmap='gray')
        top.contour(lab[sl] == int(row['cell'][4:]), [0.5], colors='cyan', linewidths=1)
        t = row['otsu_threshold']
        top.set_title(f"{row['cell']} red" + (f" (Otsu t={t:.0f})" if t == t else ''), fontsize=10)
        bot.imshow(rgb[sl]); bot.set_title(
            f"networks {row['n_networks']}, branches {row['n_branches']}, donuts {row['donuts']}", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(outdir, f'{prefix}_mina_cells.png'))


if __name__ == '__main__':
    main()
