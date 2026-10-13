"""Cell / nucleus / cytosol ROIs as ImageJ RoiSets, in the layout Cytosol Viewer (cytosol-analysis) writes.

Per cell: `cellNNN` (whole-cell polygon), `nucNNN` (nucleus polygon) and `cytoNNN` (cell minus nucleus,
a composite ROI whose holes are the nucleus). `read_roiset` turns such a file (or loose `.roi` files) back into
label images for the analysis; `write_used_rois` writes the cell and mito ROIs an analysis used, per image.
"""
import os
import re
import zipfile

import numpy as np
from skimage import draw, measure

ROI_EXT = ('.zip', '.roi')
_ROISET_TAIL = re.compile(r'([_\- ]?(cell|cyto|cytosol)?[_\- ]?(rois?|roiset)(_filtered)?)$', re.I)


def is_roi_file(path):
    return str(path).lower().endswith(ROI_EXT)


def roi_key(path):
    """'New-01_RoiSet.zip' -> 'new-01': the image name a Cytosol Viewer RoiSet belongs to."""
    stem = os.path.splitext(os.path.basename(str(path)))[0]
    return _ROISET_TAIL.sub('', stem).strip().lower()


def _keys(job):
    """Names an image set may be known by: sample name, CZI / TIFF file stems (with and without channel tags)."""
    from . import imgio, pipeline
    out = {str(job.get('name') or '').lower()}
    for k in ('czi', 'red', 'green', 'blue'):
        p = job.get(k)
        if not p:
            continue
        f = imgio.split_ref(p)[0]
        stem = os.path.splitext(os.path.basename(f))[0]
        out |= {stem.lower(), pipeline.sample_name(f).lower()}
    return out - {''}


def match_roi_files(files, jobs):
    """{job index: roi file} by name, and the files that matched no job."""
    by_key = {}
    for f in files:
        by_key.setdefault(roi_key(f), []).append(f)
    found, used = {}, set()
    for i, j in enumerate(jobs):
        keys = _keys(j)
        hit = next((by_key[k][0] for k in keys if k in by_key), None)
        if hit is None:   # 'Processed_New01' vs 'New01' and the like: one name ends with the other
            hit = next((fs[0] for k, fs in by_key.items() if k and any(
                n.endswith(k) or k.endswith(n) for n in keys if len(n) >= 3 and len(k) >= 3)), None)
        if hit is not None:
            found[i] = hit
            used.add(hit)
    return found, [f for f in files if f not in used]


def collect_roi_files(paths):
    """RoiSet .zip / .roi files among `paths` (folders are searched one level deep)."""
    out = []
    for p in paths:
        if os.path.isdir(p):
            out += sorted(os.path.join(p, f) for f in os.listdir(p) if is_roi_file(f) and not f.startswith('.'))
        elif is_roi_file(p):
            out.append(p)
    return out


# ---- reading ---------------------------------------------------------------

def _polygon_mask(xy, shape):
    """Pixel (i, j) is inside when its centre is, as in ImageJ."""
    xy = np.asarray(xy, float) - 0.5
    rr, cc = draw.polygon(xy[:, 1], xy[:, 0], shape)
    m = np.zeros(shape, bool)
    m[rr, cc] = True
    return m


def _area(xy):
    x, y = xy[:, 0], xy[:, 1]
    return abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1))) / 2


def _roi_polygons(r):
    """Polygons of an ImageJ ROI in image coordinates, largest (outer boundary) first."""
    import roifile
    if r.roitype == roifile.ROI_TYPE.RECT and not r.composite:
        x0, y0, x1, y1 = r.left, r.top, r.right, r.bottom
        return [np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]], float)]
    if r.roitype == roifile.ROI_TYPE.OVAL and not r.composite:
        cx, cy = (r.left + r.right) / 2, (r.top + r.bottom) / 2
        rx, ry = (r.right - r.left) / 2, (r.bottom - r.top) / 2
        t = np.linspace(0, 2 * np.pi, 64, endpoint=False)
        return [np.c_[cx + rx * np.cos(t), cy + ry * np.sin(t)]]
    c = r.coordinates(multi=True) if r.composite else [r.coordinates()]
    polys = [np.asarray(p, float) for p in c if p is not None and len(p) >= 3]
    polys.sort(key=lambda p: -_area(p))
    return polys


def _read_rois(path):
    import roifile
    if str(path).lower().endswith('.zip') or zipfile.is_zipfile(path):
        r = roifile.roiread(path)
    else:
        r = roifile.ImagejRoi.fromfile(path)
    return r if isinstance(r, list) else [r]


def read_roiset(path, shape, detected_nuclei=None):
    """Cell ROIs of one image from a Cytosol Viewer RoiSet (.zip) or .roi file.

    `cell*` (and any other unnamed polygon) -> cells, `nuc*` -> nuclei (matched by number, else by overlap),
    `cyto*` -> cells only when there is no cell ROI of that number (its holes are then the nucleus).
    Cells without a nucleus ROI take the nuclei found in the blue channel inside them (`detected_nuclei`).
    Returns dict(labels, nuclei, seed) as used for reviewed ROIs; seed 'imported' for every cell."""
    rois = _read_rois(path)
    labels = np.zeros(shape, np.int32)
    nuclei = np.zeros(shape, np.int32)
    cells, nucs, cytos = [], [], []
    for r in rois:
        name = (r.name or '').lower()
        (nucs if name.startswith('nuc') else cytos if name.startswith('cyto') else cells).append(r)

    def num(r):
        d = ''.join(ch for ch in (r.name or '') if ch.isdigit())
        return int(d) if d else None

    cell_nums = {num(x) for x in cells}
    used = set()
    for r in cells + [c for c in cytos if num(c) is None or num(c) not in cell_nums]:
        parts = _roi_polygons(r)
        if not parts:
            continue
        m = _polygon_mask(parts[0], shape) & (labels == 0)
        if not m.any():
            continue
        l = num(r) if num(r) and num(r) not in used else max(used | {0}) + 1
        used.add(l)
        labels[m] = l
        if r in cytos:
            for hole in parts[1:]:   # the holes of a cytosol ROI are the nucleus
                nuclei[_polygon_mask(hole, shape) & (labels == l)] = l
    for r in nucs:
        parts = _roi_polygons(r)
        if not parts:
            continue
        m = _polygon_mask(parts[0], shape)
        l = num(r)
        if l not in used:
            ov = labels[m]
            ov = ov[ov > 0]
            l = int(np.bincount(ov).argmax()) if len(ov) else None
        if l:
            nuclei[m & (labels == l)] = l
    if not used:
        raise ValueError(f'No cell ROIs inside the {shape[1]}x{shape[0]} image in {os.path.basename(path)}')
    if detected_nuclei is not None:
        has = set(np.unique(nuclei)) - {0}
        fill = (np.asarray(detected_nuclei) > 0) & (labels > 0) & ~np.isin(labels, list(has)) & (nuclei == 0)
        nuclei[fill] = labels[fill]
    return dict(labels=labels, nuclei=nuclei, seed={l: 'imported' for l in used})


# ---- writing ---------------------------------------------------------------

def _contour(mask, tolerance=0.5):
    """Outer boundary of a mask as an (N, 2) x, y polygon on the pixel edges (ImageJ fills it back exactly)."""
    if not mask.any():
        return None
    ys, xs = np.nonzero(mask)
    y0, x0 = ys.min(), xs.min()
    sub = np.pad(mask[y0:ys.max() + 1, x0:xs.max() + 1], 1).astype(float)
    cs = measure.find_contours(sub, 0.5)
    if not cs:
        return None
    c = max(cs, key=len)
    if tolerance:
        c = measure.approximate_polygon(c, tolerance)
    xy = np.c_[c[:, 1] + x0 - 1 + 0.5, c[:, 0] + y0 - 1 + 0.5]
    if len(xy) > 1 and np.allclose(xy[0], xy[-1]):
        xy = xy[:-1]
    return xy if len(xy) >= 3 else None


def _polygon_roi(xy, name):
    import roifile
    roi = roifile.ImagejRoi.frompoints(np.asarray(xy, np.float32), name=name)
    roi.roitype = roifile.ROI_TYPE.POLYGON
    return roi


def _composite_roi(polys, name):
    """ImageJ ShapeRoi (cell with the nucleus as a hole). Path ops: 0 MOVETO, 1 LINETO, 4 CLOSE."""
    import roifile
    path = []
    for xy in polys:
        xy = np.asarray(xy, np.float32)
        path += [0, float(xy[0, 0]), float(xy[0, 1])]
        for x, y in xy[1:]:
            path += [1, float(x), float(y)]
        path.append(4)
    allxy = np.concatenate([np.asarray(p) for p in polys])
    roi = roifile.ImagejRoi()
    roi.roitype = roifile.ROI_TYPE.RECT
    roi.version = 228
    roi.name = name
    roi.left, roi.top = int(np.floor(allxy[:, 0].min())), int(np.floor(allxy[:, 1].min()))
    roi.right, roi.bottom = int(np.ceil(allxy[:, 0].max())), int(np.ceil(allxy[:, 1].max()))
    roi.multi_coordinates = np.asarray(path, np.float32)
    roi.shape_roi_size = len(path)
    return roi


def _write(path, rois):
    import roifile
    if os.path.exists(path):
        os.remove(path)
    if rois:
        roifile.roiwrite(path, rois)


def cell_rois(lab, nuc, cells):
    """cellNNN / nucNNN / cytoNNN ImageJ ROIs of the given cell labels."""
    out = []
    for l in cells:
        cm = lab == l
        c = _contour(cm)
        if c is None:
            continue
        out.append(_polygon_roi(c, f'cell{l:03d}'))
        nm = (nuc == l) & cm
        if nm.any():
            nl = measure.label(nm)
            n = _contour(nl == np.bincount(nl.ravel())[1:].argmax() + 1)
            if n is not None:
                out.append(_polygon_roi(n, f'nuc{l:03d}'))
                out.append(_composite_roi([c, n[::-1]], f'cyto{l:03d}'))
                continue
        out.append(_polygon_roi(c, f'cyto{l:03d}'))
    return out


def mito_rois(mito, lab, cells):
    """One polygon per mito object inside each cell, named cellNNN_mNNN."""
    out = []
    for l in cells:
        ml = measure.label(mito & (lab == l))
        for r in measure.regionprops(ml):
            c = _contour(np.pad(r.image, 0))
            if c is None:
                continue
            c = c + np.array([r.bbox[1], r.bbox[0]])
            out.append(_polygon_roi(c, f'cell{l:03d}_m{r.label:03d}'))
    return out


def write_used_rois(outdir, name, lab, nuc, mito, cells, red=None):
    """Cell and mito ROIs one analysis used, in <outdir>/<name>_used_ROIs/:
    <name>_cell_RoiSet.zip (cellNNN / nucNNN / cytoNNN, opens in ImageJ and Cytosol Viewer and can be loaded
    back as cell ROIs), <name>_mito_RoiSet.zip (cellNNN_mNNN), label TIFFs for cells / cytosol / nuclei / mito,
    and a PNG overlay. Returns the folder."""
    import tifffile
    from skimage import io, segmentation
    d = os.path.join(outdir, f'{name}_used_ROIs')
    os.makedirs(d, exist_ok=True)
    cells = sorted(int(c) for c in cells)
    keep = np.isin(lab, cells)
    clab = np.where(keep, lab, 0).astype(np.uint16)
    nlab = np.where(keep & (nuc == lab), nuc, 0).astype(np.uint16)
    mlab = np.where(keep & mito, lab, 0).astype(np.uint16)
    _write(os.path.join(d, f'{name}_cell_RoiSet.zip'), cell_rois(lab, nuc, cells))
    _write(os.path.join(d, f'{name}_mito_RoiSet.zip'), mito_rois(mito, lab, cells))
    for tag, a in (('cell', clab), ('cytosol', np.where(nlab > 0, 0, clab).astype(np.uint16)),
                   ('nuclei', nlab), ('mito', mlab)):
        tifffile.imwrite(os.path.join(d, f'{name}_{tag}_labels.tif'), a, imagej=True,
                         metadata={'Info': f'{tag} pixels labelled with their cell number (0 = not used)'})
    if red is not None:
        r = red.astype(float)
        r = np.clip(r / max(np.percentile(r, 99.5), 1e-9), 0, 1) ** 0.7
        ov = np.dstack([r, r * 0.25, r * 0.25])
        ov[segmentation.find_boundaries(mlab, mode='inner')] = (1, 0.2, 1)
        ov[segmentation.find_boundaries(nlab, mode='inner')] = (0.3, 0.8, 1)
        ov[segmentation.find_boundaries(clab, mode='inner')] = (1, 1, 1)
        io.imsave(os.path.join(d, f'{name}_used_rois.png'), (ov * 255).astype(np.uint8), check_contrast=False)
    return d
