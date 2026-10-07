"""Threshold preview: the image each threshold is applied to, its automatic value and the resulting mask,
for one image set. Used by the GUI's threshold dialog; the pipeline applies the same formulas
(nuclei.threshold_image, cell_roi.cell_density, mito_objects.preprocess, green_cells.bright_image,
green.tophat), so a value tuned here gives the same mask in the analysis.

Every manual threshold is a value of that image (not a factor), so one value applies to every image of a
batch; the automatic value of the image being previewed is shown next to it for reference."""
import numpy as np
from skimage import morphology

from . import cell_roi, green as green_q, green_cells, mina, mito_objects, nuclei

KEYS = ('thr_nuclei', 'cell_fg_level', 'thr_mito', 'thr_green_bright', 'thr_puncta')


class Preview:
    """Lazily computed threshold images of one image set (each is computed once, then reused)."""

    def __init__(self, red_path, green_path, blue_path, pixel_size_um=0.0, nuclei_method='texture_merge',
                 sensitivity=1.0):
        a, g, px = cell_roi.load_green(green_path)
        self.px = pixel_size_um if pixel_size_um > 0 else (px or cell_roi.REF_PX_UM)
        self.k = float(np.clip(cell_roi.REF_PX_UM / self.px, 0.2, 5))
        self.gn, self.norm = cell_roi.normalize_green(g, self.k)
        self.blue, _ = cell_roi.load_nuclei(blue_path, self.k)
        self.red, _ = mina.load_channel(red_path, 0)
        self.green, _ = mina.load_channel(green_path, 1)
        self.method, self.sens = nuclei_method, sensitivity
        self._img, self._auto = {}, {}

    def image(self, key):
        """(image the threshold is applied to, automatic threshold)."""
        if key not in self._img:
            if key == 'thr_nuclei':
                im, t = nuclei.threshold_image(self.blue, self.k, self.method)
            elif key == 'cell_fg_level':
                im, t = cell_roi.cell_density(self.red, self.px), 0.12
            elif key == 'thr_mito':
                p0 = dict(mito_objects.P, sigma=mito_objects.P['sigma'] * self.k)
                im = mito_objects.preprocess(self.red, self.px, p0)
                fg = cell_roi.cell_density(self.red, self.px) > 0.12
                # reference only: the automatic mito threshold is per cell (local mean AND 0.5 x Otsu of the cell)
                from skimage import filters
                t = float(0.5 * filters.threshold_otsu(im[fg])) if fg.any() else float('nan')
            elif key == 'thr_green_bright':
                im = green_cells.bright_image(self.gn, self.k); t = green_cells.auto_threshold(im)
            elif key == 'thr_puncta':
                g = np.clip(self.green.astype(float) - self.norm['offset'], 0, None)
                im = green_q.tophat(g, self.k)
                fg = cell_roi.cell_density(self.red, self.px) > 0.12
                t = green_q.auto_puncta_threshold(im, fg, self.sens)
            else:
                raise KeyError(key)
            self._img[key], self._auto[key] = im, float(t)
        return self._img[key], self._auto[key]

    def mask(self, key, value):
        im, t = self.image(key)
        v = value if value > 0 else t
        m = im > v
        if key == 'thr_puncta':
            m = morphology.remove_small_objects(m, max_size=max(1, round(3 * self.k * self.k)))
        return m

    def range(self, key):
        """Sensible slider range (low, high) for the threshold of `key`."""
        im, t = self.image(key)
        if key == 'cell_fg_level':
            return 0.0, 1.0
        lo, hi = np.percentile(im, [1, 99.9])
        hi = max(hi, (t if np.isfinite(t) else 0) * 2, lo + 1e-6)
        return float(max(lo, 0.0)), float(hi)

    def display(self, key):
        """Grey image to draw the mask on (0-1)."""
        base = {'thr_nuclei': self.blue, 'cell_fg_level': self.red, 'thr_mito': self.red,
                'thr_green_bright': self.green, 'thr_puncta': self.green}[key].astype(float)
        hi = np.percentile(base, 99.7) or 1.0
        return np.clip(base / hi, 0, 1)
