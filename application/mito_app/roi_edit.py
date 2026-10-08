"""Editable cell ROIs for the manual review step (model only; the window is gui.RoiReviewDialog).

A cell is one filled piece of a label image; its nucleus (if any) carries the same label in `nuclei`.
Edits: add a drawn area to a cell, subtract it, draw a new cell, delete, merge two cells, with undo / redo.
"""
import numpy as np
from scipy import ndimage as ndi
from skimage import draw, measure


def _largest(m):
    cc = measure.label(m)
    if cc.max() <= 1:
        return m
    return cc == np.argmax(np.bincount(cc.ravel())[1:]) + 1


class RoiEdit:
    def __init__(self, labels, nuclei, seed, max_undo=40):
        self.labels = np.asarray(labels, np.int32).copy()
        self.nuclei = np.asarray(nuclei, np.int32).copy()
        self.seed = dict(seed)
        self.orig = (self.labels.copy(), self.nuclei.copy(), dict(self.seed))
        self._undo, self._redo, self.max_undo = [], [], max_undo
        self.changed = False

    # ---- state ----
    @property
    def shape(self):
        return self.labels.shape

    def ids(self):
        return [int(i) for i in np.unique(self.labels) if i > 0]

    def result(self):
        return dict(labels=self.labels, nuclei=self.nuclei, seed=dict(self.seed))

    def _push(self):
        self._undo.append((self.labels.copy(), self.nuclei.copy(), dict(self.seed)))
        del self._undo[:-self.max_undo]
        self._redo.clear()
        self.changed = True

    def undo(self):
        if not self._undo:
            return False
        self._redo.append((self.labels.copy(), self.nuclei.copy(), dict(self.seed)))
        self.labels, self.nuclei, self.seed = self._undo.pop()
        return True

    def redo(self):
        if not self._redo:
            return False
        self._undo.append((self.labels.copy(), self.nuclei.copy(), dict(self.seed)))
        self.labels, self.nuclei, self.seed = self._redo.pop()
        return True

    def reset(self):
        self._push()
        self.labels, self.nuclei, self.seed = (self.orig[0].copy(), self.orig[1].copy(), dict(self.orig[2]))

    # ---- geometry ----
    def polygon_mask(self, xy):
        """Mask of a polygon given as (N, 2) x, y image coordinates (pixel centres inside)."""
        xy = np.asarray(xy, float)
        if len(xy) < 3:
            return np.zeros(self.shape, bool)
        rr, cc = draw.polygon(xy[:, 1] - 0.5, xy[:, 0] - 0.5, self.shape)
        m = np.zeros(self.shape, bool)
        m[rr, cc] = True
        return m

    def label_at(self, x, y):
        H, W = self.shape
        xi, yi = int(x), int(y)
        return int(self.labels[yi, xi]) if 0 <= xi < W and 0 <= yi < H else 0

    def _erase(self, l):
        self.labels[self.labels == l] = 0
        self.nuclei[self.nuclei == l] = 0
        self.seed.pop(l, None)

    def _tidy(self, l, prefer=None):
        """Keep a cell as one filled piece: the piece overlapping `prefer`, else its nucleus, else the largest."""
        m = self.labels == l
        if not m.any():
            self._erase(l)
            return
        cc = measure.label(m)
        keep = _largest(m)
        if cc.max() > 1:
            for ref in (prefer, self.nuclei == l):
                if ref is not None and (ref & m).any():
                    votes = np.bincount(cc[ref & m], minlength=cc.max() + 1); votes[0] = 0
                    keep = cc == np.argmax(votes)
                    break
        self.labels[m & ~keep] = 0
        self.nuclei[(self.nuclei == l) & ~keep] = 0
        self.labels[ndi.binary_fill_holes(keep) & (self.labels == 0)] = l

    # ---- edits ----
    def add(self, l, mask):
        """Add the drawn area to cell l (taking it from other cells). Returns pixels gained."""
        before = int((self.labels == l).sum())
        prev = self.labels == l
        self._push()
        others = set(np.unique(self.labels[mask]).tolist()) - {0, l}
        self.labels[mask] = l
        for o in others:
            self.nuclei[(self.nuclei == o) & mask] = 0
            self._tidy(o)
        self._tidy(l, prefer=prev)
        return int((self.labels == l).sum()) - before

    def subtract(self, l, mask):
        """Remove the drawn area from cell l (from every cell it touches if l is 0). Returns pixels removed."""
        targets = [l] if l else [int(i) for i in np.unique(self.labels[mask]) if i > 0]
        if not targets:
            return 0
        before = int((self.labels > 0).sum())
        self._push()
        for t in targets:
            self.labels[mask & (self.labels == t)] = 0
            self.nuclei[mask & (self.nuclei == t)] = 0
            self._tidy(t)
        removed = before - int((self.labels > 0).sum())
        if not removed:
            self.undo(); self._redo.clear()
        return removed

    def new_cell(self, mask):
        """The drawn area becomes a new cell where no other cell is. Returns its label (0 if nothing left)."""
        m = mask & (self.labels == 0)
        if not m.any():
            return 0
        self._push()
        l = int(self.labels.max()) + 1
        self.labels[_largest(m)] = l
        self.seed[l] = 'manual'
        return l

    def delete(self, l):
        if l and (self.labels == l).any():
            self._push()
            self._erase(l)
            return True
        return False

    def touching(self, a, b):
        ma = self.labels == a
        return bool((ndi.binary_dilation(ma, iterations=2) & (self.labels == b)).any())

    def merge(self, a, b):
        """Merge cell b into cell a (b's nucleus becomes part of a). Only touching cells can be merged."""
        if not a or not b or a == b or not self.touching(a, b):
            return False
        self._push()
        self.labels[self.labels == b] = a
        self.nuclei[self.nuclei == b] = a
        self.seed.pop(b, None)
        self._tidy(a)
        return True

    # ---- display ----
    def touches_frame(self, l, margin=3):
        m = self.labels == l
        e = margin + 1
        return bool(m[:e].any() or m[-e:].any() or m[:, :e].any() or m[:, -e:].any())

    def centroid(self, l):
        ys, xs = np.nonzero(self.labels == l)
        return (float(xs.mean()), float(ys.mean())) if len(xs) else (0.0, 0.0)
