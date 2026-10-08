"""Manual review of the cell ROIs between detection and analysis (window of the Run tab).

Tools (as in Cytosol Viewer): select (V), add a drawn area to the selected cell (A), subtract a drawn
area (S), draw a new cell (N), delete the selected cell (Delete), merge the selected cell with the next
one clicked (M), undo / redo (Cmd+Z / Cmd+Shift+Z). Zoom and pan with the toolbar or the mouse wheel.
"""
import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.widgets import LassoSelector
from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QDialog, QHBoxLayout, QLabel, QListWidget, QListWidgetItem,
                               QMessageBox, QPushButton, QSplitter, QVBoxLayout, QWidget)
from scipy import ndimage as ndi
from skimage import segmentation

from .roi_edit import RoiEdit

TOOLS = (('select', 'Select', 'V', 'Click a cell to select it'),
         ('add', 'Add', 'A', 'Draw around an area to add it to the selected cell'),
         ('sub', 'Subtract', 'S', 'Draw around an area to remove it from the selected cell '
                                  '(from every cell it touches if none is selected)'),
         ('new', 'New cell', 'N', 'Draw the outline of a cell that was missed'))


class RoiReviewDialog(QDialog):
    """`items`: list of dict(name, review) where review = pipeline.detect_rois(...) output.
    After exec(), `results()` gives the edited ROIs (labels, nuclei, seed) per item."""

    def __init__(self, items, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Review cell ROIs')
        self.resize(1400, 900)
        self.items = items
        self.edits = [RoiEdit(it['review']['labels'], it['review']['nuclei'], it['review']['seed']) for it in items]
        self.seen = [False] * len(items)
        self.cur = 0; self.sel = 0; self.tool = 'select'; self.merge_from = 0

        # left: image list
        self.list = QListWidget(); self.list.setMinimumWidth(220)
        self.list.currentRowChanged.connect(self.show_item)
        left = QWidget(); lv = QVBoxLayout(left); lv.setContentsMargins(0, 0, 0, 0)
        lv.addWidget(QLabel('Image sets')); lv.addWidget(self.list, 1)
        self.count = QLabel(); lv.addWidget(self.count)

        # toolbar
        self.tool_group = QButtonGroup(self); self.tool_btns = {}
        tools = QHBoxLayout()
        for key, label, sc, tip in TOOLS:
            b = QPushButton(f'{label} ({sc})'); b.setCheckable(True); b.setToolTip(tip)
            self.tool_group.addButton(b); self.tool_btns[key] = b
            b.clicked.connect(lambda _=False, k=key: self.set_tool(k))
            tools.addWidget(b)
            QShortcut(QKeySequence(sc), self, activated=lambda k=key: self.set_tool(k))
        self.tool_btns['select'].setChecked(True)
        for text, sc, fn, tip in (('Delete cell', 'Delete', self.delete_cell, 'Delete the selected cell'),
                                  ('Merge…', 'M', self.start_merge, 'Merge the selected cell with the next one you click'),
                                  ('Undo', QKeySequence.Undo, self.undo, ''), ('Redo', QKeySequence.Redo, self.redo, ''),
                                  ('Reset image', None, self.reset, 'Back to the detected ROIs of this image')):
            b = QPushButton(text); b.clicked.connect(fn)
            if tip:
                b.setToolTip(tip)
            tools.addWidget(b)
            if sc is not None:
                QShortcut(QKeySequence(sc), self, activated=fn)
        QShortcut(QKeySequence(Qt.Key_Backspace), self, activated=self.delete_cell)
        QShortcut(QKeySequence(Qt.Key_Escape), self, activated=lambda: (self.set_tool('select'), self.select(0)))
        tools.addStretch(1)
        view = QHBoxLayout()
        self.ch = {}
        for key, label in (('r', 'Mito'), ('g', 'Green'), ('b', 'Nuclei')):
            cb = QCheckBox(label); cb.setChecked(True); cb.toggled.connect(self.redraw); self.ch[key] = cb
            view.addWidget(cb)
        self.show_nuc = QCheckBox('Nucleus outlines'); self.show_nuc.setChecked(True); self.show_nuc.toggled.connect(self.redraw)
        self.show_num = QCheckBox('Numbers'); self.show_num.setChecked(True); self.show_num.toggled.connect(self.redraw)
        view.addWidget(self.show_nuc); view.addWidget(self.show_num); view.addStretch(1)

        # canvas
        self.fig = Figure(figsize=(9, 9)); self.canvas = FigureCanvasQTAgg(self.fig)
        self.ax = self.fig.add_axes([0, 0, 1, 1]); self.ax.set_axis_off()
        self.nav = NavigationToolbar2QT(self.canvas, self)
        self.canvas.mpl_connect('button_press_event', self.on_click)
        self.canvas.mpl_connect('scroll_event', self.on_scroll)
        self.lasso = None
        self.hint = QLabel(); self.hint.setWordWrap(True)
        self.msg = QLabel(); self.msg.setStyleSheet('color: #2a6ad8;')

        center = QWidget(); cv = QVBoxLayout(center); cv.setContentsMargins(0, 0, 0, 0)
        cv.addLayout(tools); cv.addLayout(view); cv.addWidget(self.nav); cv.addWidget(self.canvas, 1)
        cv.addWidget(self.hint); cv.addWidget(self.msg)
        split = QSplitter(); split.addWidget(left); split.addWidget(center); split.setSizes([240, 1160])

        prev_b = QPushButton('◀ Previous image'); next_b = QPushButton('Next image ▶')
        prev_b.clicked.connect(lambda: self.list.setCurrentRow(max(0, self.cur - 1)))
        next_b.clicked.connect(lambda: self.list.setCurrentRow(min(len(self.items) - 1, self.cur + 1)))
        ok = QPushButton('Confirm ROIs and continue analysis'); ok.setDefault(True); ok.setMinimumHeight(32)
        ok.clicked.connect(self.confirm)
        cancel = QPushButton('Cancel'); cancel.clicked.connect(self.reject)
        bottom = QHBoxLayout()
        for w in (prev_b, next_b):
            bottom.addWidget(w)
        bottom.addStretch(1); bottom.addWidget(cancel); bottom.addWidget(ok)
        lay = QVBoxLayout(self); lay.addWidget(split, 1); lay.addLayout(bottom)

        for it in items:
            self.list.addItem(QListWidgetItem(it['name']))
        self.set_tool('select')
        self.list.setCurrentRow(0)

    # ---- state ----
    @property
    def ed(self):
        return self.edits[self.cur]

    def results(self):
        return [e.result() for e in self.edits]

    def refresh_list(self):
        for i, (it, e) in enumerate(zip(self.items, self.edits)):
            n = len(e.ids()); n_ok = sum(not e.touches_frame(l) for l in e.ids())
            mark = '✎ ' if e.changed else ('✓ ' if self.seen[i] else '   ')
            self.list.item(i).setText(f"{mark}{it['name']}  ({n_ok}/{n} cells)")
        self.count.setText(f'{sum(self.seen)} of {len(self.items)} images viewed, '
                           f'{sum(e.changed for e in self.edits)} edited')

    def show_item(self, i):
        if i < 0:
            return
        self.cur = i; self.sel = 0; self.merge_from = 0; self.seen[i] = True
        self.ax.clear(); self.ax.set_axis_off()
        self._img = None
        self.redraw(reset_view=True)
        self.refresh_list()

    # ---- tools ----
    def set_tool(self, tool):
        self.tool = tool; self.merge_from = 0
        self.tool_btns[tool].setChecked(True)
        if self.lasso is not None:
            self.lasso.set_active(tool != 'select')
        self.hint.setText({'select': 'Click a cell to select it. Delete removes it; Merge… joins it with the next cell '
                                     'you click. Cells with grey outlines and numbers touch the image border and are not analysed; * = drawn by hand.',
                           'add': 'Draw around the area to add to the selected cell (it takes the area from other cells).',
                           'sub': 'Draw around the area to cut from the selected cell (or from every cell it touches).',
                           'new': 'Draw the outline of a missed cell (only where no other cell is).'}[tool])

    def _make_lasso(self):
        if self.lasso is not None:
            self.lasso.disconnect_events()
        self.lasso = LassoSelector(self.ax, self.on_lasso, useblit=True,
                                   props=dict(color='#ffd400', linewidth=1.5))
        self.lasso.set_active(self.tool != 'select')

    def on_lasso(self, verts):
        if self.nav.mode or self.tool == 'select':
            return
        mask = self.ed.polygon_mask(verts)
        if not mask.any():
            return
        e = self.ed
        if self.tool == 'new':
            l = e.new_cell(mask)
            self.say(f'cell {l} drawn' if l else 'The drawn area lies inside other cells; nothing added.')
            if l:
                self.sel = l
        elif self.tool == 'add':
            l = self.sel
            if not l:
                ov = e.labels[mask]; ov = ov[ov > 0]
                l = int(np.bincount(ov).argmax()) if len(ov) else 0
            if not l:
                l = e.new_cell(mask); self.sel = l
                self.say(f'No cell selected: drawn as new cell {l}' if l else 'Nothing to add.')
            else:
                gained = e.add(l, mask); self.sel = l
                self.say(f'cell {l}: +{gained} px' if gained > 0 else
                         f'cell {l}: the area must touch the cell to be added.')
        elif self.tool == 'sub':
            removed = e.subtract(self.sel, mask)
            if self.sel and self.sel not in e.ids():
                self.say(f'cell {self.sel} removed completely'); self.sel = 0
            else:
                self.say(f'−{removed} px' if removed else
                         'Nothing removed: draw across the cell border (a cell cannot have a hole).')
        self.redraw(); self.refresh_list()

    def on_click(self, ev):
        if ev.inaxes is not self.ax or self.nav.mode or ev.button != 1 or self.tool != 'select':
            return
        l = self.ed.label_at(ev.xdata, ev.ydata)
        if self.merge_from:
            a, self.merge_from = self.merge_from, 0
            if l and l != a and self.ed.merge(a, l):
                self.say(f'cell {l} merged into cell {a}'); self.sel = a
                self.redraw(); self.refresh_list()
            elif l and l != a:
                self.say(f'cells {a} and {l} do not touch, so they cannot be one cell. Use Add (A) to draw the '
                         'part between them first.')
            else:
                self.say('Merge cancelled (click a different cell).')
            return
        self.select(l)

    def on_scroll(self, ev):
        if ev.inaxes is not self.ax:
            return
        f = 0.8 if ev.button == 'up' else 1.25
        x0, x1 = self.ax.get_xlim(); y0, y1 = self.ax.get_ylim()
        cx, cy = ev.xdata, ev.ydata
        self.ax.set_xlim(cx - (cx - x0) * f, cx + (x1 - cx) * f)
        self.ax.set_ylim(cy - (cy - y0) * f, cy + (y1 - cy) * f)
        self.canvas.draw_idle()

    def select(self, l):
        self.sel = l
        if l:
            area = int((self.ed.labels == l).sum()); px = self.items[self.cur]['review'].get('px') or 0
            kind = self.ed.seed.get(l, 'manual')
            self.say(f'cell {l}: {area * px * px:.0f} µm², seed: {kind}'
                     + (', touches the border (not analysed)' if self.ed.touches_frame(l) else ''))
        self.redraw()

    def start_merge(self):
        if not self.sel:
            self.say('Select the first cell, then press Merge… and click the second cell.'); return
        self.set_tool('select'); self.merge_from = self.sel
        self.say(f'Click the cell to merge into cell {self.sel} (Esc to cancel).')

    def delete_cell(self):
        if self.sel and self.ed.delete(self.sel):
            self.say(f'cell {self.sel} deleted'); self.sel = 0
            self.redraw(); self.refresh_list()

    def undo(self):
        if self.ed.undo():
            self.sel = 0; self.redraw(); self.refresh_list(); self.say('undone')

    def redo(self):
        if self.ed.redo():
            self.redraw(); self.refresh_list(); self.say('redone')

    def reset(self):
        self.ed.reset(); self.sel = 0; self.redraw(); self.refresh_list(); self.say('detected ROIs restored')

    def say(self, text):
        self.msg.setText(text)

    # ---- drawing ----
    def redraw(self, *_, reset_view=False):
        e = self.ed
        img = self.items[self.cur]['review']['image'].astype(float) / 255
        for i, k in enumerate('rgb'):
            if not self.ch[k].isChecked():
                img[..., i] = 0
        lab = e.labels
        edge = {l for l in e.ids() if e.touches_frame(l)}
        bnd = segmentation.find_boundaries(lab, mode='thick')
        if lab.shape[0] > 700:   # keep outlines visible when a large image is shown scaled down
            bnd = ndi.binary_dilation(bnd)
        edge_m = np.isin(lab, list(edge)) if edge else np.zeros(lab.shape, bool)
        img[bnd & ~edge_m] = (1.0, 1.0, 1.0)
        img[bnd & edge_m] = (0.55, 0.55, 0.55)
        if self.show_nuc.isChecked():
            img[segmentation.find_boundaries(e.nuclei, mode='thick')] = (0.3, 0.8, 1.0)
        if self.sel:
            sb = ndi.binary_dilation(segmentation.find_boundaries(lab == self.sel, mode='thick'), iterations=2)
            img[sb] = (1.0, 0.85, 0.0)
        if reset_view or not self.ax.images:
            self.ax.clear(); self.ax.set_axis_off()
            self.ax.imshow(img, interpolation='nearest')
            self.ax.set_xlim(-0.5, lab.shape[1] - 0.5); self.ax.set_ylim(lab.shape[0] - 0.5, -0.5)
            self._make_lasso()
        else:
            self.ax.images[0].set_data(img)
        for t in list(self.ax.texts):
            t.remove()
        if self.show_num.isChecked():
            for l in e.ids():
                cx, cy = e.centroid(l)
                manual = e.seed.get(l) == 'manual'
                self.ax.text(cx, cy, f'{l}' + ('*' if manual else ''), color='#ffd400' if l == self.sel else
                             ('#bbbbbb' if l in edge else 'white'), fontsize=9, ha='center', va='center',
                             weight='bold', bbox=dict(facecolor='black', alpha=0.45, pad=1, edgecolor='none'))
        self.canvas.draw_idle()

    def confirm(self):
        unseen = [self.items[i]['name'] for i, s in enumerate(self.seen) if not s]
        if unseen and QMessageBox.question(
                self, 'Review cell ROIs', f'{len(unseen)} image set(s) were not opened '
                f"({', '.join(unseen[:4])}{' …' if len(unseen) > 4 else ''}).\n"
                'Use their detected ROIs as they are and continue?') != QMessageBox.Yes:
            return
        for e in self.edits:
            if not e.ids():
                if QMessageBox.question(self, 'Review cell ROIs', 'An image set has no cells left. Continue?') \
                        != QMessageBox.Yes:
                    return
                break
        self.accept()
