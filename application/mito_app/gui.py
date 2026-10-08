"""PySide6 desktop GUI for the per-cell MiNA + green-on-mito pipeline."""
import math, os, sys, traceback

os.environ.setdefault('QT_API', 'pyside6')

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from scipy.stats import t as scipy_t
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QObject, QSortFilterProxyModel, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QAction, QColor, QDesktopServices, QIcon, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QDialog, QDialogButtonBox, QHBoxLayout, QHeaderView, QInputDialog, QLabel, QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea,
    QSplitter, QTableView, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__, imgio, nuclei, pipeline, stats

TIFF_FILTER = 'Images (*.tif *.tiff *.czi);;TIFF images (*.tif *.tiff);;Zeiss CZI (*.czi);;All files (*)'
GREEN_COLORS = {'positive': '#2CA02C', 'negative': '#7F7F7F'}
CELL_COLORS = ['#4C9BE8', '#F28E2B', '#59A14F', '#E15759', '#B07AA1', '#EDC948', '#76B7B2', '#FF9DA7',
               '#9C755F', '#BAB0AC']


def tight(fig):
    """tight_layout that never raises (a hidden canvas has no size yet)."""
    try:
        fig.tight_layout()
    except (ValueError, np.linalg.LinAlgError):
        pass


class AspectBox(QWidget):
    """Holds one widget (a plot canvas) at a fixed width:height ratio, as large as fits and centred,
    so a plot keeps its proportions however wide or tall the window is."""

    def __init__(self, child, ratio=4 / 3):
        super().__init__()
        self.child, self.ratio = child, ratio
        child.setParent(self)
        self.setMinimumSize(320, 240)

    def resizeEvent(self, e):
        super().resizeEvent(e)
        w, h = self.width(), self.height()
        cw, ch = (w, int(w / self.ratio)) if w / max(h, 1) < self.ratio else (int(h * self.ratio), h)
        self.child.setGeometry((w - cw) // 2, (h - ch) // 2, cw, ch)


def short(text, n=18):
    """Shorten a label in the middle ('Condition_01_…_ctrl') so tick labels do not run into each other."""
    text = str(text)
    return text if len(text) <= n else text[:n // 2 + 2] + '…' + text[-(n - n // 2 - 3):]


def short_groups(groups, n=18):
    """Short labels for a set of group names: a prefix shared by all names (cut at _ - space) is dropped
    ('Condition_01_ctrl', 'Condition_02_ctrl' -> '01_ctrl', '02_ctrl'), then each is shortened in the middle."""
    groups = [str(g) for g in groups]
    pre = os.path.commonprefix(groups) if len(groups) > 1 else ''
    cut = max(pre.rfind(c) for c in '_- ') + 1
    if cut and all(len(g) > cut for g in groups) and any(len(g) > n for g in groups):
        return {g: short(g[cut:], n) for g in groups}
    return {g: short(g, n) for g in groups}


def group_list(groups, n=6):
    groups = list(groups)
    return ', '.join(groups[:n]) + (f', … (+{len(groups) - n})' if len(groups) > n else '')


def safe_draw(fn):
    """Matplotlib cannot lay out a canvas that has no size yet (a hidden tab right after loading results):
    skip the drawing then and redo it when the tab is shown (see redraw_when_shown)."""
    def wrapper(self, *a, **k):
        try:
            return fn(self, *a, **k)
        except (ValueError, np.linalg.LinAlgError):
            if self.isVisible():
                raise
            self._stale = True
    wrapper.__name__ = fn.__name__
    return wrapper


class RedrawWhenShown:
    """Mixin for tabs with safe_draw methods: redraw everything once the tab becomes visible."""
    _stale = False

    def showEvent(self, e):
        super().showEvent(e)
        if self._stale:
            self._stale = False
            self.redraw()


def resource(name):
    base = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(base, 'resources', name)


class Worker(QObject):
    """Runs the image sets one after another; a failed set is reported and the batch goes on."""
    log = Signal(str)
    started = Signal(int)
    job_done = Signal(int, object)
    job_failed = Signal(int, str)
    finished = Signal(int, int)

    def __init__(self, jobs, params, base_out='', mode='analyse'):
        super().__init__()
        self.jobs, self.params, self.base_out, self.mode = jobs, params, base_out, mode
        self.stop = False

    def run(self):
        results, failed = [], 0
        detect = self.mode == 'detect'
        for n, (row, job) in enumerate(self.jobs, 1):
            if self.stop:
                break
            self.started.emit(row)
            self.log.emit(f'=== [{n}/{len(self.jobs)}] {job["dataset"]} / {job["name"]}'
                          + (' (cell ROIs for review) ===' if detect else ' ==='))
            try:
                ov = dict(job.get('overrides') or {})
                if job.get('px') and not job.get('px_tiff') and not self.params.pixel_size_um:  # entered per image
                    ov['pixel_size_um'] = job['px']
                if detect:
                    res = pipeline.detect_rois(job['red'], job['green'], job['blue'], self.params, log=self.log.emit,
                                               overrides=ov)
                    if job.get('roi_start'):   # start from the ROIs reviewed in an earlier run
                        start = job['roi_start']
                        if np.shape(start['labels']) == res['labels'].shape:
                            res.update(labels=start['labels'], nuclei=start['nuclei'], seed=start['seed'])
                            self.log.emit('     starting from the ROIs reviewed in an earlier run')
                else:
                    res = pipeline.run(job['red'], job['green'], job['blue'], job['outdir'], self.params, job['name'],
                                       log=self.log.emit, overrides=ov, roi_edit=job.get('roi_edit'))
                results.append(res)
                self.job_done.emit(row, res)
            except Exception as e:  # report to the user instead of crashing the app
                failed += 1
                self.log.emit(traceback.format_exc())
                self.job_failed.emit(row, f'{type(e).__name__}: {e}')
        try:
            if detect:
                self.finished.emit(len(results), failed)
                return
            pipeline.write_group_tables(results)
            if results and self.base_out:
                pipeline.write_group_exports(self.base_out, log=self.log.emit)
        except Exception:
            self.log.emit(traceback.format_exc())
        self.finished.emit(len(results), failed)


class JobTable(QTableWidget):
    """Input table, one image set per row (three TIFFs, or one CZI file with a channel per role).
    Accepts dropped TIFFs / CZIs / folders."""
    COLS = ('Group', 'Dataset', 'Preset', 'Sample', 'µm/px', 'Red', 'Green', 'Blue', 'Thresholds', 'Status')
    GROUP, DATASET, PRESET, NAME, PX, RED, GREEN, BLUE, THR, STATUS = range(10)
    CH_COL = {'red': 5, 'green': 6, 'blue': 7}
    dropped = Signal(list)

    def __init__(self):
        super().__init__(0, len(self.COLS))
        self.setHorizontalHeaderLabels(self.COLS)
        self.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed
                             | QAbstractItemView.SelectedClicked)
        self.setAlternatingRowColors(True)
        self.verticalHeader().setDefaultSectionSize(22)
        h = self.horizontalHeader(); h.setSectionResizeMode(QHeaderView.Interactive); h.setStretchLastSection(True)
        for c, w in zip(range(10), (70, 70, 50, 80, 58, 84, 84, 84, 70, 60)):
            self.setColumnWidth(c, w)
        self.setAcceptDrops(True)

    def dragEnterEvent(self, e):
        if e.mimeData().hasUrls():
            e.acceptProposedAction()
        else:
            super().dragEnterEvent(e)

    dragMoveEvent = dragEnterEvent

    def dropEvent(self, e):
        if e.mimeData().hasUrls():
            self.dropped.emit([u.toLocalFile() for u in e.mimeData().urls() if u.isLocalFile()])
            e.acceptProposedAction()
        else:
            super().dropEvent(e)

    def add_set(self, s, preset):
        r = self.rowCount(); self.insertRow(r)
        dataset = s.get('dataset') or os.path.basename(os.path.dirname(s['red']))
        texts = (dataset, dataset, preset, s.get('name') or pipeline.sample_name(s['red']))
        for c, t in enumerate(texts):
            self.setItem(r, c, QTableWidgetItem(t))
        self.item(r, self.DATASET).setData(Qt.UserRole, s.get('root') or os.path.dirname(s['red']))
        self.item(r, self.GROUP).setToolTip('Group for the group comparison in the Analysis tab (default: the dataset). '
                                            'Select rows and press "Set group…" to group several sets.')
        for k, c in self.CH_COL.items():
            self.set_path(r, k, s[k])
        self.item(r, self.NAME).setData(Qt.UserRole, s.get('czi'))
        px, src = pipeline.cell_roi.read_pixel_size(s['green'])
        self.setItem(r, self.PX, QTableWidgetItem(''))
        self.set_px(r, px, src)
        th = QTableWidgetItem('auto'); th.setFlags(th.flags() & ~Qt.ItemIsEditable); th.setData(Qt.UserRole, {})
        th.setToolTip('Per-image manual thresholds (set them in Thresholds → Preview / adjust…, "This image only")')
        self.setItem(r, self.THR, th)
        st = QTableWidgetItem('ready'); st.setFlags(st.flags() & ~Qt.ItemIsEditable)
        self.setItem(r, self.STATUS, st)

    def set_px(self, r, px, src=''):
        """Pixel size cell: the value read from the TIFF (or entered); '?' in red when the TIFF has none."""
        it = self.item(r, self.PX)
        it.setData(Qt.UserRole, (px or 0.0, src))
        if px:
            it.setText(f'{px:.4f}'); it.setBackground(QColor(0, 0, 0, 0))
            it.setToolTip(f'{px:.5f} µm per pixel ({src}). Double-click to change.')
        else:
            it.setText('?'); it.setBackground(QColor(255, 120, 120, 90))
            it.setToolTip(f'No pixel size in the TIFF ({src}). Double-click and type the µm per pixel; '
                          'you are asked for it when you press Run all.')

    def px(self, r):
        """Pixel size for row r (0 = missing) and whether it is the TIFF's own value."""
        stored, src = self.item(r, self.PX).data(Qt.UserRole) or (0.0, '')
        try:
            v = float(self.text(r, self.PX).replace(',', '.'))
        except ValueError:
            return 0.0, False
        if stored and src != 'entered' and abs(v - round(stored, 4)) < 1e-9:
            return stored, True          # unchanged TIFF value, full precision
        return (v, False) if v > 0 else (0.0, False)

    def overrides(self, r):
        it = self.item(r, self.THR)
        return dict(it.data(Qt.UserRole) or {}) if it else {}

    def set_overrides(self, r, d):
        it = self.item(r, self.THR)
        d = {k: v for k, v in d.items() if v}
        it.setData(Qt.UserRole, d)
        labels = {k: lab for k, lab, *_ in pipeline.THRESHOLDS}
        it.setText(', '.join(f'{labels[k].split()[0]} {v:g}' for k, v in d.items()) or 'auto')
        it.setToolTip('This image uses: ' + (', '.join(f'{labels[k]} = {v:g}' for k, v in d.items()) or
                                             'the Options thresholds') + '. Change it in Thresholds → Preview / adjust…')

    def set_path(self, r, key, path):
        it = QTableWidgetItem(imgio.display_name(path)); it.setFlags(it.flags() & ~Qt.ItemIsEditable)
        it.setData(Qt.UserRole, path)
        czi, ch = imgio.split_ref(path)
        it.setToolTip(f'{czi}, channel {ch}. Double-click (or Channel roles…) to change the channel roles.'
                      if ch is not None else path)
        self.setItem(r, self.CH_COL[key], it)

    def czi(self, r):
        """CZI file of row r, or None for a TIFF set."""
        it = self.item(r, self.NAME)
        return it.data(Qt.UserRole) if it else None

    def roles(self, r):
        """{'red': channel, 'green': channel, 'blue': channel} of a CZI row."""
        return {k: imgio.split_ref(self.item(r, c).data(Qt.UserRole))[1] for k, c in self.CH_COL.items()}

    def set_roles(self, r, roles):
        for k in self.CH_COL:
            self.set_path(r, k, imgio.channel_ref(self.czi(r), roles[k]))

    def text(self, r, c):
        it = self.item(r, c)
        return it.text().strip() if it else ''

    def job(self, r):
        d = {k: self.item(r, c).data(Qt.UserRole) for k, c in self.CH_COL.items()}
        d.update(dataset=self.text(r, self.DATASET) or 'dataset', preset=self.text(r, self.PRESET),
                 name=self.text(r, self.NAME) or pipeline.sample_name(d['red']),
                 root=self.item(r, self.DATASET).data(Qt.UserRole) or os.path.dirname(d['red']))
        d['group'] = self.text(r, self.GROUP) or d['dataset']
        d['px'], d['px_tiff'] = self.px(r)
        d['overrides'] = self.overrides(r)
        d['czi'] = self.czi(r)
        return d

    def groups(self):
        return sorted({self.text(r, self.GROUP) for r in range(self.rowCount())} - {''})

    def set_status(self, r, text, result=None, tip=''):
        it = self.item(r, self.STATUS); it.setText(text); it.setToolTip(tip)
        if result is not None:
            it.setData(Qt.UserRole, result)

    def result(self, r):
        it = self.item(r, self.STATUS) if 0 <= r < self.rowCount() else None
        return it.data(Qt.UserRole) if it else None


class DictTableModel(QAbstractTableModel):
    """Read-only table over a list of dicts; sorts numerically via UserRole."""

    def __init__(self, rows=None):
        super().__init__()
        self.set_rows(rows or [])

    def set_rows(self, rows):
        self.beginResetModel()
        self.rows = rows
        self.cols = list(rows[0]) if rows else []
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.cols)

    def data(self, idx, role=Qt.DisplayRole):
        v = self.rows[idx.row()].get(self.cols[idx.column()])
        if role == Qt.DisplayRole:
            if isinstance(v, float):
                return '' if math.isnan(v) else f'{v:.4g}'
            return str(v)
        if role == Qt.UserRole:
            if isinstance(v, (bool, np.bool_)):
                return int(v)
            return v if isinstance(v, (int, float)) and not (isinstance(v, float) and math.isnan(v)) else (
                -math.inf if isinstance(v, float) else str(v))
        if role == Qt.TextAlignmentRole and isinstance(v, (int, float)) and not isinstance(v, bool):
            return int(Qt.AlignRight | Qt.AlignVCenter)

    def headerData(self, s, orient, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self.cols[s] if orient == Qt.Horizontal else str(s + 1)


class TableTab(QWidget):
    """Sortable table with an optional per-cell filter."""

    def __init__(self, cell_filter=True):
        super().__init__()
        self.model = DictTableModel()
        self.proxy = QSortFilterProxyModel(); self.proxy.setSourceModel(self.model)
        self.proxy.setSortRole(Qt.UserRole)
        self.view = QTableView(); self.view.setModel(self.proxy); self.view.setSortingEnabled(True)
        self.view.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.view.setAlternatingRowColors(True)
        lay = QVBoxLayout(self)
        self.cell = QComboBox()
        self.group = green_group_combo()
        self.count = QLabel()
        self.all_rows = []
        if cell_filter:
            row = QHBoxLayout()
            for w in (QLabel('Cells'), self.group, QLabel('Cell'), self.cell):
                row.addWidget(w)
            row.addStretch(1); row.addWidget(self.count); lay.addLayout(row)
            self.cell.currentTextChanged.connect(self.apply_filter)
            self.group.currentIndexChanged.connect(self.group_changed)
        lay.addWidget(self.view)

    def set_rows(self, rows):
        self.all_rows = rows
        self.group_changed()

    def group_rows(self):
        return filter_group(self.all_rows, self.group.currentData())

    def group_changed(self, *_):
        self.cell.blockSignals(True); self.cell.clear()
        self.cell.addItems(['All'] + sorted({r.get('cell', '') for r in self.group_rows()}))
        self.cell.blockSignals(False)
        self.apply_filter('All')

    def apply_filter(self, text):
        rows = self.group_rows()
        if text not in ('All', ''):
            rows = [r for r in rows if r.get('cell') == text]
        self.model.set_rows(rows)
        self.view.sortByColumn(-1, Qt.AscendingOrder); self.proxy.sort(-1)  # original order until a header is clicked
        self.view.resizeColumnsToContents()
        self.count.setText(f'{len(rows)} rows')


GREEN_GROUPS = (('All cells', ''), ('Green+ only', 'positive'), ('Green− only', 'negative'))


def green_group_combo():
    c = QComboBox()
    for label, key in GREEN_GROUPS:
        c.addItem(label, key)
    c.setToolTip('Green-positive and green-negative cells are measured and exported separately '
                 '(…_green_pos / …_green_neg files)')
    return c


def filter_group(rows, key):
    return rows if not key else [r for r in rows if r.get('green_status') == key]


class ImageView(QScrollArea):
    """Scrollable image that fits the width of the panel."""

    def __init__(self):
        super().__init__()
        self.label = QLabel('Run an analysis to see the result here.')
        self.label.setAlignment(Qt.AlignCenter)
        self.setWidget(self.label)
        self.setWidgetResizable(True)
        self.pix = None

    def set_image(self, path):
        self.pix = QPixmap(path) if path and os.path.exists(path) else None
        self.label.setText('' if self.pix else 'No image')
        self._fit()

    def resizeEvent(self, e):
        super().resizeEvent(e)
        self._fit()

    def _fit(self):
        if self.pix and not self.pix.isNull():
            w = max(self.viewport().width() - 4, 100)
            self.label.setPixmap(self.pix.scaledToWidth(w, Qt.SmoothTransformation))


class CorrelationTab(RedrawWhenShown, QWidget):
    """Heatmap of all mito x green metric pairs (click a square) + regression scatter of the chosen pair.

    The scatter is split into quadrants at the mean (or median) of X and Y, so positive (I/III) and
    negative (II/IV) relations are visible at a glance. With `pooled=True` (Analysis tab) the rows come
    from many samples and can be filtered by dataset / preset / sample and coloured by any of them."""

    def __init__(self, pooled=False):
        super().__init__()
        self.pooled = pooled
        self.data = {'cell': [], 'mito': []}
        self.level = QComboBox(); self.level.addItem('Per cell', 'cell'); self.level.addItem('Per mito object', 'mito')
        self.group = green_group_combo()
        self.x = QComboBox(); self.y = QComboBox(); self.cell = QComboBox()
        self.grp = QComboBox(); self.dataset = QComboBox(); self.preset = QComboBox(); self.sample = QComboBox()
        self.hstat = QComboBox(); self.hstat.addItem('Pearson r', 'r'); self.hstat.addItem('Spearman ρ', 'rho')
        self.center = QComboBox(); self.center.addItem('mean', 'mean'); self.center.addItem('median', 'median')
        self.color = QComboBox()
        for label, key in (('green status', 'green_status'), ('group', 'group'), ('sample', 'sample'), ('dataset', 'dataset'),
                           ('preset', 'preset'), ('none', '')):
            self.color.addItem(label, key)
        self.logx = QCheckBox('log X'); self.logy = QCheckBox('log Y')
        row1, row2 = QHBoxLayout(), QHBoxLayout()
        for w in (QLabel('Level'), self.level, QLabel('Cells'), self.group):
            row1.addWidget(w)
        if pooled:
            for w in (QLabel('Group'), self.grp, QLabel('Dataset'), self.dataset, QLabel('Preset'), self.preset,
                      QLabel('Sample'), self.sample):
                row1.addWidget(w)
        else:
            row1.addWidget(QLabel('Cell')); row1.addWidget(self.cell)
        row1.addStretch(1)
        for w in (QLabel('X (mito)'), self.x, QLabel('Y (green)'), self.y, QLabel('Heatmap'), self.hstat,
                  QLabel('Quadrants at'), self.center):
            row2.addWidget(w)
        if pooled:
            row2.addWidget(QLabel('Colour by')); row2.addWidget(self.color)
        row2.addWidget(self.logx); row2.addWidget(self.logy); row2.addStretch(1)
        self.export_btn = QPushButton('Export tables…')
        if pooled:
            row2.addWidget(self.export_btn)
        self.heat_fig = Figure(figsize=(6, 4)); self.heat = FigureCanvasQTAgg(self.heat_fig)
        self.sc_fig = Figure(figsize=(6, 4)); self.scatter = FigureCanvasQTAgg(self.sc_fig)
        self.heat.mpl_connect('button_press_event', self.on_heat_click)
        self.heat.mpl_connect('motion_notify_event', self.on_heat_hover)
        self.vals = None; self.reasons = {}
        self.stat = QLabel(); self.stat.setTextInteractionFlags(Qt.TextSelectableByMouse); self.stat.setWordWrap(True)
        split = QSplitter(Qt.Horizontal); split.addWidget(self.heat); split.addWidget(self.scatter)
        split.setSizes([500, 600])
        lay = QVBoxLayout(self); lay.addLayout(row1); lay.addLayout(row2); lay.addWidget(split, 1); lay.addWidget(self.stat)
        note = QLabel('Line: least-squares regression (Pearson r, R², p of the slope); with log X / log Y it is fitted '
                      'on log10 values. Spearman ρ is the rank version. Mito objects in the same cell are not '
                      'independent, so per-object p values are optimistic.')
        note.setWordWrap(True); lay.addWidget(note)
        self.level.currentIndexChanged.connect(self.level_changed)
        for w in (self.group, self.grp, self.dataset, self.preset, self.sample):
            w.currentIndexChanged.connect(self.filters_changed)
        for w in (self.x, self.y, self.cell, self.center, self.color):
            w.currentIndexChanged.connect(self.draw_scatter)
        for w in (self.logx, self.logy):
            w.toggled.connect(self.draw_heat)
            w.toggled.connect(self.draw_scatter)
        for w in (self.x, self.y):
            w.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.cell.currentIndexChanged.connect(self.draw_heat)
        self.hstat.currentIndexChanged.connect(self.draw_heat)
        self.export_btn.clicked.connect(self.export)

    def set_data(self, cell_rows, mito_rows):
        self.data = {'cell': cell_rows, 'mito': mito_rows}
        for combo, key in ((self.grp, 'group'), (self.dataset, 'dataset'), (self.preset, 'preset'),
                           (self.sample, 'sample')):
            combo.blockSignals(True); combo.clear()
            combo.addItems(['All'] + sorted({str(r.get(key, '')) for r in cell_rows + mito_rows}))
            combo.blockSignals(False)
        self.level_changed()

    def metrics(self):
        lv = self.level.currentData()
        xs, ys = (stats.CELL_X, stats.CELL_Y) if lv == 'cell' else (stats.MITO_X, stats.MITO_Y)
        rows = self.data[lv]
        return stats.available(xs, rows), stats.available(ys, rows)

    def base_rows(self, level=None):
        rows = filter_group(self.data[level or self.level.currentData()], self.group.currentData())
        if self.pooled:
            for combo, key in ((self.grp, 'group'), (self.dataset, 'dataset'), (self.preset, 'preset'),
                           (self.sample, 'sample')):
                v = combo.currentText()
                if v not in ('All', ''):
                    rows = [r for r in rows if str(r.get(key, '')) == v]
        return rows

    def rows(self):
        rows = self.base_rows()
        c = self.cell.currentText()
        if self.pooled or self.level.currentData() == 'cell' or c in ('All', ''):
            return rows
        return [r for r in rows if r['cell'] == c]

    def level_changed(self, *_):
        xs, ys = self.metrics()
        for combo, items in ((self.x, xs), (self.y, ys)):
            combo.blockSignals(True); combo.clear()
            for k, label in items:
                combo.addItem(label, k)
            combo.blockSignals(False)
        mito = self.level.currentData() == 'mito'
        # object sizes span orders of magnitude; puncta counts include 0, which log would drop
        for cb, on in ((self.logx, mito), (self.logy, False)):
            cb.blockSignals(True); cb.setChecked(on); cb.blockSignals(False)
        self.x.setCurrentIndex(0); self.y.setCurrentIndex(3 if mito else 1)  # length vs puncta / footprint vs green on mito
        self.filters_changed()

    def redraw(self):
        self.draw_heat(); self.draw_scatter()

    def filters_changed(self, *_):
        self.cell.blockSignals(True); self.cell.clear()
        self.cell.addItems(['All'] + sorted({r['cell'] for r in self.base_rows('mito')}))
        self.cell.setEnabled(self.level.currentData() == 'mito')
        self.cell.blockSignals(False)
        self.draw_heat(); self.draw_scatter()

    @safe_draw
    def draw_heat(self, *_):
        xs, ys = self.metrics(); rows = self.rows()
        use_r = self.hstat.currentData() == 'r'
        vals = np.full((len(ys), len(xs)), np.nan)
        for j, (xk, _) in enumerate(xs):
            for i, (yk, _) in enumerate(ys):
                x, y = self.values(rows, xk, self.logx), self.values(rows, yk, self.logy)
                vals[i, j] = stats.regression(x, y)['r'] if use_r else stats.spearman(x, y)[0]
        name = 'Pearson r' if use_r else 'Spearman ρ'
        if use_r and (self.logx.isChecked() or self.logy.isChecked()):
            name += ' (' + ', '.join(a for a, cb in (('log X', self.logx), ('log Y', self.logy)) if cb.isChecked()) + ')'
        f = self.heat_fig; f.clear(); ax = f.add_subplot(111)
        im = ax.imshow(vals, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        ax.set_xticks(range(len(xs)), [l for _, l in xs], rotation=60, ha='right', fontsize=7)
        ax.set_yticks(range(len(ys)), [l for _, l in ys], fontsize=7)
        self.vals = vals
        self.reasons = {}
        for j, (xk, xl) in enumerate(xs):
            for i, (yk, yl) in enumerate(ys):
                if not np.isfinite(vals[i, j]):
                    self.reasons[i, j] = stats.why_nan(self.values(rows, xk, self.logx), self.values(rows, yk, self.logy),
                                                       xl, yl)
        for i in range(vals.shape[0]):
            for j in range(vals.shape[1]):
                if np.isfinite(vals[i, j]) and vals.size <= 100:  # numbers only where they fit; hover shows the rest
                    ax.text(j, i, f'{vals[i, j]:.2f}', ha='center', va='center', fontsize=6,
                            color='white' if abs(vals[i, j]) > 0.6 else 'black')
                elif not np.isfinite(vals[i, j]) and vals.size <= 400:
                    ax.text(j, i, '–', ha='center', va='center', fontsize=6, color='#999999')
        f.colorbar(im, ax=ax, fraction=0.04, label=name)
        unit = 'cells' if self.level.currentData() == 'cell' else 'mito objects'
        if len(rows) < 3:
            note = f'\nnot computable: at least 3 {unit} are needed (filter fewer groups or add images)'
        else:
            note = (' (too few for a reliable value)' if len(rows) < 8 else '') + '\nclick a square to plot it; – = not computable'
        ax.set_title(f'{name}, n = {len(rows)} {unit}' + note, fontsize=9)
        tight(f); self.heat.draw_idle()

    def on_heat_hover(self, ev):
        if ev.xdata is None or ev.ydata is None or self.vals is None:
            return
        j, i = int(round(ev.xdata)), int(round(ev.ydata))
        if 0 <= i < self.vals.shape[0] and 0 <= j < self.vals.shape[1]:
            v = self.vals[i, j]
            self.heat.setToolTip(f'{self.y.itemText(i)} vs {self.x.itemText(j)}: ' + (
                f'{self.hstat.currentText()} = {v:.3f}' if np.isfinite(v) else
                f'not computable ({self.reasons.get((i, j), "")})'))

    def on_heat_click(self, ev):
        if ev.xdata is None or ev.ydata is None:
            return
        j, i = int(round(ev.xdata)), int(round(ev.ydata))
        if 0 <= j < self.x.count() and 0 <= i < self.y.count():
            self.x.blockSignals(True); self.x.setCurrentIndex(j); self.x.blockSignals(False)
            self.y.setCurrentIndex(i)
            self.draw_scatter()

    @staticmethod
    def values(rows, key, log_cb):
        v = stats.numeric(rows, key)
        if log_cb.isChecked():
            with np.errstate(divide='ignore', invalid='ignore'):
                v = np.where(v > 0, np.log10(v), np.nan)  # values <= 0 have no log and are left out
        return v

    def color_key(self):
        if not self.pooled:
            return 'cell'
        return self.color.currentData()

    @safe_draw
    def draw_scatter(self, *_):
        xk, yk = self.x.currentData(), self.y.currentData()
        f = self.sc_fig; f.clear(); ax = f.add_subplot(111)
        rows = self.rows()
        if not xk or not yk or not rows:
            self.scatter.draw_idle(); self.stat.setText(''); return
        x, y = self.values(rows, xk, self.logx), self.values(rows, yk, self.logy)
        xl = ('log10 ' if self.logx.isChecked() else '') + self.x.currentText()
        yl = ('log10 ' if self.logy.isChecked() else '') + self.y.currentText()
        ok = np.isfinite(x) & np.isfinite(y)
        key = self.color_key()
        cats = sorted({str(r.get(key, '')) for r in rows}) if key else ['']
        many = len(cats) > 12
        small = self.level.currentData() == 'mito' or len(rows) > 200
        for n, c in enumerate(cats):
            m = ok & (np.array([str(r.get(key, '')) == c for r in rows]) if key else ok)
            color = GREEN_COLORS.get(c) if key == 'green_status' else CELL_COLORS[n % len(CELL_COLORS)]
            ax.scatter(x[m], y[m], s=12 if small else 40, alpha=0.6 if small else 0.8, color=color or '#888888',
                       label=None if many or not key else {'positive': 'green+', 'negative': 'green−'}.get(c, c),
                       edgecolors='none')
            if not self.pooled and self.level.currentData() == 'cell':
                for xi, yi in zip(x[m], y[m]):
                    ax.annotate(c[4:], (xi, yi), textcoords='offset points', xytext=(4, 4), fontsize=8)
        reg = stats.regression(x, y)
        rho, p_s, _ = stats.spearman(x, y)
        if ok.sum():
            med = self.center.currentData() == 'median'
            cx = np.median(x[ok]) if med else np.mean(x[ok])
            cy = np.median(y[ok]) if med else np.mean(y[ok])
            ax.axvline(cx, color='#555555', lw=0.9, ls='--'); ax.axhline(cy, color='#555555', lw=0.9, ls='--')
            q = {'I': (x > cx) & (y > cy), 'II': (x < cx) & (y > cy), 'III': (x < cx) & (y < cy), 'IV': (x > cx) & (y < cy)}
            tot = max(int(ok.sum()), 1)
            for name, (hx, hy, ha, va) in {'I': (0.98, 0.98, 'right', 'top'), 'II': (0.02, 0.98, 'left', 'top'),
                                           'III': (0.02, 0.02, 'left', 'bottom'), 'IV': (0.98, 0.02, 'right', 'bottom')}.items():
                k_ = int((q[name] & ok).sum())
                ax.text(hx, hy, f'{name}: {k_} ({100 * k_ / tot:.0f}%)', transform=ax.transAxes, ha=ha, va=va,
                        fontsize=8, color='#333333', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=1.5))
        if np.isfinite(reg['slope']):
            xs_ = np.linspace(np.nanmin(x[ok]), np.nanmax(x[ok]), 100)
            ys_ = reg['slope'] * xs_ + reg['intercept']
            # 95 % confidence band of the fitted line
            n_ = reg['n']; xm = x[ok].mean()
            resid = y[ok] - (reg['slope'] * x[ok] + reg['intercept'])
            se = np.sqrt(np.sum(resid ** 2) / (n_ - 2)) * np.sqrt(1 / n_ + (xs_ - xm) ** 2 / np.sum((x[ok] - xm) ** 2))
            tq = scipy_t.ppf(0.975, n_ - 2)
            ax.fill_between(xs_, ys_ - tq * se, ys_ + tq * se, color='#E15759', alpha=0.15, lw=0)
            ax.plot(xs_, ys_, color='#E15759', lw=1.8)
        ax.set_xlabel(xl); ax.set_ylabel(yl)
        if not many and key:
            ax.legend(fontsize=7, frameon=False, ncol=2, loc='upper center')
        ax.spines[['top', 'right']].set_visible(False)
        why = stats.why_nan(x, y, self.x.currentText(), self.y.currentText())
        if why:
            unit = 'cells' if self.level.currentData() == 'cell' else 'mito objects'
            ax.set_title(f'regression not computable: {why.replace("values", unit)}', fontsize=9, color='#B03A2E')
            tight(f); self.scatter.draw_idle()
            self.stat.setText(f'No regression / correlation for {yl} vs {xl}: {why.replace("values", unit)}. '
                              + ('Cell-level statistics need at least 3 cells in the current selection; choose '
                                 '"All" cells or a wider dataset filter, or analyse more images of this dataset.'
                                 if reg['n'] < 3 else 'A metric that does not vary has no correlation.'))
            return
        ax.set_title(f"r = {reg['r']:.3f}, R² = {reg['r2']:.3f}, p = {reg['p']:.2g}, n = {reg['n']}"
                     f"   (Spearman ρ = {rho:.3f})", fontsize=9)
        tight(f); self.scatter.draw_idle()
        self.stat.setText(f"{yl} = {reg['slope']:.4g} × {xl} + {reg['intercept']:.4g};  Pearson r = {reg['r']:.3f}, "
                          f"R² = {reg['r2']:.3f}, p = {reg['p']:.3g}, n = {reg['n']};  Spearman ρ = {rho:.3f}, "
                          f"p = {p_s:.3g}")

    def export(self):
        if not (self.data['cell'] or self.data['mito']):
            return
        d = QFileDialog.getExistingDirectory(self, 'Export pooled tables to', getattr(self, 'export_dir', ''))
        if not d:
            return
        tag = self.group.currentData() or 'all'
        tag = {'positive': 'green_pos', 'negative': 'green_neg'}.get(tag, tag)
        cells, mito = self.base_rows('cell'), self.base_rows('mito')
        pipeline.write_csv(os.path.join(d, f'pooled_cells_{tag}.csv'), cells)
        pipeline.write_csv(os.path.join(d, f'pooled_mito_{tag}.csv'), mito)
        reg = stats.correlation_table(cells, mito)
        pipeline.write_csv(os.path.join(d, f'pooled_regression_{tag}.csv'), reg)
        pipeline.write_xlsx(os.path.join(d, f'pooled_{tag}.xlsx'),
                            [('regression', reg), ('cells', cells), ('mito', mito)])
        QMessageBox.information(self, 'Exported', f'{len(cells)} cells, {len(mito)} mito objects and the regression '
                                f'table written to {d} (pooled_*_{tag}.*)')


def group_combo():
    c = QComboBox(); c.setSizeAdjustPolicy(QComboBox.AdjustToContents)
    c.setMaximumWidth(240)   # long group names must not widen the window; the list itself shows them in full
    c.view().setMinimumWidth(240); c.view().setTextElideMode(Qt.ElideNone)
    return c


def fmt_p(p):
    """p value for display; values that underflow to 0 are shown as < 1e-300."""
    return 'n/a' if not np.isfinite(p) else ('< 1e-300' if p < 1e-300 else f'{p:.2g}')


class GroupPicker(QPushButton):
    """Button with a checkable menu of groups (which groups to show); emits `changed`."""
    changed = Signal()

    def __init__(self, text='Groups shown'):
        super().__init__(text)
        from PySide6.QtWidgets import QMenu
        self.menu = QMenu(self); self.setMenu(self.menu); self.groups = []; self.label = text

    def set_groups(self, groups):
        old = {a.text(): a.isChecked() for a in self.menu.actions() if a.isCheckable()}
        self.menu.clear(); self.groups = list(groups)
        a_all = self.menu.addAction('Show all'); a_all.triggered.connect(lambda: self._set_all(True))
        a_none = self.menu.addAction('Hide all'); a_none.triggered.connect(lambda: self._set_all(False))
        self.menu.addSeparator()
        for g in self.groups:
            a = QAction(g, self.menu); a.setCheckable(True); a.setChecked(old.get(g, True))
            a.toggled.connect(self._toggled); self.menu.addAction(a)
        self._update_text()

    def _set_all(self, on):
        for a in self.menu.actions():
            if a.isCheckable():
                a.blockSignals(True); a.setChecked(on); a.blockSignals(False)
        self._toggled()

    def _toggled(self, *_):
        self._update_text(); self.changed.emit()

    def _update_text(self):
        sel = self.selected()
        self.setText(f'{self.label}: {len(sel)}/{len(self.groups)}')

    def selected(self):
        return [a.text() for a in self.menu.actions() if a.isCheckable() and a.isChecked()]


class GroupCorrTab(RedrawWhenShown, QWidget):
    """Correlation heatmaps per group.

    A vs B: two groups side by side and their difference (B − A) with a Fisher z test per square
    (* p < 0.05, ** p < 0.01); ◀ / ▶ step group B through the other groups.
    All groups: one heatmap per shown group (small multiples) and a map of Cochran's Q test that the
    correlation is the same in every shown group."""

    def __init__(self):
        super().__init__()
        self.data = {'cell': [], 'mito': []}
        self.cache = {}
        self.mode = QComboBox(); self.mode.addItem('A vs B', 'pair'); self.mode.addItem('All groups', 'all')
        self.level = QComboBox(); self.level.addItem('Per cell', 'cell'); self.level.addItem('Per mito object', 'mito')
        self.green = green_group_combo()
        self.a, self.b = group_combo(), group_combo()
        self.prev_b, self.next_b = QPushButton('◀'), QPushButton('▶')
        for btn in (self.prev_b, self.next_b):
            btn.setFixedWidth(30); btn.setToolTip('Step group B through the groups')
        self.picker = GroupPicker()
        self.hstat = QComboBox(); self.hstat.addItem('Pearson r', 'r'); self.hstat.addItem('Spearman ρ', 'rho')
        self.export_btn = QPushButton('Export…')
        row = QHBoxLayout()
        for w in (QLabel('Mode'), self.mode, QLabel('Level'), self.level, QLabel('Cells'), self.green,
                  QLabel('Group A'), self.a, QLabel('Group B'), self.prev_b, self.b, self.next_b, self.picker,
                  QLabel('Statistic'), self.hstat):
            row.addWidget(w)
        row.addStretch(1); row.addWidget(self.export_btn)
        self.fig = Figure(figsize=(8, 6)); self.canvas = FigureCanvasQTAgg(self.fig)
        self.canvas.mpl_connect('motion_notify_event', self.on_hover)
        self.top = QLabel(); self.top.setTextInteractionFlags(Qt.TextSelectableByMouse); self.top.setWordWrap(True)
        self.top.setStyleSheet('font-family: Menlo, Consolas, monospace; font-size: 11px;')
        self.top.setMaximumHeight(130)
        lay = QVBoxLayout(self); lay.addLayout(row); lay.addWidget(AspectBox(self.canvas), 1); lay.addWidget(self.top)
        note = QLabel('Each heatmap is computed from the cells (or mito objects) of one group. A vs B: Δ = B − A, '
                      '* p < 0.05, ** p < 0.01 (Fisher z test for two independent correlations). All groups: '
                      "Cochran's Q test that the correlation is equal in all shown groups. Hover a square for the "
                      'values. Many squares are tested at once: about 5 % reach p < 0.05 by chance.')
        note.setWordWrap(True); lay.addWidget(note)
        for w in (self.mode, self.level, self.green, self.a, self.b, self.hstat):
            w.currentIndexChanged.connect(self.draw)
        self.picker.changed.connect(self.draw)
        self.prev_b.clicked.connect(lambda: self.step_b(-1)); self.next_b.clicked.connect(lambda: self.step_b(1))
        self.export_btn.clicked.connect(self.export)
        self.res = None

    def set_data(self, cell_rows, mito_rows):
        self.data = {'cell': cell_rows, 'mito': mito_rows}; self.cache = {}
        groups = sorted({str(r.get('group', '')) for r in cell_rows + mito_rows})
        for n, c in enumerate((self.a, self.b)):
            cur = c.currentText()
            c.blockSignals(True); c.clear(); c.addItems(groups)
            c.setCurrentIndex(groups.index(cur) if cur in groups else min(n, len(groups) - 1))
            c.blockSignals(False)
        self.picker.blockSignals(True); self.picker.set_groups(groups); self.picker.blockSignals(False)
        self.draw()

    def step_b(self, d):
        n = self.b.count()
        if n < 2:
            return
        i = self.b.currentIndex()
        for _ in range(n):
            i = (i + d) % n
            if self.b.itemText(i) != self.a.currentText():
                break
        self.b.setCurrentIndex(i)

    def axes_metrics(self):
        lv = self.level.currentData()
        rows = filter_group(self.data[lv], self.green.currentData())
        xs, ys = (stats.CELL_X, stats.CELL_Y) if lv == 'cell' else (stats.MITO_X, stats.MITO_Y)
        return rows, stats.available(xs, rows), stats.available(ys, rows)

    def matrix(self, g):
        """Correlation matrix of one group, cached per level / green filter / statistic."""
        key = (self.level.currentData(), self.green.currentData(), self.hstat.currentData(), g)
        if key not in self.cache:
            rows, xs, ys = self.axes_metrics()
            rg = [r for r in rows if str(r.get('group', '')) == g]
            v, n = stats.corr_matrix(rg, xs, ys, self.hstat.currentData())
            self.cache[key] = (v, n, len(rg))
        return self.cache[key]

    def compute(self):
        _, xs, ys = self.axes_metrics()
        ga, gb = self.a.currentText(), self.b.currentText()
        va, na, nra = self.matrix(ga); vb, nb, nrb = self.matrix(gb)
        return dict(mode='pair', xs=xs, ys=ys, ga=ga, gb=gb, va=va, vb=vb, na=na, nb=nb, d=vb - va,
                    p=stats.corr_diff_p(va, na, vb, nb), nra=nra, nrb=nrb)

    def compute_all(self):
        _, xs, ys = self.axes_metrics()
        gs = self.picker.selected()
        mats = [self.matrix(g) for g in gs]
        vals = np.array([m[0] for m in mats]) if mats else np.zeros((0, len(ys), len(xs)))
        ns = np.array([m[1] for m in mats]) if mats else np.zeros((0, len(ys), len(xs)))
        q = stats.corr_heterogeneity_p(vals, ns) if len(gs) >= 2 else np.full((len(ys), len(xs)), np.nan)
        return dict(mode='all', xs=xs, ys=ys, groups=gs, vals=vals, ns=ns, nrows=[m[2] for m in mats], q=q)

    def redraw(self):
        self.draw()

    @safe_draw
    def draw(self, *_):
        f = self.fig; f.clear()
        pair = self.mode.currentData() == 'pair'
        for w in (self.a, self.b, self.prev_b, self.next_b):
            w.setEnabled(pair)
        self.picker.setEnabled(not pair)
        if not (self.data['cell'] or self.data['mito']) or not self.a.count():
            self.canvas.draw_idle(); return
        (self.draw_pair if pair else self.draw_all)(f)

    def draw_pair(self, f):
        self.res = R = self.compute()
        name = self.hstat.currentText()
        unit = 'cells' if self.level.currentData() == 'cell' else 'mito objects'
        grid = f.subplots(2, 2)            # 4:3 canvas: A | B on top, Δ below A (the fourth cell holds the key)
        axs = [grid[0, 0], grid[0, 1], grid[1, 0]]
        dmax = max(0.5, float(np.nanmax(np.abs(R['d']))) if np.isfinite(R['d']).any() else 0.5)
        for k, (ax, v, title, cmap, lim) in enumerate((
                (axs[0], R['va'], f"A: {short(R['ga'], 24)} (n = {R['nra']} {unit})", 'RdBu_r', 1),
                (axs[1], R['vb'], f"B: {short(R['gb'], 24)} (n = {R['nrb']} {unit})", 'RdBu_r', 1),
                (axs[2], R['d'], f'Δ {name} (B − A)', 'PuOr_r', dmax))):
            im = ax.imshow(v, cmap=cmap, vmin=-lim, vmax=lim, aspect='auto')
            ax.set_xticks(range(len(R['xs'])), [short(l, 22) for _, l in R['xs']] if k == 2 else [],
                          rotation=60, ha='right', fontsize=5)
            ax.set_yticks(range(len(R['ys'])), [short(l, 24) for _, l in R['ys']] if k != 1 else [], fontsize=5)
            ax.set_title(title, fontsize=8)
            f.colorbar(im, ax=ax, fraction=0.05, pad=0.02).ax.tick_params(labelsize=6)
        key = grid[1, 1]; key.set_axis_off()
        key.text(0.02, 0.95, 'Columns (X, mito metrics):\n' + '\n'.join(short(l, 34) for _, l in R['xs']),
                 fontsize=5, va='top', transform=key.transAxes)
        for i in range(R['d'].shape[0]):
            for j in range(R['d'].shape[1]):
                p = R['p'][i, j]
                if np.isfinite(p) and p < 0.05:
                    axs[2].text(j, i, '**' if p < 0.01 else '*', ha='center', va='center', fontsize=8, color='black')
        tight(f); self.canvas.draw_idle()
        order = [(i, j) for i in range(R['d'].shape[0]) for j in range(R['d'].shape[1]) if np.isfinite(R['d'][i, j])]
        order.sort(key=lambda ij: -abs(R['d'][ij]))
        lines = [f'Largest differences in {name} (B − A):']
        for i, j in order[:8]:
            p = R['p'][i, j]
            lines.append(f"  {R['ys'][i][1]} vs {R['xs'][j][1]}:  A {R['va'][i, j]:+.2f} (n={R['na'][i, j]})  "
                         f"B {R['vb'][i, j]:+.2f} (n={R['nb'][i, j]})  Δ {R['d'][i, j]:+.2f}  "
                         f"p {fmt_p(p) if fmt_p(p).startswith('<') else '= ' + fmt_p(p)}"
                         + (' *' if np.isfinite(p) and p < 0.05 else ''))
        if not order:
            lines.append('  none computable: each group needs at least 3 ' + unit)
        n_sig = int(np.nansum(R['p'] < 0.05)); n_tot = int(np.isfinite(R['p']).sum())
        lines.append(f'  {n_sig} of {n_tot} squares at p < 0.05 (about {0.05 * n_tot:.0f} expected by chance)')
        if self.level.currentData() == 'mito':
            lines.append('  Mito objects of one cell are not independent: these p values are optimistic.')
        self.top.setText('\n'.join(lines))

    def draw_all(self, f):
        self.res = R = self.compute_all()
        gs = R['groups']; name = self.hstat.currentText()
        if not gs:
            self.top.setText('No group selected (Groups shown).'); self.canvas.draw_idle(); return
        n = len(gs) + 1
        cols = max(1, int(np.ceil(np.sqrt(n * 4 / 3))))   # grid close to the 4:3 canvas
        rows_ = int(np.ceil(n / cols)); small = n > 8
        axs = np.atleast_1d(f.subplots(rows_, cols, squeeze=False)).ravel()
        for k, ax in enumerate(axs):
            ax.set_xticks([]); ax.set_yticks([])
            if k >= n:
                ax.set_axis_off(); continue
            if k < len(gs):
                im = ax.imshow(R['vals'][k], cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
                ax.set_title(f"{short_groups(gs, 16 if small else 26)[gs[k]]} (n={R['nrows'][k]})", fontsize=5 if small else 8)
            else:
                q = R['q']
                ax.imshow(-np.log10(np.clip(q, 1e-12, 1)), cmap='Greys', vmin=0, vmax=4, aspect='auto')
                for i in range(q.shape[0]):
                    for j in range(q.shape[1]):
                        if np.isfinite(q[i, j]) and q[i, j] < 0.05:
                            ax.text(j, i, '**' if q[i, j] < 0.01 else '*', ha='center', va='center',
                                    fontsize=6 if small else 8, color='#d62728')
                ax.set_title("Q test (p)" if small else
                             "Q test: differs between groups\n(dark = small p; * p<0.05, ** p<0.01)",
                             fontsize=5 if small else 8)
            if not small:
                if k % cols == 0:
                    ax.set_yticks(range(len(R['ys'])), [l for _, l in R['ys']], fontsize=5)
                if k >= n - cols:
                    ax.set_xticks(range(len(R['xs'])), [l for _, l in R['xs']], rotation=70, ha='right', fontsize=5)
        f.colorbar(im, ax=list(axs[:len(gs)]), fraction=0.02, pad=0.01, label=name)
        q = R['q']; n_sig = int(np.nansum(q < 0.05)); n_tot = int(np.isfinite(q).sum())
        order = sorted([(q[i, j], i, j) for i in range(q.shape[0]) for j in range(q.shape[1]) if np.isfinite(q[i, j])])
        lines = [f"{len(gs)} groups. Pairs whose {name} differs most between groups (Cochran's Q):"]
        for p, i, j in order[:6]:
            sg = short_groups(gs, 14)
            per = '  '.join(f'{sg[g]} {R["vals"][k][i, j]:+.2f}' for k, g in enumerate(gs[:6])) + (' …' if len(gs) > 6 else '')
            lines.append(f"  {R['ys'][i][1]} vs {R['xs'][j][1]}: p {fmt_p(p) if fmt_p(p).startswith('<') else '= ' + fmt_p(p)}   {per}")
        lines.append(f'  {n_sig} of {n_tot} squares at p < 0.05 (about {0.05 * n_tot:.0f} expected by chance)')
        self.top.setText('\n'.join(lines))
        self.canvas.draw_idle()

    def on_hover(self, ev):
        R = self.res
        if R is None or ev.inaxes is None or ev.xdata is None or not ev.inaxes.images:   # heatmaps only
            return
        j, i = int(round(ev.xdata)), int(round(ev.ydata))
        if R['mode'] == 'pair':
            if 0 <= i < R['d'].shape[0] and 0 <= j < R['d'].shape[1]:
                self.canvas.setToolTip(f"{R['ys'][i][1]} vs {R['xs'][j][1]}\nA ({R['ga']}): {R['va'][i, j]:.3f}, "
                                       f"n = {R['na'][i, j]}\nB ({R['gb']}): {R['vb'][i, j]:.3f}, n = {R['nb'][i, j]}\n"
                                       f"Δ = {R['d'][i, j]:.3f}, p = {fmt_p(R['p'][i, j])}")
        elif 0 <= i < R['q'].shape[0] and 0 <= j < R['q'].shape[1]:
            per = '\n'.join(f"{g}: {R['vals'][k][i, j]:.3f} (n = {R['ns'][k][i, j]})" for k, g in enumerate(R['groups']))
            self.canvas.setToolTip(f"{R['ys'][i][1]} vs {R['xs'][j][1]}\n{per}\nQ test p = {fmt_p(R['q'][i, j])}")

    def export(self):
        R = self.res
        if R is None:
            return
        st = self.hstat.currentData()
        if R['mode'] == 'pair':
            fname = f"group_corr_{R['ga']}_vs_{R['gb']}.csv"
            rows = [dict(y=R['ys'][i][0], x=R['xs'][j][0], **{f"{st}_A ({R['ga']})": R['va'][i, j], 'n_A': R['na'][i, j],
                                                              f"{st}_B ({R['gb']})": R['vb'][i, j], 'n_B': R['nb'][i, j],
                                                              'delta_B_minus_A': R['d'][i, j], 'p_fisher_z': R['p'][i, j]})
                    for i in range(len(R['ys'])) for j in range(len(R['xs']))]
        else:
            fname = 'group_corr_all_groups.csv'
            rows = []
            for i in range(len(R['ys'])):
                for j in range(len(R['xs'])):
                    d = dict(y=R['ys'][i][0], x=R['xs'][j][0])
                    for k, g in enumerate(R['groups']):
                        d[f'{st} ({g})'] = R['vals'][k][i, j]; d[f'n ({g})'] = R['ns'][k][i, j]
                    d['p_cochran_q'] = R['q'][i, j]
                    rows.append(d)
        path, _ = QFileDialog.getSaveFileName(self, 'Export group correlation comparison',
                                              os.path.join(getattr(self, 'export_dir', ''), fname), 'CSV (*.csv)')
        if path:
            pipeline.write_csv(path, rows)


ERRORS = (('mean ± SD', 'mean', 'sd'), ('mean ± SEM', 'mean', 'sem'), ('mean ± 95% CI', 'mean', 'ci95'),
          ('median ± IQR', 'median', 'iqr'))


class GroupStatsTab(RedrawWhenShown, QWidget):
    """Every metric compared across groups: mean / median with error bars over the individual values, and a
    long table (one row per metric x group) with n, mean, SD, SEM, 95 % CI, median, IQR, the across-group tests
    and, when a control group is chosen, each group vs the control (Welch, Mann-Whitney, Holm-adjusted)."""

    def __init__(self):
        super().__init__()
        self.data = {'cell': [], 'mito': []}
        self.level = QComboBox(); self.level.addItem('Per cell', 'cell'); self.level.addItem('Per mito object', 'mito')
        self.green = green_group_combo()
        self.unit = QComboBox()
        self.unit.addItem('cells / objects', 'row'); self.unit.addItem('images (mean per image)', 'image')
        self.unit.setToolTip('Images: each image counts once (the mean of its cells), so a group with many cells '
                             'from few images is not over-weighted')
        self.metric = QComboBox(); self.metric.setSizeAdjustPolicy(QComboBox.AdjustToContents)
        self.err = QComboBox()
        for label, *_ in ERRORS:
            self.err.addItem(label)
        self.control = QComboBox(); self.control.setToolTip('Compare every group with this one (Holm-adjusted)')
        self.picker = GroupPicker()
        self.table_metric = QCheckBox('Table: this metric only')
        self.export_btn = QPushButton('Export…')
        row1, row2 = QHBoxLayout(), QHBoxLayout()
        for w in (QLabel('Level'), self.level, QLabel('Cells'), self.green, QLabel('Unit'), self.unit,
                  QLabel('Metric'), self.metric):
            row1.addWidget(w)
        for w in (QLabel('Show'), self.err, QLabel('Control'), self.control, self.picker, self.table_metric):
            row2.addWidget(w)
        row1.addStretch(1); row2.addStretch(1); row2.addWidget(self.export_btn)
        row = QVBoxLayout(); row.addLayout(row1); row.addLayout(row2)
        self.fig = Figure(figsize=(8, 6)); self.canvas = FigureCanvasQTAgg(self.fig)
        self.table = TableTab(cell_filter=False)
        self.table.view.clicked.connect(self.row_clicked)
        split = QSplitter(Qt.Vertical); split.addWidget(AspectBox(self.canvas)); split.addWidget(self.table)
        split.setSizes([560, 260])
        lay = QVBoxLayout(self); lay.addLayout(row); lay.addWidget(split, 1)
        note = QLabel('Points = individual cells / objects (or image means); bar = the chosen centre and error. '
                      'All groups: Welch t and Mann-Whitney U for two groups, one-way ANOVA and Kruskal-Wallis for '
                      'more. With a control group, stars mark groups that differ from it (Mann-Whitney, Holm-adjusted: '
                      '* p < 0.05, ** p < 0.01, *** p < 0.001). Click a table row to plot that metric.')
        note.setWordWrap(True); lay.addWidget(note)
        self.level.currentIndexChanged.connect(self.level_changed)
        for w in (self.green, self.unit, self.control):
            w.currentIndexChanged.connect(self.refresh)
        self.picker.changed.connect(self.refresh)
        self.table_metric.toggled.connect(self.fill_table)
        for w in (self.metric, self.err):
            w.currentIndexChanged.connect(self.draw)
        self.metric.currentIndexChanged.connect(lambda *_: self.table_metric.isChecked() and self.fill_table())
        from PySide6.QtWidgets import QMenu
        menu = QMenu(self)
        menu.addAction('Plotted values (CSV)…', self.export_plotted).setToolTip(
            'Only the dots in the plot: the shown groups and the plotted metric')
        menu.addAction('Statistics table (CSV)…', self.export)
        self.export_btn.setMenu(menu)
        self.summary = []

    def set_data(self, cell_rows, mito_rows):
        self.data = {'cell': cell_rows, 'mito': mito_rows}
        groups = sorted({str(r.get('group', '')) for r in cell_rows + mito_rows})
        cur = self.control.currentText()
        self.control.blockSignals(True); self.control.clear(); self.control.addItem('(none)', '')
        for g in groups:
            self.control.addItem(g, g)
        i = self.control.findText(cur); self.control.setCurrentIndex(max(i, 0)); self.control.blockSignals(False)
        self.picker.blockSignals(True); self.picker.set_groups(groups); self.picker.blockSignals(False)
        self.level_changed()

    def metrics(self):
        lv = self.level.currentData()
        xs, ys = (stats.CELL_X, stats.CELL_Y) if lv == 'cell' else (stats.MITO_X, stats.MITO_Y)
        return stats.available(xs + ys, self.data[lv])

    def level_changed(self, *_):
        cur = self.metric.currentData()
        self.metric.blockSignals(True); self.metric.clear()
        for k, label in self.metrics():
            self.metric.addItem(label, k)
        i = self.metric.findData(cur)
        self.metric.setCurrentIndex(max(i, 0)); self.metric.blockSignals(False)
        self.refresh()

    def rows(self):
        rows = filter_group(self.data[self.level.currentData()], self.green.currentData())
        shown = set(self.groups())
        rows = [r for r in rows if str(r.get('group', '')) in shown]
        if self.unit.currentData() == 'image':
            rows = stats.sample_means(rows, [k for k, _ in self.metrics()])
        return rows

    def groups(self):
        return self.picker.selected()

    def refresh(self, *_):
        ctrl = self.control.currentData() or None
        self.summary = stats.group_summary_long(self.rows(), self.metrics(), self.groups(),
                                                ctrl if ctrl in self.groups() else None)
        self.fill_table()
        self.draw()

    def fill_table(self, *_):
        rows = self.summary
        if self.table_metric.isChecked():
            rows = [r for r in rows if r.get('column') == self.metric.currentData()]
        self.table.set_rows(rows)

    def row_clicked(self, idx):
        r = self.table.proxy.mapToSource(idx).row()
        if 0 <= r < len(self.table.model.rows):
            i = self.metric.findData(self.table.model.rows[r].get('column'))
            if i >= 0:
                self.metric.setCurrentIndex(i)

    def redraw(self):
        self.draw()

    @safe_draw
    def draw(self, *_):
        f = self.fig; f.clear(); ax = f.add_subplot(111)
        mk = self.metric.currentData(); groups = self.groups()
        if not mk or not groups:
            self.canvas.draw_idle(); return
        rows = self.rows()
        _, center, err = ERRORS[self.err.currentIndex()]
        rng = np.random.default_rng(0)
        arrays = []
        many = len(groups) > 6
        for n, g in enumerate(groups):
            v = stats.numeric([r for r in rows if str(r.get('group', '')) == g], mk); v = v[np.isfinite(v)]
            arrays.append(v)
            color = CELL_COLORS[n % len(CELL_COLORS)]
            ax.scatter(n + rng.uniform(-0.18, 0.18, len(v)), v, s=6 if (len(v) > 100 or many) else 22, alpha=0.5,
                       color=color, edgecolors='none')
            d = stats.describe(v)
            if not d['n']:
                continue
            c = d[center]
            lo, hi = (c - d['q1'], d['q3'] - c) if err == 'iqr' else (d[err], d[err])
            ax.errorbar(n, c, yerr=[[np.nan_to_num(lo)], [np.nan_to_num(hi)]], fmt='_', color='black',
                        ms=14 if many else 28, mew=2.0, elinewidth=1.2, capsize=4 if many else 8, capthick=1.2, zorder=5)
        ctrl = self.control.currentData()
        if ctrl in groups:
            top = max((a.max() for a in arrays if len(a)), default=0)
            span = top - min((a.min() for a in arrays if len(a)), default=0) or 1
            for r in self.summary:
                if r.get('column') == mk and r['group'] in groups and r['group'] != ctrl:
                    p = r.get('p MWU vs control (Holm)', np.nan)
                    if np.isfinite(p) and p < 0.05:
                        ax.text(groups.index(r['group']), top + 0.04 * span,
                                '***' if p < 0.001 else '**' if p < 0.01 else '*', ha='center', fontsize=10)
            ax.axvspan(groups.index(ctrl) - 0.4, groups.index(ctrl) + 0.4, color='#eeeeee', zorder=0)
        if many:   # names shortened in the middle so rotated labels never run into each other
            ng = len(groups); sg = short_groups(groups, 24)
            ax.set_xticks(range(ng), [f'{sg[g]} (n={len(a)})' for g, a in zip(groups, arrays)],
                          rotation=60 if ng <= 12 else 90, ha='right' if ng <= 12 else 'center',
                          fontsize=7 if ng <= 12 else 6 if ng <= 24 else 5)
        else:
            sg = short_groups(groups, 18)
            ax.set_xticks(range(len(groups)), [f'{sg[g]}\nn = {len(a)}' for g, a in zip(groups, arrays)],
                          fontsize=8)
        pad = max(0.0, (4 - len(groups)) / 2)  # keep few groups close together on a wide axis
        ax.set_xlim(-0.6 - pad, len(groups) - 0.4 + pad)
        ax.set_ylabel(self.metric.currentText())
        ax.spines[['top', 'right']].set_visible(False)
        p1, p2, n1, n2 = stats.group_test(arrays)
        unit = {'row': 'cells' if self.level.currentData() == 'cell' else 'mito objects',
                'image': 'images'}[self.unit.currentData()]
        ax.set_title(f'{self.err.currentText()}, unit = {unit}, {len(groups)} groups\n{n1} p {fmt_p(p1) if fmt_p(p1).startswith("<") else "= " + fmt_p(p1)},'
                     f'  {n2} p {fmt_p(p2) if fmt_p(p2).startswith("<") else "= " + fmt_p(p2)}'
                     + (f';  * vs {short(ctrl, 24)}' if ctrl in groups else ''), fontsize=8)
        tight(f); self.canvas.draw_idle()

    def plotted_rows(self):
        """The dots of the current plot: shown groups, current level / cells filter / unit, finite values of the
        plotted metric only, with just the columns that identify each dot."""
        mk = self.metric.currentData(); groups = self.groups()
        if not mk or not groups:
            return []
        ident = ['group', 'dataset', 'preset', 'sample']
        ident += ['n_rows'] if self.unit.currentData() == 'image' else (
            ['cell', 'green_status'] if self.level.currentData() == 'cell' else ['mito', 'cell', 'green_status'])
        out, rows = [], self.rows()
        for g in groups:
            for r in rows:
                if str(r.get('group', '')) != g:
                    continue
                v = stats.numeric([r], mk)[0]
                if np.isfinite(v):
                    out.append({**{k: r.get(k, '') for k in ident}, mk: float(v)})
        return out

    def export_plotted(self):
        rows = self.plotted_rows()
        if not rows:
            QMessageBox.information(self, 'Export', 'Nothing is plotted.'); return
        tag = ('cells' if self.level.currentData() == 'cell' else 'mito') + (
            '_per_image' if self.unit.currentData() == 'image' else '')
        green = self.green.currentData()
        tag += f'_{green}' if green and green != 'all' else ''
        path, _ = QFileDialog.getSaveFileName(self, 'Export plotted values',
                                              os.path.join(getattr(self, 'export_dir', ''),
                                                           f'{self.metric.currentData()}_{tag}.csv'), 'CSV (*.csv)')
        if path:
            pipeline.write_csv(path, rows)

    def export(self):
        if not self.summary:
            return
        tag = 'cells' if self.level.currentData() == 'cell' else 'mito'
        tag += '_per_image' if self.unit.currentData() == 'image' else ''
        path, _ = QFileDialog.getSaveFileName(self, 'Export group statistics',
                                              os.path.join(getattr(self, 'export_dir', ''), f'group_stats_{tag}.csv'),
                                              'CSV (*.csv)')
        if path:
            pipeline.write_csv(path, self.summary)


class GroupEditor(QDialog):
    """Edit the group of each sample of a results folder (saved to <folder>/groups.csv)."""

    def __init__(self, samples, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Groups'); self.resize(640, 480)
        self.samples = samples
        self.table = QTableWidget(len(samples), 4)
        self.table.setHorizontalHeaderLabels(('Group', 'Dataset', 'Preset', 'Sample'))
        self.table.horizontalHeader().setStretchLastSection(True)
        for r, s in enumerate(samples):
            for c, k in enumerate(('group', 'dataset', 'preset', 'sample')):
                it = QTableWidgetItem(str(s[k]))
                if c:
                    it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(r, c, it)
        setb = QPushButton('Set group for selected…'); setb.clicked.connect(self.set_selected)
        bb = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.accept); bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel('Edit the Group column (double-click), or select rows and set one group for them.'))
        lay.addWidget(self.table, 1); lay.addWidget(setb); lay.addWidget(bb)

    def set_selected(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return
        names = sorted({self.table.item(r, 0).text() for r in range(self.table.rowCount())})
        name, ok = QInputDialog.getItem(self, 'Set group', f'Group for {len(rows)} sample(s):', names, 0, True)
        if ok and name.strip():
            for r in rows:
                self.table.item(r, 0).setText(name.strip())

    def entries(self):
        return [dict(outdir=s['folder'], dataset=s['dataset'], sample=s['sample'],
                     group=self.table.item(r, 0).text().strip() or s['dataset']) for r, s in enumerate(self.samples)]



class AnalysisTab(QWidget):
    """Pooled analysis of everything under an output folder (every <sample>_per_cell.csv / _per_mito.csv)."""

    def __init__(self):
        super().__init__()
        self.folder = QLineEdit(); self.folder.setPlaceholderText('Output folder of earlier runs (searched recursively)')
        browse = QPushButton('Browse…'); browse.clicked.connect(self.browse)
        load = QPushButton('Load'); load.clicked.connect(self.load)
        self.groups_btn = QPushButton('Groups…'); self.groups_btn.clicked.connect(self.edit_groups)
        self.groups_btn.setToolTip('Change which group each sample belongs to (saved in groups.csv of the results folder)')
        self.groups_btn.setEnabled(False)
        self.export_groups_btn = QPushButton('Export by group…')
        self.export_groups_btn.setToolTip('Write one folder per group (cells, mito objects, per-image means, '
                                          'correlations, green+ / green−) and a group comparison table')
        self.export_groups_btn.setEnabled(False); self.export_groups_btn.clicked.connect(self.export_groups)
        self.samples = []
        self.folder.returnPressed.connect(self.load)
        self.info = QLabel(''); self.info.setWordWrap(True)
        top = QHBoxLayout(); top.addWidget(QLabel('Results folder')); top.addWidget(self.folder, 1)
        top.addWidget(browse); top.addWidget(load); top.addWidget(self.groups_btn); top.addWidget(self.export_groups_btn)
        self.corr = CorrelationTab(pooled=True)
        self.cells = TableTab(); self.mito = TableTab()
        tabs = QTabWidget()
        self.gcorr = GroupCorrTab(); self.gstats = GroupStatsTab()
        tabs.addTab(self.corr, 'Regression / correlation')
        tabs.addTab(self.gcorr, 'Groups: correlation')
        tabs.addTab(self.gstats, 'Groups: mean / median')
        tabs.addTab(self.cells, 'Cells (all samples)')
        tabs.addTab(self.mito, 'Mito objects (all samples)')
        lay = QVBoxLayout(self); lay.addLayout(top); lay.addWidget(self.info); lay.addWidget(tabs, 1)

    def browse(self):
        d = QFileDialog.getExistingDirectory(self, 'Results folder', self.folder.text() or os.path.expanduser('~'))
        if d:
            self.folder.setText(d); self.load()

    def load(self, *_):
        d = os.path.expanduser(self.folder.text().strip())
        if not os.path.isdir(d):
            QMessageBox.warning(self, 'Not a folder', f'{d or "(empty)"} is not a folder.'); return
        cells, mito, samples = pipeline.load_results(d)
        if not samples:
            QMessageBox.warning(self, 'Nothing found', f'No <sample>_per_cell.csv found under {d}.'); return
        self.samples = samples; self.root = d; self.groups_btn.setEnabled(True); self.export_groups_btn.setEnabled(True)
        for t in (self.corr, self.gcorr, self.gstats):
            t.export_dir = d
            t.set_data(cells, mito)
        self.cells.set_rows(cells); self.mito.set_rows(mito)
        n_pos = sum(r.get('green_status') == 'positive' for r in cells)
        n_ds = len({r['dataset'] for r in cells}); n_pr = len({(r['dataset'], r['preset']) for r in cells})
        groups = sorted({s['group'] for s in samples})
        self.info.setToolTip('Groups: ' + ', '.join(groups))
        self.info.setText(f'{len(samples)} samples in {len(groups)} group(s) ({group_list(groups)}); '
                          f'{n_ds} dataset(s), {n_pr} preset-dataset folder(s): '
                          f'{len(cells)} cells ({n_pos} green+, {len(cells) - n_pos} green−), {len(mito)} mito objects')


    def export_groups(self):
        if not self.samples:
            return
        groups = sorted({s['group'] for s in self.samples})
        ctrl, ok = QInputDialog.getItem(self, 'Export by group', 'Control group for the "vs control" tests:',
                                        ['(none)'] + groups, 0, False)
        if not ok:
            return
        d = QFileDialog.getExistingDirectory(self, 'Export group tables to (a "groups" folder is used by default)',
                                             os.path.join(self.root, pipeline.GROUP_EXPORT_DIR)
                                             if os.path.isdir(os.path.join(self.root, pipeline.GROUP_EXPORT_DIR))
                                             else self.root)
        if not d:
            return
        out = d if os.path.basename(d) == pipeline.GROUP_EXPORT_DIR else os.path.join(d, pipeline.GROUP_EXPORT_DIR)
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            out = pipeline.write_group_exports(self.root, out, None if ctrl == '(none)' else ctrl, log=lambda s: None)
        except OSError as e:
            QMessageBox.warning(self, 'Export by group', f'Cannot write to {d}:\n{e}'); return
        finally:
            QApplication.restoreOverrideCursor()
        QMessageBox.information(self, 'Export by group', f'{len(groups)} group folder(s) and group_comparison.xlsx '
                                f'written to\n{out}')
        QDesktopServices.openUrl(QUrl.fromLocalFile(out))

    def edit_groups(self):
        if not self.samples:
            return
        dlg = GroupEditor(self.samples, self)
        if dlg.exec() == QDialog.Accepted:
            try:
                pipeline.write_groups(self.root, dlg.entries())
            except OSError as e:
                QMessageBox.warning(self, 'Groups', f'Cannot write groups.csv to {self.root}:\n{e}'); return
            self.load()


class ThresholdDialog(QDialog):
    """Live preview of one threshold on one image set: slider / value box, the automatic value of that image,
    and the mask drawn on the image (zoom and pan with the toolbar). "Apply to all images" sets the manual
    value in Options; "This image only" stores it for that row of the table."""

    def __init__(self, win, row=0):
        super().__init__(win)
        from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
        from PySide6.QtCore import QTimer
        from PySide6.QtWidgets import QSlider
        from . import thresholds
        self.win, self.T = win, thresholds
        self.setWindowTitle('Thresholds: preview and adjust'); self.resize(1100, 860)
        self.sets = QComboBox()
        for r in range(win.table.rowCount()):
            j = win.table.job(r)
            self.sets.addItem(f"{j['group']} / {j['name']}", r)
        self.sets.setCurrentIndex(min(row, self.sets.count() - 1))
        self.key = QComboBox()
        for key, label, unit, img, auto in pipeline.THRESHOLDS:
            self.key.addItem(label, key)
        self.manual = QCheckBox('Manual')
        self.slider = QSlider(Qt.Horizontal); self.slider.setRange(0, 1000)
        self.value = QDoubleSpinBox(); self.value.setDecimals(3); self.value.setRange(0, 1e6)
        self.auto_lbl = QLabel(); self.info = QLabel(); self.info.setWordWrap(True)
        self.fig = Figure(figsize=(8, 7)); self.canvas = FigureCanvasQTAgg(self.fig)
        self.ax = self.fig.add_axes([0, 0, 1, 1]); self.ax.set_axis_off()
        self.view_mode = QComboBox(); self.view_mode.addItems(['Mask on image', 'Outline on image', 'Image only'])
        top = QHBoxLayout()
        for w in (QLabel('Image set'), self.sets, QLabel('Threshold'), self.key, QLabel('Show'), self.view_mode):
            top.addWidget(w)
        top.addStretch(1)
        mid = QHBoxLayout()
        for w in (self.manual, self.slider, self.value, self.auto_lbl):
            mid.addWidget(w, 3 if w is self.slider else 0)
        self.all_btn = QPushButton('Apply to all images'); self.one_btn = QPushButton('This image only')
        self.clear_btn = QPushButton("Clear this image's value"); close = QPushButton('Close')
        self.all_btn.setToolTip('Use this manual value for every image (Options → Thresholds)')
        self.one_btn.setToolTip('Use this value only for the selected image set (Thresholds column of the table)')
        bot = QHBoxLayout()
        for b in (self.all_btn, self.one_btn, self.clear_btn):
            bot.addWidget(b)
        bot.addStretch(1); bot.addWidget(close)
        lay = QVBoxLayout(self); lay.addLayout(top); lay.addLayout(mid); lay.addWidget(self.info)
        lay.addWidget(NavigationToolbar2QT(self.canvas, self)); lay.addWidget(self.canvas, 1); lay.addLayout(bot)
        self.timer = QTimer(self); self.timer.setSingleShot(True); self.timer.setInterval(80)
        self.timer.timeout.connect(self.draw)
        self.prev = {}; self.lo, self.hi = 0.0, 1.0; self.im = None; self.zoom = None
        self.sets.currentIndexChanged.connect(self.load)
        self.key.currentIndexChanged.connect(self.key_changed)
        self.slider.valueChanged.connect(self.slider_moved)
        self.value.valueChanged.connect(self.value_changed)
        self.manual.toggled.connect(self.manual_toggled)
        self.view_mode.currentIndexChanged.connect(self.timer.start)
        self.all_btn.clicked.connect(self.apply_all); self.one_btn.clicked.connect(self.apply_one)
        self.clear_btn.clicked.connect(self.clear_one); close.clicked.connect(self.accept)
        self.load()

    def preview(self):
        r = self.sets.currentData()
        if r not in self.prev:
            j = self.win.table.job(r)
            QApplication.setOverrideCursor(Qt.WaitCursor)
            try:
                self.prev[r] = self.T.Preview(j['red'], j['green'], j['blue'], self.win.px.value() or j['px'],
                                              self.win.nuclei_method.currentData(), self.win.sens.value())
            finally:
                QApplication.restoreOverrideCursor()
        return self.prev[r]

    def load(self, *_):
        self.zoom = None
        self.key_changed()

    def current_value(self, key):
        """Value in use for this image: its own override, else the Options manual value, else 0 (auto)."""
        ov = self.win.table.overrides(self.sets.currentData())
        if key in ov:
            return ov[key], True
        cb, sp = self.win.thr[key]
        return (sp.value(), True) if cb.isChecked() else (0.0, False)

    def key_changed(self, *_):
        key = self.key.currentData()
        P = self.preview()
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self.im, auto = P.image(key)
            self.lo, self.hi = P.range(key)
        finally:
            QApplication.restoreOverrideCursor()
        self.auto = auto
        v, man = self.current_value(key)
        for w in (self.manual, self.value, self.slider):
            w.blockSignals(True)
        self.manual.setChecked(man)
        self.value.setValue(v if man else (auto if np.isfinite(auto) else 0))
        self.slider.setValue(self.to_slider(self.value.value()))
        for w in (self.manual, self.value, self.slider):
            w.blockSignals(False)
        self.value.setEnabled(man); self.slider.setEnabled(man)
        self.auto_lbl.setText(f'auto for this image: {auto:.3g}' + {
            'thr_mito': ' (reference; the automatic mito threshold is set per cell)',
            'thr_puncta': ' (≈; the analysis computes it over the analysed cells only)'}.get(key, ''))
        info = {k: (img, a) for k, _, _, img, a in pipeline.THRESHOLDS}[key]
        self.info.setText(f'Applied to: {info[0]}.   Auto: {info[1]}.   Pixels above the value are coloured.')
        self.zoom = None
        self.draw()

    def to_slider(self, v):
        return int(round(1000 * (v - self.lo) / max(self.hi - self.lo, 1e-12)))

    def slider_moved(self, s):
        self.value.blockSignals(True); self.value.setValue(self.lo + (self.hi - self.lo) * s / 1000)
        self.value.blockSignals(False); self.timer.start()

    def value_changed(self, v):
        self.slider.blockSignals(True); self.slider.setValue(self.to_slider(v)); self.slider.blockSignals(False)
        self.timer.start()

    def manual_toggled(self, on):
        self.value.setEnabled(on); self.slider.setEnabled(on)
        if not on:
            self.value.blockSignals(True); self.value.setValue(self.auto if np.isfinite(self.auto) else 0)
            self.value.blockSignals(False); self.slider.setValue(self.to_slider(self.value.value()))
        self.timer.start()

    def draw(self):
        key = self.key.currentData(); P = self.preview()
        if self.ax.images:
            self.zoom = (self.ax.get_xlim(), self.ax.get_ylim())
        v = self.value.value() if self.manual.isChecked() else 0.0
        m = P.mask(key, v)
        base = P.display(key)
        rgb = np.dstack([base, base, base])
        mode = self.view_mode.currentIndex()
        col = np.array({'thr_nuclei': (0.2, 0.5, 1.0), 'cell_fg_level': (1.0, 0.8, 0.2), 'thr_mito': (1.0, 0.25, 0.8),
                        'thr_green_bright': (0.2, 1.0, 0.3), 'thr_puncta': (1.0, 1.0, 0.1)}[key])
        if mode == 0:
            rgb[m] = rgb[m] * 0.45 + col * 0.55
        elif mode == 1:
            from skimage.segmentation import find_boundaries
            rgb[find_boundaries(m, mode='inner')] = col
        self.ax.clear(); self.ax.set_axis_off(); self.ax.imshow(rgb, interpolation='nearest')
        if self.zoom:
            self.ax.set_xlim(self.zoom[0]); self.ax.set_ylim(self.zoom[1])
        used = v if v > 0 else self.auto
        self.ax.set_title(f"{self.key.currentText()} {'manual' if v > 0 else 'auto'} = {used:.3g}: "
                          f"{100 * m.mean():.1f} % of the image above", fontsize=9, color='white',
                          backgroundcolor='black', loc='left', y=0.97)
        self.canvas.draw_idle()

    def apply_all(self):
        key = self.key.currentData()
        if not self.manual.isChecked():
            cb, sp = self.win.thr[key]; cb.setChecked(False)
            self.win.statusBar().showMessage(f'{self.key.currentText()}: automatic for all images')
            return
        self.win.set_threshold(key, self.value.value())
        self.win.statusBar().showMessage(f'{self.key.currentText()} = {self.value.value():g} for all images '
                                         '(Options → Thresholds)')

    def apply_one(self):
        r = self.sets.currentData(); key = self.key.currentData()
        ov = self.win.table.overrides(r)
        if self.manual.isChecked():
            ov[key] = self.value.value()
        else:
            ov.pop(key, None)
        self.win.table.set_overrides(r, ov)

    def clear_one(self):
        r = self.sets.currentData(); ov = self.win.table.overrides(r); ov.pop(self.key.currentData(), None)
        self.win.table.set_overrides(r, ov); self.key_changed()


class ChannelRolesDialog(QDialog):
    """Which CZI channel is mitochondria, protein (POI) and nucleus, set once for many CZI files.
    Roles are given per channel number, so they apply to every chosen file that has the same channel order."""
    CHOICES = (('red', 'Mito'), ('green', 'Protein (POI)'), ('blue', 'Nucleus'), ('', 'Not used'))

    def __init__(self, files, roles, parent=None):
        super().__init__(parent)
        self.setWindowTitle('CZI channel roles'); self.resize(620, 300)
        infos = [imgio.czi_info(f) for f in files]
        self.n_min = min(len(i['channels']) for i in infos)
        chans = max((i['channels'] for i in infos), key=len)
        names = [tuple(c['name'] for c in i['channels']) for i in infos]
        self.table = QTableWidget(len(chans), 5)
        self.table.setHorizontalHeaderLabels(('Channel', 'Name', 'Dye', 'Emission (nm)', 'Role'))
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.combos = []
        role_of = {v: k for k, v in roles.items()}
        for i, c in enumerate(chans):
            for col, t in enumerate((str(i), c['name'], c['fluor'], f"{c['emission']:g}" if c['emission'] else '')):
                it = QTableWidgetItem(t); it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                self.table.setItem(i, col, it)
            cb = QComboBox()
            for k, lab in self.CHOICES:
                cb.addItem(lab, k)
            cb.setCurrentIndex(next(n for n, (k, _) in enumerate(self.CHOICES) if k == role_of.get(i, '')))
            self.table.setCellWidget(i, 4, cb); self.combos.append(cb)
        self.table.resizeColumnsToContents()
        msg = (f'Applies to {len(files)} CZI file(s) (' + ', '.join(os.path.basename(f) for f in files[:3])
               + (', …' if len(files) > 3 else '') + ').\nChoose the role of each channel; each role must be '
               'given to exactly one channel.')
        if len(set(names)) > 1:
            msg += (f'\nNote: the channel names differ between these files ({len(set(names))} different sets); '
                    'roles are applied by channel number.')
        lab = QLabel(msg); lab.setWordWrap(True)
        self.guess_btn = QPushButton('Guess from channel names')
        self.guess_btn.clicked.connect(lambda: self.set_roles(imgio.guess_roles(chans)))
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(self.check); bb.rejected.connect(self.reject)
        lay = QVBoxLayout(self)
        lay.addWidget(lab); lay.addWidget(self.table, 1)
        row = QHBoxLayout(); row.addWidget(self.guess_btn); row.addStretch(1); row.addWidget(bb)
        lay.addLayout(row)

    def set_roles(self, roles):
        role_of = {v: k for k, v in roles.items()}
        for i, cb in enumerate(self.combos):
            cb.setCurrentIndex(cb.findData(role_of.get(i, '')))

    def roles(self):
        out = {}
        for i, cb in enumerate(self.combos):
            if cb.currentData():
                out.setdefault(cb.currentData(), []).append(i)
        return out

    def check(self):
        got = self.roles()
        bad = [lab for k, lab in self.CHOICES[:3] if len(got.get(k, [])) != 1]
        if bad:
            QMessageBox.warning(self, 'Channel roles', 'Give each of these roles to exactly one channel: '
                                + ', '.join(bad) + '.')
            return
        if max(v[0] for v in got.values()) >= self.n_min:
            QMessageBox.warning(self, 'Channel roles', f'Some of the chosen files have only {self.n_min} channels.')
            return
        self.accept()

    def result_roles(self):
        return {k: v[0] for k, v in self.roles().items()}


class MainWindow(QMainWindow):
    OUT_HINT = 'Default: <common folder of the input folders>/dataset'

    def __init__(self):
        super().__init__()
        self.setWindowTitle(f'Mito Analyzer {__version__}')
        self.resize(1400, 900)
        self.thread = self.worker = None
        self.result = None
        self.base_out = ''
        self._mode = 'analyse'; self._review = None

        # inputs: one row per image set (red/green/blue of the same field)
        inp = QGroupBox('Image sets')
        iv = QVBoxLayout(inp)
        self.table = JobTable()
        self.table.dropped.connect(self.add_paths)
        self.table.cellDoubleClicked.connect(self.cell_double_clicked)
        self.table.currentCellChanged.connect(lambda r, *_: self.show_result(self.table.result(r)))
        iv.addWidget(self.table, 1)
        self.add_btn = QPushButton('Add files…'); self.add_btn.clicked.connect(self.browse_files)
        self.addf_btn = QPushButton('Add folder…'); self.addf_btn.clicked.connect(self.browse_folder)
        self.del_btn = QPushButton('Remove'); self.del_btn.clicked.connect(self.remove_rows)
        self.clear_btn = QPushButton('Clear'); self.clear_btn.clicked.connect(self.clear_rows)
        self.group_btn = QPushButton('Set group…'); self.group_btn.clicked.connect(self.set_group)
        self.group_btn.setToolTip('Put the selected image sets in one group; groups are compared in the Analysis tab')
        self.roles_btn = QPushButton('Channel roles…'); self.roles_btn.clicked.connect(self.edit_roles)
        self.roles_btn.setToolTip('CZI files: choose which channel is mito, protein (POI) and nucleus, for the '
                                  'selected rows (or all CZI rows if none is selected) at once')
        self.czi_roles = {}  # channel names of a CZI -> roles chosen for it; used for CZIs added later
        row = QHBoxLayout()
        for b in (self.add_btn, self.addf_btn, self.group_btn, self.roles_btn, self.del_btn, self.clear_btn):
            row.addWidget(b)
        iv.addLayout(row)
        hint = QLabel('Drop TIFFs, CZIs or folders on the table; TIFFs are grouped into red/green/blue sets by name '
                      '(…Ch1_Red / …Ch2_Green / …Ch3_Blue), a CZI file is one set (its channel roles are guessed '
                      'from the channel names; change them for many files at once with Channel roles…). '
                      'Group, dataset, preset and sample can be edited; '
                      'select rows and press Set group… to group image sets; double-click a channel to change its file.')
        hint.setWordWrap(True); hint.setStyleSheet('color: gray; font-size: 11px;')
        iv.addWidget(hint)
        form = QFormLayout()
        self.outdir = QLineEdit(); self.outdir.setPlaceholderText(self.OUT_HINT)
        ob = QPushButton('Browse…'); ob.clicked.connect(self.browse_out)
        row = QHBoxLayout(); row.addWidget(self.outdir, 1); row.addWidget(ob)
        form.addRow('Output folder', row)
        form.addRow('', QLabel('Saved as <output>/<dataset>/<preset>-<dataset>/<sample>/'))
        iv.addLayout(form)

        # options
        opt = QGroupBox('Options')
        of = QFormLayout(opt)
        self.mito_method = QComboBox()
        self.mito_method.addItem('Split objects: adaptive threshold + watershed', 'split')
        self.mito_method.addItem('MiNA classic: Otsu per cell', 'otsu')
        self.mito_method.currentIndexChanged.connect(self.method_changed)
        self.preset_default = self.mito_method.currentData()
        self.mito_method.setToolTip('Split objects separates touching mitochondria only where the contact is both dark '
                                    '(< 0.75 of the dimmer peak) and narrow (shorter than the thinner width), '
                                    'with a 1-px gap, before MiNA. MiNA classic uses one Otsu threshold per cell.')
        self.excl_binuc = QCheckBox('Exclude cells that share a weak border (look binucleate)')
        self.review_rois = QCheckBox('Review and edit the cell ROIs before the analysis (pause after ROI detection)')
        self.review_rois.setChecked(True)
        self.review_rois.setToolTip('Run all first finds the cell ROIs of every image set, then opens a window where you '
                                    'check them and fix them by hand (add, subtract, draw, delete, merge). The analysis '
                                    'continues with the ROIs you confirm; they are saved with the results and offered '
                                    'again the next time you run the same images.')
        self.trim_edge = QCheckBox('Keep cells whose tip touches the border (cut the tip off at a mito-free line)')
        self.trim_edge.setChecked(True)
        self.trim_edge.setToolTip('A cell whose nucleus is well inside the image but whose ROI reaches the border '
                                  'is cut along the mito-free line (and its shape) between its body and the border, '
                                  'when such a line exists, and analysed. The cut-off part is hatched in the summary.')
        self.incl_edge = QCheckBox('Also analyse cells touching the image border')
        self.min_area = self._spin(0, 1e5, 250, 1, ' µm²')
        self.binuc_tau = self._spin(-1, 1, 0.0, 3, '', 0.01,
                                    'Relative border contrast below which two neighbouring cells count as a weak-border pair')
        self.px = self._spin(0, 10, 0.0, 4, ' µm/px', 0.001, '0 = read the pixel size from the TIFF')
        self.sens = self._spin(0.2, 5, 1.0, 2, ' ×', 0.05,
                               'Multiplies the automatic green puncta threshold: >1 keeps fewer, brighter puncta')
        self.min_mito = self._spin(0, 10, 0.05, 3, ' µm²', 0.01, 'Mito pieces smaller than this are left out of the per-mito table')
        self.nuclei_method = QComboBox()
        for key, label in nuclei.METHODS:
            self.nuclei_method.addItem(label, key)
        self.nuclei_method.setToolTip('Texture: nuclei are told from diffuse cytoplasmic blue by their chromatin '
                                      'texture, each with its own threshold; merge keeps curved nuclei whole while '
                                      'touching nuclei with a dark, narrow contact stay apart. Global Otsu is the '
                                      'old one-threshold method.')
        self.dim_nuclei = QCheckBox('Dim (out-of-focus) nuclei also get their own cell')
        self.dim_nuclei.setChecked(True)
        self.dim_nuclei.setToolTip('Faint, compact nuclei clear of the in-focus ones seed a cell too, so their '
                                   'mitochondria are not given to a neighbour')
        self.green_pos = self._spin(0, 100, 2.0, 1, ' %', 0.5,
                                    'A cell is green-positive when bright green covers at least this share of its '
                                    'cytoplasm (cell minus nucleus)')
        of.addRow('Mito segmentation', self.mito_method)
        of.addRow(self.review_rois)
        of.addRow(self.excl_binuc); of.addRow(self.trim_edge); of.addRow(self.incl_edge)
        of.addRow('Minimum cell area', self.min_area)
        of.addRow('Weak-border threshold', self.binuc_tau)
        of.addRow('Pixel size (0 = auto)', self.px)
        of.addRow('Green puncta threshold', self.sens)
        of.addRow('Minimum mito object', self.min_mito)
        of.addRow('Nucleus detection', self.nuclei_method)
        of.addRow(self.dim_nuclei)
        of.addRow('Green+ cell: bright green ≥', self.green_pos)

        # thresholds: automatic per image, or one manual value for every image (per-image values in the table)
        thr = QGroupBox('Thresholds (auto per image, or manual for all images)')
        tf = QFormLayout(thr)
        self.thr = {}
        for key, label, unit, img, auto in pipeline.THRESHOLDS:
            cb = QCheckBox('Manual'); sp = QDoubleSpinBox(); sp.setDecimals(3); sp.setRange(0, 1e6)
            sp.setValue(0.12 if key == 'cell_fg_level' else 0.0); sp.setEnabled(False); sp.setSuffix(f'  ({unit})')
            sp.setSingleStep(0.01 if key == 'cell_fg_level' else 1.0)
            cb.toggled.connect(sp.setEnabled)
            tip = f'Applied to: {img}.\nAuto: {auto}.\nManual: this value is used for every image of the batch.'
            cb.setToolTip(tip); sp.setToolTip(tip)
            row = QHBoxLayout(); row.addWidget(cb); row.addWidget(sp, 1)
            tf.addRow(label, row)
            self.thr[key] = (cb, sp)
        self.preview_btn = QPushButton('Preview / adjust thresholds…')
        self.preview_btn.setToolTip('Tune each threshold on any image set with a live preview, then apply it to '
                                    'all images or to that image only')
        self.preview_btn.clicked.connect(self.open_threshold_preview)
        tf.addRow(self.preview_btn)
        of.addRow(thr)

        self.run_btn = QPushButton('Run all'); self.run_btn.setDefault(True)
        self.run_btn.setMinimumHeight(36); self.run_btn.clicked.connect(self.start)
        self.open_btn = QPushButton('Open output folder'); self.open_btn.setEnabled(False)
        self.open_btn.clicked.connect(self.open_out)
        btns = QHBoxLayout(); btns.addWidget(self.run_btn, 2); btns.addWidget(self.open_btn, 1)

        self.logbox = QPlainTextEdit(); self.logbox.setReadOnly(True); self.logbox.setMaximumHeight(200)
        self.logbox.setPlaceholderText('Progress appears here.')

        left = QWidget(); lv = QVBoxLayout(left)
        opt_scroll = QScrollArea(); opt_scroll.setWidget(opt); opt_scroll.setWidgetResizable(True)
        opt_scroll.setFrameShape(QScrollArea.NoFrame); opt_scroll.setMaximumHeight(300); opt_scroll.setMinimumHeight(200)
        lv.addWidget(inp, 5); lv.addWidget(opt_scroll, 2); lv.addLayout(btns); lv.addWidget(QLabel('Log')); lv.addWidget(self.logbox, 1)

        # results
        self.tabs = QTabWidget()
        self.cells_tab = TableTab(cell_filter=False)
        self.mito_tab = TableTab(); self.puncta_tab = TableTab()
        self.corr_tab = CorrelationTab()
        self.views = {k: ImageView() for k in ('green_overlay', 'green_cells', 'overlay', 'cells', 'summary')}
        self.tabs.addTab(self.cells_tab, 'Cells')
        self.tabs.addTab(self.mito_tab, 'Mito objects')
        self.tabs.addTab(self.puncta_tab, 'Green puncta')
        self.tabs.addTab(self.corr_tab, 'Correlation')
        self.tabs.addTab(self.views['green_overlay'], 'Green on mito')
        self.tabs.addTab(self.views['green_cells'], 'Green+ / − cells')
        self.tabs.addTab(self.views['overlay'], 'MiNA overlay')
        self.tabs.addTab(self.views['cells'], 'Cells (zoom)')
        self.tabs.addTab(self.views['summary'], 'Cell ROIs')

        split = QSplitter(); split.addWidget(left); split.addWidget(self.tabs)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1); split.setSizes([700, 800])
        self.analysis = AnalysisTab()
        top = QTabWidget()
        top.addTab(split, 'Run'); top.addTab(self.analysis, 'Analysis (all samples)')
        self.top_tabs = top
        self.setCentralWidget(top)
        self.statusBar().showMessage('Add image sets (files or folders), then press Run all.')

    @staticmethod
    def _spin(lo, hi, val, dec, suffix='', step=None, tip=None):
        s = QDoubleSpinBox(); s.setRange(lo, hi); s.setDecimals(dec); s.setValue(val); s.setSuffix(suffix)
        if step:
            s.setSingleStep(step)
        if tip:
            s.setToolTip(tip)
        return s

    # inputs
    def last_dir(self):
        n = self.table.rowCount()
        return os.path.dirname(self.table.job(n - 1)['red']) if n else ''

    def browse_files(self):
        paths, _ = QFileDialog.getOpenFileNames(self, 'Choose TIFF images (all channels) or CZI files',
                                                self.last_dir(), TIFF_FILTER)
        self.add_paths(paths)

    def browse_folder(self):
        d = QFileDialog.getExistingDirectory(self, 'Choose a folder of TIFF / CZI images', self.last_dir())
        if d:
            self.add_paths([d])

    def add_paths(self, paths):
        if not paths:
            return
        sets, leftover = pipeline.find_sets(paths)
        for s in sets:  # CZIs with the channel names of one whose roles were set: same roles
            if 'czi' in s:
                roles = self.czi_roles.get(tuple(c['name'] for c in imgio.czi_info(s['czi'])['channels']))
                if roles:
                    s.update({k: imgio.channel_ref(s['czi'], v) for k, v in roles.items()}, roles=roles)
        have = {self.table.czi(r) or self.table.job(r)['red'] for r in range(self.table.rowCount())}
        new = [s for s in sets if (s.get('czi') or s['red']) not in have]
        for s in new:
            self.table.add_set(s, self.preset_default)
        self.update_out_hint()
        msg = f'Added {len(new)} image set(s)'
        if len(sets) > len(new):
            msg += f', {len(sets) - len(new)} already in the table'
        n_czi = sum('czi' in s for s in new)
        if n_czi:
            msg += (f' ({n_czi} CZI: channel roles guessed from the channel names or reused; check the Red / Green / Blue '
                    'columns and use Channel roles… to change them)')
        self.statusBar().showMessage(msg)
        if leftover:
            bad_czi = [f for f in leftover if imgio.is_czi(f)]
            tif = [f for f in leftover if not imgio.is_czi(f)]
            text = ''
            if tif:
                text += ('No complete red/green/blue set was found for:\n'
                         + '\n'.join(os.path.basename(f) for f in tif[:20]) + ('\n…' if len(tif) > 20 else '')
                         + '\n\nAdd them with names ending in _Red / _Green / _Blue (optionally with Ch1/Ch2/Ch3).')
            if bad_czi:
                text += ('\n\n' if text else '') + ('These CZI files could not be read or have fewer than 3 channels:\n'
                         + '\n'.join(os.path.basename(f) for f in bad_czi[:20]))
            QMessageBox.information(self, 'Some files were not added', text)

    def cell_double_clicked(self, r, c):
        key = next((k for k, col in JobTable.CH_COL.items() if col == c), None)
        if not key or self.busy():
            return
        if self.table.czi(r):
            rows = sorted({i.row() for i in self.table.selectedIndexes()} | {r})
            self.edit_roles(rows)
            return
        old = self.table.item(r, c).data(Qt.UserRole)
        path, _ = QFileDialog.getOpenFileName(self, f'Choose {key} image', os.path.dirname(old), TIFF_FILTER)
        if path:
            self.table.set_path(r, key, path)

    def edit_roles(self, rows=None):
        """Set the channel roles of the selected CZI rows (all CZI rows if none is selected) in one go."""
        if self.busy():
            return
        if not rows:
            rows = sorted({i.row() for i in self.table.selectedIndexes()}) or range(self.table.rowCount())
        rows = [r for r in rows if self.table.czi(r)]
        if not rows:
            QMessageBox.information(self, 'Channel roles', 'Channel roles are set for CZI files; add CZI files first '
                                    '(TIFF sets take their roles from the _Red / _Green / _Blue file names).')
            return
        files = [self.table.czi(r) for r in rows]
        try:
            dlg = ChannelRolesDialog(files, self.table.roles(rows[0]), self)
        except Exception as e:
            QMessageBox.warning(self, 'Channel roles', f'Cannot read the CZI files: {e}')
            return
        if dlg.exec() != QDialog.Accepted:
            return
        roles = dlg.result_roles()
        for r in rows:
            self.table.set_roles(r, roles)
        for f in files:
            self.czi_roles[tuple(c['name'] for c in imgio.czi_info(f)['channels'])] = roles
        self.statusBar().showMessage(f'Channel roles set for {len(rows)} CZI file(s): '
                                     + ', '.join(f'{lab} = ch{roles[k]}' for k, lab in imgio.ROLES))

    def remove_rows(self):
        for r in sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)
        self.update_out_hint()

    def clear_rows(self):
        self.table.setRowCount(0); self.update_out_hint()

    def default_out(self):
        jobs = [self.table.job(r) for r in range(self.table.rowCount())]
        return pipeline.default_outdir([j['root'] for j in jobs]) or (
            os.path.join(os.path.dirname(os.path.dirname(jobs[0]['red'])), 'results') if jobs else '')

    def update_out_hint(self):
        d = self.default_out()
        self.outdir.setPlaceholderText(f'Default: {d}' if d else self.OUT_HINT)

    def set_group(self):
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            QMessageBox.information(self, 'Set group', 'Select the image sets (rows) to put in one group first.')
            return
        cur = self.table.text(rows[0], JobTable.GROUP)
        names = self.table.groups()
        name, ok = QInputDialog.getItem(self, 'Set group', f'Group for {len(rows)} image set(s):', names,
                                        names.index(cur) if cur in names else 0, True)
        name = name.strip()
        if ok and name:
            for r in rows:
                self.table.item(r, JobTable.GROUP).setText(name)

    def method_changed(self):
        """Rows still showing the old default preset follow the new mito method."""
        new = self.mito_method.currentData()
        for r in range(self.table.rowCount()):
            if self.table.text(r, JobTable.PRESET) == self.preset_default:
                self.table.item(r, JobTable.PRESET).setText(new)
        self.preset_default = new

    def browse_out(self):
        d = QFileDialog.getExistingDirectory(self, 'Choose output folder', self.outdir.text())
        if d:
            self.outdir.setText(d)

    def open_out(self):
        res = self.table.result(self.table.currentRow())
        d = os.path.dirname(os.path.dirname(res['outdir'])) if res else self.base_out
        if d and os.path.isdir(d):
            QDesktopServices.openUrl(QUrl.fromLocalFile(d))

    def busy(self):
        return bool(self.thread and self.thread.isRunning())

    # run
    def start(self):
        if self.busy():  # the button reads "Stop" while a batch runs
            self.worker.stop = True
            self.run_btn.setEnabled(False); self.run_btn.setText('Stopping after this set…')
            return
        n = self.table.rowCount()
        if not n:
            QMessageBox.warning(self, 'No input', 'Add at least one red/green/blue image set.')
            return
        jobs = [self.table.job(r) for r in range(n)]
        # a CZI channel ('<file>.czi::chN') exists when its CZI file does; report each missing file once
        missing = list(dict.fromkeys(f'row {r + 1}: {os.path.basename(imgio.split_ref(j[k])[0]) or k}'
                                     for r, j in enumerate(jobs) for k in ('red', 'green', 'blue')
                                     if not os.path.isfile(imgio.split_ref(j[k])[0])))
        if missing:
            QMessageBox.warning(self, 'Missing input', 'These files do not exist:\n' + '\n'.join(missing))
            return
        if not self.px.value() and not self.ask_pixel_size(jobs):
            return
        jobs = [self.table.job(r) for r in range(n)]
        self.base_out = os.path.expanduser(self.outdir.text().strip()) or self.default_out()
        for j in jobs:
            j['preset'] = j['preset'] or self.mito_method.currentData()
            j['outdir'] = pipeline.job_outdir(self.base_out, j['dataset'], j['preset'], j['name'])
        dup = sorted({j['outdir'] for j in jobs if sum(k['outdir'] == j['outdir'] for k in jobs) > 1})
        if dup:
            QMessageBox.warning(self, 'Duplicate names', 'Several rows would write to the same folder; give them '
                                'different sample names:\n' + '\n'.join(dup))
            return
        try:
            pipeline.write_groups(self.base_out, [dict(j, sample=j['name']) for j in jobs])
        except OSError as e:
            QMessageBox.warning(self, 'Output folder', f'Cannot write to {self.base_out}:\n{e}')
            return
        params = pipeline.Params(min_area_um2=self.min_area.value(), binuc_tau=self.binuc_tau.value(),
                                 exclude_binucleate=self.excl_binuc.isChecked(),
                                 include_edge_cells=self.incl_edge.isChecked(), trim_edge_cells=self.trim_edge.isChecked(),
                                 pixel_size_um=self.px.value(),
                                 puncta_sensitivity=self.sens.value(), min_mito_area_um2=self.min_mito.value(),
                                 mito_method=self.mito_method.currentData(),
                                 green_pos_percent=self.green_pos.value(), dim_nuclei=self.dim_nuclei.isChecked(),
                                 nuclei_method=self.nuclei_method.currentData(), **self.threshold_values())
        self.logbox.clear()
        if self.review_rois.isChecked():
            saved = {i: pipeline.read_roi_edit(j['outdir'], j['name']) for i, j in enumerate(jobs)}
            saved = {i: s for i, s in saved.items() if s is not None}
            if saved and QMessageBox.question(
                    self, 'Review cell ROIs', f'ROIs you reviewed in an earlier run were found for {len(saved)} of the '
                    f'{n} image set(s). Start from them?\n(No: start from newly detected ROIs.)') == QMessageBox.Yes:
                for i, s in saved.items():
                    jobs[i]['roi_start'] = s
            self._review = dict(jobs=jobs, params=params, items={})
            self._launch(jobs, params, 'detect')
        else:
            self._launch(jobs, params, 'analyse')

    def _launch(self, jobs, params, mode):
        """Run `jobs` in the worker thread: 'detect' = cell ROIs for review, 'analyse' = the full analysis."""
        self._mode = mode
        for r in range(len(jobs)):
            self.table.set_status(r, 'queued')
        self.set_running(True)
        self.statusBar().showMessage(f'Finding cell ROIs in {len(jobs)} image set(s) for review…' if mode == 'detect'
                                     else f'Running {len(jobs)} image set(s)…')
        self.thread = QThread()
        self.worker = Worker(list(enumerate(jobs)), params, self.base_out, mode)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.logbox.appendPlainText)
        self.worker.started.connect(lambda r: (self.table.set_status(r, 'ROIs…' if mode == 'detect' else 'running…'),
                                               self.table.selectRow(r)))
        self.worker.job_done.connect(self.job_done)
        self.worker.job_failed.connect(lambda r, msg: self.table.set_status(r, 'failed', tip=msg))
        self.worker.finished.connect(self.done)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def threshold_values(self):
        """Params values of the threshold options (0 = auto; cell area level 0.12 when not manual)."""
        out = {}
        for key, (cb, sp) in self.thr.items():
            out[key] = sp.value() if cb.isChecked() else (0.12 if key == 'cell_fg_level' else 0.0)
        return out

    def set_threshold(self, key, value):
        cb, sp = self.thr[key]
        sp.setValue(value); cb.setChecked(True)

    def open_threshold_preview(self):
        if not self.table.rowCount():
            QMessageBox.information(self, 'Thresholds', 'Add image sets first; the preview runs on one of them.')
            return
        r = self.table.currentRow()
        ThresholdDialog(self, max(r, 0)).exec()

    def ask_pixel_size(self, jobs):
        """Rows whose TIFF has no pixel size: ask once for the µm per pixel and fill them in. False = cancelled."""
        rows = [r for r, j in enumerate(jobs) if not j['px']]
        if not rows:
            return True
        names = ', '.join(jobs[r]['name'] for r in rows[:6]) + (' …' if len(rows) > 6 else '')
        v, ok = QInputDialog.getDouble(
            self, 'Pixel size', f'{len(rows)} image set(s) have no pixel size in their TIFF ({names}).\n'
            'µm per pixel for these images (from the microscope / acquisition settings):',
            cell_roi_ref_px(), 0.0001, 100.0, 5)
        if not ok:
            return False
        for r in rows:
            self.table.set_px(r, v, 'entered')
        return True

    def set_running(self, on):
        for b in (self.add_btn, self.addf_btn, self.group_btn, self.del_btn, self.clear_btn):
            b.setEnabled(not on)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers if on else
                                   QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed
                                   | QAbstractItemView.SelectedClicked)
        self.run_btn.setText('Stop' if on else 'Run all'); self.run_btn.setEnabled(True)

    def job_done(self, r, res):
        if self._mode == 'detect':
            self._review['items'][r] = res
            self.table.set_status(r, f"{int(len(np.unique(res['labels'])) - 1)} ROIs to review")
            return
        self.table.set_status(r, f"{len(res['rows'])} cells", res, res['outdir'])
        self.open_btn.setEnabled(True)
        self.show_result(res)

    def done(self, ok, failed):
        if self._mode == 'detect':
            self.thread.quit(); self.thread.wait()   # the review window opens only once the worker has stopped
            self.review_and_continue(failed)
            return
        self.set_running(False)
        for r in range(self.table.rowCount()):
            if self.table.text(r, JobTable.STATUS) == 'queued':
                self.table.set_status(r, 'stopped')
        self.statusBar().showMessage(f'{ok} image set(s) done, {failed} failed. Results in {self.base_out}')
        if ok and os.path.isdir(self.base_out):  # the pooled view follows the latest batch
            self.analysis.folder.setText(self.base_out)
            self.analysis.load()
        if failed:
            QMessageBox.warning(self, 'Some sets failed', f'{failed} image set(s) failed; hover the status cell '
                                'or see the log for details.')

    def review_and_continue(self, failed):
        """After ROI detection: open the review window, then analyse with the confirmed ROIs."""
        from .roi_review import RoiReviewDialog
        rv = self._review; jobs = rv['jobs']
        rows = sorted(rv['items'])
        if self.worker.stop or not rows:
            self.set_running(False)
            for r in range(self.table.rowCount()):
                if self.table.text(r, JobTable.STATUS) in ('queued',) or 'ROIs' in self.table.text(r, JobTable.STATUS):
                    self.table.set_status(r, 'stopped')
            self.statusBar().showMessage('Stopped before the analysis.' if rows else 'No cell ROIs found.')
            return
        self.statusBar().showMessage('Review the cell ROIs, then confirm to continue the analysis.')
        dlg = RoiReviewDialog([dict(name=f"{jobs[r]['group']} / {jobs[r]['name']}", review=rv['items'][r])
                               for r in rows], self)
        if dlg.exec() != QDialog.Accepted:
            self.set_running(False)
            for r in rows:
                self.table.set_status(r, 'ROIs not confirmed')
            self.statusBar().showMessage('Analysis cancelled at the ROI review. Press Run all to start again.')
            return
        todo = []
        for r, edit in zip(rows, dlg.results()):
            jobs[r]['roi_edit'] = edit
            todo.append(jobs[r])
        if failed:
            self.logbox.appendPlainText(f'{failed} image set(s) failed at ROI detection and are skipped.')
        self._launch_rows(rows, todo, rv['params'])

    def _launch_rows(self, rows, jobs, params):
        """Analyse the reviewed rows (keeping their row numbers for the table)."""
        self._mode = 'analyse'
        for r in rows:
            self.table.set_status(r, 'queued')
        self.set_running(True)
        self.statusBar().showMessage(f'Running {len(jobs)} image set(s) with the reviewed ROIs…')
        self.thread = QThread()
        self.worker = Worker(list(zip(rows, jobs)), params, self.base_out, 'analyse')
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.logbox.appendPlainText)
        self.worker.started.connect(lambda r: (self.table.set_status(r, 'running…'), self.table.selectRow(r)))
        self.worker.job_done.connect(self.job_done)
        self.worker.job_failed.connect(lambda r, msg: self.table.set_status(r, 'failed', tip=msg))
        self.worker.finished.connect(self.done)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def show_result(self, res):
        if not res:
            return
        self.result = res
        self.cells_tab.set_rows(transpose(res['rows']))
        self.mito_tab.set_rows(res['mito_rows'])
        self.puncta_tab.set_rows(res['puncta_rows'])
        self.corr_tab.set_data(res['rows'], res['mito_rows'])
        for k, v in self.views.items():
            v.set_image(res[k])
        self.tabs.setTabText(0, f"Cells ({res['name']})")

    def closeEvent(self, e):
        if self.busy():
            if QMessageBox.question(self, 'Quit', 'An analysis is still running. Quit anyway?') != QMessageBox.Yes:
                e.ignore(); return
        e.accept()


def cell_roi_ref_px():
    return pipeline.cell_roi.REF_PX_UM


def transpose(rows):
    """Cells as columns, metrics as rows: easier to read with ~45 metrics and a handful of cells."""
    if not rows:
        return []
    head = {r['cell']: r['cell'] + {'positive': ' (+)', 'negative': ' (−)'}.get(r.get('green_status'), '') for r in rows}
    return [dict(metric=k, **{head[r['cell']]: r[k] for r in rows}) for k in rows[0] if k != 'cell']


def main():
    app = QApplication(sys.argv)
    app.setApplicationName('Mito Analyzer')
    if os.path.exists(resource('icon.png')):
        app.setWindowIcon(QIcon(resource('icon.png')))
    w = MainWindow(); w.show()
    return app.exec()
