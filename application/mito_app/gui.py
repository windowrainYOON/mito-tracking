"""PySide6 desktop GUI for the per-cell MiNA + green-on-mito pipeline."""
import math, os, sys, traceback

os.environ.setdefault('QT_API', 'pyside6')

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from scipy.stats import t as scipy_t
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QObject, QSortFilterProxyModel, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QIcon, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea,
    QSplitter, QTableView, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__, nuclei, pipeline, stats

TIFF_FILTER = 'TIFF images (*.tif *.tiff);;All files (*)'
GREEN_COLORS = {'positive': '#2CA02C', 'negative': '#7F7F7F'}
CELL_COLORS = ['#4C9BE8', '#F28E2B', '#59A14F', '#E15759', '#B07AA1', '#EDC948', '#76B7B2', '#FF9DA7',
               '#9C755F', '#BAB0AC']


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

    def __init__(self, jobs, params):
        super().__init__()
        self.jobs, self.params = jobs, params
        self.stop = False

    def run(self):
        results, failed = [], 0
        for n, (row, job) in enumerate(self.jobs, 1):
            if self.stop:
                break
            self.started.emit(row)
            self.log.emit(f'=== [{n}/{len(self.jobs)}] {job["dataset"]} / {job["name"]} ===')
            try:
                res = pipeline.run(job['red'], job['green'], job['blue'], job['outdir'], self.params, job['name'],
                                   log=self.log.emit)
                results.append(res)
                self.job_done.emit(row, res)
            except Exception as e:  # report to the user instead of crashing the app
                failed += 1
                self.log.emit(traceback.format_exc())
                self.job_failed.emit(row, f'{type(e).__name__}: {e}')
        try:
            pipeline.write_group_tables(results)
        except Exception:
            self.log.emit(traceback.format_exc())
        self.finished.emit(len(results), failed)


class JobTable(QTableWidget):
    """Input table, one image set per row. Accepts dropped TIFFs / folders."""
    COLS = ('Dataset', 'Preset', 'Sample', 'Red', 'Green', 'Blue', 'Status')
    DATASET, PRESET, NAME, RED, GREEN, BLUE, STATUS = range(7)
    CH_COL = {'red': 3, 'green': 4, 'blue': 5}
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
        for c, w in zip(range(7), (70, 50, 80, 92, 92, 92, 60)):
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
        texts = (s.get('dataset') or os.path.basename(os.path.dirname(s['red'])), preset,
                 s.get('name') or pipeline.sample_name(s['red']))
        for c, t in enumerate(texts):
            self.setItem(r, c, QTableWidgetItem(t))
        for k, c in self.CH_COL.items():
            self.set_path(r, k, s[k])
        st = QTableWidgetItem('ready'); st.setFlags(st.flags() & ~Qt.ItemIsEditable)
        self.setItem(r, self.STATUS, st)

    def set_path(self, r, key, path):
        it = QTableWidgetItem(os.path.basename(path)); it.setFlags(it.flags() & ~Qt.ItemIsEditable)
        it.setData(Qt.UserRole, path); it.setToolTip(path)
        self.setItem(r, self.CH_COL[key], it)

    def text(self, r, c):
        it = self.item(r, c)
        return it.text().strip() if it else ''

    def job(self, r):
        d = {k: self.item(r, c).data(Qt.UserRole) for k, c in self.CH_COL.items()}
        d.update(dataset=self.text(r, self.DATASET) or 'dataset', preset=self.text(r, self.PRESET),
                 name=self.text(r, self.NAME) or pipeline.sample_name(d['red']))
        return d

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


class CorrelationTab(QWidget):
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
        self.dataset = QComboBox(); self.preset = QComboBox(); self.sample = QComboBox()
        self.hstat = QComboBox(); self.hstat.addItem('Pearson r', 'r'); self.hstat.addItem('Spearman ρ', 'rho')
        self.center = QComboBox(); self.center.addItem('mean', 'mean'); self.center.addItem('median', 'median')
        self.color = QComboBox()
        for label, key in (('green status', 'green_status'), ('sample', 'sample'), ('dataset', 'dataset'),
                           ('preset', 'preset'), ('none', '')):
            self.color.addItem(label, key)
        self.logx = QCheckBox('log X'); self.logy = QCheckBox('log Y')
        row1, row2 = QHBoxLayout(), QHBoxLayout()
        for w in (QLabel('Level'), self.level, QLabel('Cells'), self.group):
            row1.addWidget(w)
        if pooled:
            for w in (QLabel('Dataset'), self.dataset, QLabel('Preset'), self.preset, QLabel('Sample'), self.sample):
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
        for w in (self.group, self.dataset, self.preset, self.sample):
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
        for combo, key in ((self.dataset, 'dataset'), (self.preset, 'preset'), (self.sample, 'sample')):
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
            for combo, key in ((self.dataset, 'dataset'), (self.preset, 'preset'), (self.sample, 'sample')):
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

    def filters_changed(self, *_):
        self.cell.blockSignals(True); self.cell.clear()
        self.cell.addItems(['All'] + sorted({r['cell'] for r in self.base_rows('mito')}))
        self.cell.setEnabled(self.level.currentData() == 'mito')
        self.cell.blockSignals(False)
        self.draw_heat(); self.draw_scatter()

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
        f.tight_layout(); self.heat.draw_idle()

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
            f.tight_layout(); self.scatter.draw_idle()
            self.stat.setText(f'No regression / correlation for {yl} vs {xl}: {why.replace("values", unit)}. '
                              + ('Cell-level statistics need at least 3 cells in the current selection; choose '
                                 '"All" cells or a wider dataset filter, or analyse more images of this dataset.'
                                 if reg['n'] < 3 else 'A metric that does not vary has no correlation.'))
            return
        ax.set_title(f"r = {reg['r']:.3f}, R² = {reg['r2']:.3f}, p = {reg['p']:.2g}, n = {reg['n']}"
                     f"   (Spearman ρ = {rho:.3f})", fontsize=9)
        f.tight_layout(); self.scatter.draw_idle()
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


class AnalysisTab(QWidget):
    """Pooled analysis of everything under an output folder (every <sample>_per_cell.csv / _per_mito.csv)."""

    def __init__(self):
        super().__init__()
        self.folder = QLineEdit(); self.folder.setPlaceholderText('Output folder of earlier runs (searched recursively)')
        browse = QPushButton('Browse…'); browse.clicked.connect(self.browse)
        load = QPushButton('Load'); load.clicked.connect(self.load)
        self.folder.returnPressed.connect(self.load)
        self.info = QLabel('')
        top = QHBoxLayout(); top.addWidget(QLabel('Results folder')); top.addWidget(self.folder, 1)
        top.addWidget(browse); top.addWidget(load)
        self.corr = CorrelationTab(pooled=True)
        self.cells = TableTab(); self.mito = TableTab()
        tabs = QTabWidget()
        tabs.addTab(self.corr, 'Regression / correlation')
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
        self.corr.export_dir = d
        self.corr.set_data(cells, mito)
        self.cells.set_rows(cells); self.mito.set_rows(mito)
        n_pos = sum(r.get('green_status') == 'positive' for r in cells)
        n_ds = len({r['dataset'] for r in cells}); n_pr = len({(r['dataset'], r['preset']) for r in cells})
        self.info.setText(f'{len(samples)} samples ({n_ds} dataset(s), {n_pr} preset-dataset folder(s)): '
                          f'{len(cells)} cells ({n_pos} green+, {len(cells) - n_pos} green−), {len(mito)} mito objects')


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f'Mito Analyzer {__version__}')
        self.resize(1400, 900)
        self.thread = self.worker = None
        self.result = None
        self.base_out = ''

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
        self.clear_btn = QPushButton('Clear'); self.clear_btn.clicked.connect(lambda: self.table.setRowCount(0))
        row = QHBoxLayout()
        for b in (self.add_btn, self.addf_btn, self.del_btn, self.clear_btn):
            row.addWidget(b)
        iv.addLayout(row)
        hint = QLabel('Drop TIFFs or folders on the table; files are grouped into red/green/blue sets by name '
                      '(…Ch1_Red / …Ch2_Green / …Ch3_Blue). Dataset, preset and sample can be edited; '
                      'double-click a channel to change its file.')
        hint.setWordWrap(True); hint.setStyleSheet('color: gray; font-size: 11px;')
        iv.addWidget(hint)
        form = QFormLayout()
        self.outdir = QLineEdit(); self.outdir.setPlaceholderText('Defaults to a "results" folder next to the dataset folder')
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
        of.addRow(self.excl_binuc); of.addRow(self.trim_edge); of.addRow(self.incl_edge)
        of.addRow('Minimum cell area', self.min_area)
        of.addRow('Weak-border threshold', self.binuc_tau)
        of.addRow('Pixel size (0 = auto)', self.px)
        of.addRow('Green puncta threshold', self.sens)
        of.addRow('Minimum mito object', self.min_mito)
        of.addRow('Nucleus detection', self.nuclei_method)
        of.addRow(self.dim_nuclei)
        of.addRow('Green+ cell: bright green ≥', self.green_pos)

        self.run_btn = QPushButton('Run all'); self.run_btn.setDefault(True)
        self.run_btn.setMinimumHeight(36); self.run_btn.clicked.connect(self.start)
        self.open_btn = QPushButton('Open output folder'); self.open_btn.setEnabled(False)
        self.open_btn.clicked.connect(self.open_out)
        btns = QHBoxLayout(); btns.addWidget(self.run_btn, 2); btns.addWidget(self.open_btn, 1)

        self.logbox = QPlainTextEdit(); self.logbox.setReadOnly(True)
        self.logbox.setPlaceholderText('Progress appears here.')

        left = QWidget(); lv = QVBoxLayout(left)
        lv.addWidget(inp, 3); lv.addWidget(opt); lv.addLayout(btns); lv.addWidget(QLabel('Log')); lv.addWidget(self.logbox, 2)

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
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1); split.setSizes([620, 880])
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
        paths, _ = QFileDialog.getOpenFileNames(self, 'Choose TIFF images (all channels)', self.last_dir(), TIFF_FILTER)
        self.add_paths(paths)

    def browse_folder(self):
        d = QFileDialog.getExistingDirectory(self, 'Choose a folder of TIFF images', self.last_dir())
        if d:
            self.add_paths([d])

    def add_paths(self, paths):
        if not paths:
            return
        sets, leftover = pipeline.find_sets(paths)
        have = {self.table.job(r)['red'] for r in range(self.table.rowCount())}
        new = [s for s in sets if s['red'] not in have]
        for s in new:
            self.table.add_set(s, self.preset_default)
        msg = f'Added {len(new)} image set(s)'
        if len(sets) > len(new):
            msg += f', {len(sets) - len(new)} already in the table'
        self.statusBar().showMessage(msg)
        if leftover:
            QMessageBox.information(self, 'Some files were not added',
                                    'No complete red/green/blue set was found for:\n'
                                    + '\n'.join(os.path.basename(f) for f in leftover[:20])
                                    + ('\n…' if len(leftover) > 20 else '')
                                    + '\n\nAdd them with names ending in _Red / _Green / _Blue '
                                      '(optionally with Ch1/Ch2/Ch3).')

    def cell_double_clicked(self, r, c):
        key = next((k for k, col in JobTable.CH_COL.items() if col == c), None)
        if not key or self.busy():
            return
        old = self.table.item(r, c).data(Qt.UserRole)
        path, _ = QFileDialog.getOpenFileName(self, f'Choose {key} image', os.path.dirname(old), TIFF_FILTER)
        if path:
            self.table.set_path(r, key, path)

    def remove_rows(self):
        for r in sorted({i.row() for i in self.table.selectedIndexes()}, reverse=True):
            self.table.removeRow(r)

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
        missing = [f'row {r + 1}: {os.path.basename(j[k]) or k}' for r, j in enumerate(jobs)
                   for k in ('red', 'green', 'blue') if not os.path.isfile(j[k])]
        if missing:
            QMessageBox.warning(self, 'Missing input', 'These files do not exist:\n' + '\n'.join(missing))
            return
        self.base_out = os.path.expanduser(self.outdir.text().strip()) or os.path.join(
            os.path.dirname(os.path.dirname(jobs[0]['red'])), 'results')
        for j in jobs:
            j['preset'] = j['preset'] or self.mito_method.currentData()
            j['outdir'] = pipeline.job_outdir(self.base_out, j['dataset'], j['preset'], j['name'])
        dup = sorted({j['outdir'] for j in jobs if sum(k['outdir'] == j['outdir'] for k in jobs) > 1})
        if dup:
            QMessageBox.warning(self, 'Duplicate names', 'Several rows would write to the same folder; give them '
                                'different sample names:\n' + '\n'.join(dup))
            return
        params = pipeline.Params(min_area_um2=self.min_area.value(), binuc_tau=self.binuc_tau.value(),
                                 exclude_binucleate=self.excl_binuc.isChecked(),
                                 include_edge_cells=self.incl_edge.isChecked(), trim_edge_cells=self.trim_edge.isChecked(),
                                 pixel_size_um=self.px.value(),
                                 puncta_sensitivity=self.sens.value(), min_mito_area_um2=self.min_mito.value(),
                                 mito_method=self.mito_method.currentData(),
                                 green_pos_percent=self.green_pos.value(), dim_nuclei=self.dim_nuclei.isChecked(),
                                 nuclei_method=self.nuclei_method.currentData())
        for r in range(n):
            self.table.set_status(r, 'queued')
        self.logbox.clear()
        self.set_running(True)
        self.statusBar().showMessage(f'Running {n} image set(s)…')
        self.thread = QThread()
        self.worker = Worker(list(enumerate(jobs)), params)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.logbox.appendPlainText)
        self.worker.started.connect(lambda r: (self.table.set_status(r, 'running…'), self.table.selectRow(r)))
        self.worker.job_done.connect(self.job_done)
        self.worker.job_failed.connect(lambda r, msg: self.table.set_status(r, 'failed', tip=msg))
        self.worker.finished.connect(self.done)
        self.worker.finished.connect(self.thread.quit)
        self.thread.start()

    def set_running(self, on):
        for b in (self.add_btn, self.addf_btn, self.del_btn, self.clear_btn):
            b.setEnabled(not on)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers if on else
                                   QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed
                                   | QAbstractItemView.SelectedClicked)
        self.run_btn.setText('Stop' if on else 'Run all'); self.run_btn.setEnabled(True)

    def job_done(self, r, res):
        self.table.set_status(r, f"{len(res['rows'])} cells", res, res['outdir'])
        self.open_btn.setEnabled(True)
        self.show_result(res)

    def done(self, ok, failed):
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
