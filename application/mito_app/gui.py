"""PySide6 desktop GUI for the per-cell MiNA + green-on-mito pipeline."""
import math, os, sys, traceback

os.environ.setdefault('QT_API', 'pyside6')

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QObject, QSortFilterProxyModel, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QIcon, QPixmap
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox,
    QHBoxLayout, QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea,
    QSplitter, QTableView, QTableWidget, QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__, pipeline, stats

TIFF_FILTER = 'TIFF images (*.tif *.tiff);;All files (*)'
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
        texts = (os.path.basename(os.path.dirname(s['red'])), preset, pipeline.sample_name(s['red']))
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
        self.count = QLabel()
        if cell_filter:
            row = QHBoxLayout(); row.addWidget(QLabel('Cell')); row.addWidget(self.cell); row.addStretch(1)
            row.addWidget(self.count); lay.addLayout(row)
            self.cell.currentTextChanged.connect(self.apply_filter)
        lay.addWidget(self.view)

    def set_rows(self, rows):
        self.model.set_rows(rows)
        self.cell.blockSignals(True); self.cell.clear()
        self.cell.addItems(['All'] + sorted({r.get('cell', '') for r in rows}))
        self.cell.blockSignals(False)
        self.proxy.setFilterKeyColumn(self.model.cols.index('cell') if 'cell' in self.model.cols else 0)
        self.apply_filter('All')
        self.view.sortByColumn(-1, Qt.AscendingOrder); self.proxy.sort(-1)  # original order until a header is clicked
        self.view.resizeColumnsToContents()

    def apply_filter(self, text):
        self.proxy.setFilterFixedString('' if text in ('All', '') else text)
        self.count.setText(f'{self.proxy.rowCount()} rows')


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
    """Spearman heatmap of all mito x green metric pairs (click a square) + scatter of the chosen pair."""

    def __init__(self):
        super().__init__()
        self.data = {'cell': [], 'mito': []}
        self.level = QComboBox(); self.level.addItem('Per cell', 'cell'); self.level.addItem('Per mito object', 'mito')
        self.x = QComboBox(); self.y = QComboBox(); self.cell = QComboBox()
        self.logx = QCheckBox('log X'); self.logy = QCheckBox('log Y')
        top = QHBoxLayout()
        for w in (QLabel('Level'), self.level, QLabel('X (mito)'), self.x, QLabel('Y (green)'), self.y,
                  QLabel('Cell'), self.cell, self.logx, self.logy):
            top.addWidget(w)
        top.addStretch(1)
        self.heat_fig = Figure(figsize=(6, 4)); self.heat = FigureCanvasQTAgg(self.heat_fig)
        self.sc_fig = Figure(figsize=(6, 4)); self.scatter = FigureCanvasQTAgg(self.sc_fig)
        self.heat.mpl_connect('button_press_event', self.on_heat_click)
        self.heat.mpl_connect('motion_notify_event', self.on_heat_hover)
        self.rho = None
        self.stat = QLabel(); self.stat.setTextInteractionFlags(Qt.TextSelectableByMouse)
        split = QSplitter(Qt.Horizontal); split.addWidget(self.heat); split.addWidget(self.scatter)
        split.setSizes([500, 500])
        lay = QVBoxLayout(self); lay.addLayout(top); lay.addWidget(split, 1); lay.addWidget(self.stat)
        lay.addWidget(QLabel('Spearman rank correlation. Mito objects in the same cell are not independent, '
                             'so per-object p values are optimistic; compare cells or images for inference.'))
        self.level.currentIndexChanged.connect(self.level_changed)
        for w in (self.x, self.y, self.cell):
            w.currentIndexChanged.connect(self.draw_scatter)
        for w in (self.logx, self.logy):
            w.toggled.connect(self.draw_scatter)
        self.cell.currentIndexChanged.connect(self.draw_heat)

    def set_data(self, cell_rows, mito_rows):
        self.data = {'cell': cell_rows, 'mito': mito_rows}
        self.level_changed()

    def metrics(self):
        lv = self.level.currentData()
        xs, ys = (stats.CELL_X, stats.CELL_Y) if lv == 'cell' else (stats.MITO_X, stats.MITO_Y)
        rows = self.data[lv]
        return stats.available(xs, rows), stats.available(ys, rows)

    def rows(self):
        rows = self.data[self.level.currentData()]
        c = self.cell.currentText()
        return rows if self.level.currentData() == 'cell' or c in ('All', '') else [r for r in rows if r['cell'] == c]

    def level_changed(self, *_):
        xs, ys = self.metrics()
        for combo, items in ((self.x, xs), (self.y, ys)):
            combo.blockSignals(True); combo.clear()
            for k, label in items:
                combo.addItem(label, k)
            combo.blockSignals(False)
        self.cell.blockSignals(True); self.cell.clear()
        self.cell.addItems(['All'] + sorted({r['cell'] for r in self.data['mito']}))
        self.cell.setEnabled(self.level.currentData() == 'mito')
        self.cell.blockSignals(False)
        mito = self.level.currentData() == 'mito'
        for cb in (self.logx, self.logy):  # object sizes span orders of magnitude
            cb.blockSignals(True); cb.setChecked(mito); cb.blockSignals(False)
        self.x.setCurrentIndex(0); self.y.setCurrentIndex(3 if mito else 1)  # length vs puncta / footprint vs green on mito
        self.draw_heat(); self.draw_scatter()

    def draw_heat(self, *_):
        xs, ys = self.metrics(); rows = self.rows()
        rho = np.full((len(ys), len(xs)), np.nan)
        for j, (xk, _) in enumerate(xs):
            for i, (yk, _) in enumerate(ys):
                rho[i, j] = stats.spearman(stats.numeric(rows, xk), stats.numeric(rows, yk))[0]
        f = self.heat_fig; f.clear(); ax = f.add_subplot(111)
        im = ax.imshow(rho, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')
        ax.set_xticks(range(len(xs)), [l for _, l in xs], rotation=60, ha='right', fontsize=7)
        ax.set_yticks(range(len(ys)), [l for _, l in ys], fontsize=7)
        self.rho = rho
        for i in range(rho.shape[0]):
            for j in range(rho.shape[1]):
                if np.isfinite(rho[i, j]) and rho.size <= 100:  # numbers only where they fit; hover shows the rest
                    ax.text(j, i, f'{rho[i, j]:.2f}', ha='center', va='center', fontsize=6,
                            color='white' if abs(rho[i, j]) > 0.6 else 'black')
        f.colorbar(im, ax=ax, fraction=0.04, label='Spearman ρ')
        ax.set_title(f'Spearman ρ, n = {len(rows)}' + (' (too few cells for a reliable ρ)' if len(rows) < 8 else '')
                     + '\nclick a square to plot it', fontsize=9)
        f.tight_layout(); self.heat.draw_idle()

    def on_heat_hover(self, ev):
        if ev.xdata is None or ev.ydata is None or self.rho is None:
            return
        j, i = int(round(ev.xdata)), int(round(ev.ydata))
        if 0 <= i < self.rho.shape[0] and 0 <= j < self.rho.shape[1]:
            self.heat.setToolTip(f'{self.y.itemText(i)} vs {self.x.itemText(j)}: ρ = {self.rho[i, j]:.3f}')

    def on_heat_click(self, ev):
        if ev.xdata is None or ev.ydata is None:
            return
        j, i = int(round(ev.xdata)), int(round(ev.ydata))
        if 0 <= j < self.x.count() and 0 <= i < self.y.count():
            self.x.blockSignals(True); self.x.setCurrentIndex(j); self.x.blockSignals(False)
            self.y.setCurrentIndex(i)
            self.draw_scatter()

    def draw_scatter(self, *_):
        xk, yk = self.x.currentData(), self.y.currentData()
        f = self.sc_fig; f.clear(); ax = f.add_subplot(111)
        rows = self.rows()
        if not xk or not yk or not rows:
            self.scatter.draw_idle(); return
        x, y = stats.numeric(rows, xk), stats.numeric(rows, yk)
        cells = sorted({r['cell'] for r in rows})
        for n, c in enumerate(cells):
            m = np.array([r['cell'] == c for r in rows])
            ax.scatter(x[m], y[m], s=40 if self.level.currentData() == 'cell' else 12, alpha=0.75,
                       color=CELL_COLORS[n % len(CELL_COLORS)], label=c, edgecolors='none')
            if self.level.currentData() == 'cell':
                for xi, yi in zip(x[m], y[m]):
                    ax.annotate(c[4:], (xi, yi), textcoords='offset points', xytext=(4, 4), fontsize=8)
        if self.logx.isChecked():
            ax.set_xscale('symlog' if np.nanmin(x) <= 0 else 'log')
        if self.logy.isChecked():
            ax.set_yscale('symlog' if np.nanmin(y) <= 0 else 'log')
        ax.set_xlabel(self.x.currentText()); ax.set_ylabel(self.y.currentText())
        ax.legend(fontsize=7, frameon=False, ncol=2)
        ax.spines[['top', 'right']].set_visible(False)
        rho, p, n = stats.spearman(x, y)
        ax.set_title(f'ρ = {rho:.3f}, p = {p:.2g}, n = {n}', fontsize=10)
        f.tight_layout(); self.scatter.draw_idle()
        self.stat.setText(f'{self.y.currentText()} vs {self.x.currentText()}: Spearman ρ = {rho:.3f}, '
                          f'p = {p:.3g}, n = {n}')


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
        self.incl_edge = QCheckBox('Also analyse cells touching the image border')
        self.min_area = self._spin(0, 1e5, 250, 1, ' µm²')
        self.binuc_tau = self._spin(-1, 1, 0.0, 3, '', 0.01,
                                    'Relative border contrast below which two neighbouring cells count as a weak-border pair')
        self.px = self._spin(0, 10, 0.0, 4, ' µm/px', 0.001, '0 = read the pixel size from the TIFF')
        self.sens = self._spin(0.2, 5, 1.0, 2, ' ×', 0.05,
                               'Multiplies the automatic green puncta threshold: >1 keeps fewer, brighter puncta')
        self.min_mito = self._spin(0, 10, 0.05, 3, ' µm²', 0.01, 'Mito pieces smaller than this are left out of the per-mito table')
        of.addRow('Mito segmentation', self.mito_method)
        of.addRow(self.excl_binuc); of.addRow(self.incl_edge)
        of.addRow('Minimum cell area', self.min_area)
        of.addRow('Weak-border threshold', self.binuc_tau)
        of.addRow('Pixel size (0 = auto)', self.px)
        of.addRow('Green puncta threshold', self.sens)
        of.addRow('Minimum mito object', self.min_mito)

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
        self.views = {k: ImageView() for k in ('green_overlay', 'overlay', 'cells', 'summary')}
        self.tabs.addTab(self.cells_tab, 'Cells')
        self.tabs.addTab(self.mito_tab, 'Mito objects')
        self.tabs.addTab(self.puncta_tab, 'Green puncta')
        self.tabs.addTab(self.corr_tab, 'Correlation')
        self.tabs.addTab(self.views['green_overlay'], 'Green on mito')
        self.tabs.addTab(self.views['overlay'], 'MiNA overlay')
        self.tabs.addTab(self.views['cells'], 'Cells (zoom)')
        self.tabs.addTab(self.views['summary'], 'Cell ROIs')

        split = QSplitter(); split.addWidget(left); split.addWidget(self.tabs)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1); split.setSizes([620, 880])
        self.setCentralWidget(split)
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
                                 include_edge_cells=self.incl_edge.isChecked(), pixel_size_um=self.px.value(),
                                 puncta_sensitivity=self.sens.value(), min_mito_area_um2=self.min_mito.value(),
                                 mito_method=self.mito_method.currentData())
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
    return [dict(metric=k, **{r['cell']: r[k] for r in rows}) for k in rows[0] if k != 'cell']


def main():
    app = QApplication(sys.argv)
    app.setApplicationName('Mito Analyzer')
    if os.path.exists(resource('icon.png')):
        app.setWindowIcon(QIcon(resource('icon.png')))
    w = MainWindow(); w.show()
    return app.exec()
