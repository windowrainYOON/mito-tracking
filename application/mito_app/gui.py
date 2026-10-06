"""PySide6 desktop GUI for the per-cell MiNA + green-on-mito pipeline."""
import math, os, sys, traceback

os.environ.setdefault('QT_API', 'pyside6')

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QObject, QSortFilterProxyModel, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QIcon, QPixmap
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QHeaderView, QLabel, QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea, QSplitter,
    QTableView, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__, pipeline, stats

TIFF_FILTER = 'TIFF images (*.tif *.tiff);;All files (*)'
CHANNELS = (('red', 'Red (mitochondria)'), ('green', 'Green (protein of interest)'), ('blue', 'Blue (nuclei)'))
CELL_COLORS = ['#4C9BE8', '#F28E2B', '#59A14F', '#E15759', '#B07AA1', '#EDC948', '#76B7B2', '#FF9DA7',
               '#9C755F', '#BAB0AC']


def resource(name):
    base = getattr(sys, '_MEIPASS', os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    return os.path.join(base, 'resources', name)


class Worker(QObject):
    log = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, kwargs):
        super().__init__()
        self.kwargs = kwargs

    def run(self):
        try:
            self.finished.emit(pipeline.run(log=self.log.emit, **self.kwargs))
        except Exception as e:  # report to the user instead of crashing the app
            self.log.emit(traceback.format_exc())
            self.failed.emit(f'{type(e).__name__}: {e}')


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

        # inputs
        self.paths = {}
        inp = QGroupBox('Input images (same field)')
        form = QFormLayout(inp)
        for key, label in CHANNELS:
            edit = QLineEdit(); edit.setPlaceholderText('Choose a TIFF file…')
            btn = QPushButton('Browse…'); btn.clicked.connect(lambda _=False, k=key: self.browse(k))
            row = QHBoxLayout(); row.addWidget(edit, 1); row.addWidget(btn)
            form.addRow(label, row)
            self.paths[key] = edit
        self.outdir = QLineEdit(); self.outdir.setPlaceholderText('Defaults to a "results" folder next to the red image')
        ob = QPushButton('Browse…'); ob.clicked.connect(self.browse_out)
        row = QHBoxLayout(); row.addWidget(self.outdir, 1); row.addWidget(ob)
        form.addRow('Output folder', row)
        self.name = QLineEdit(); self.name.setPlaceholderText('Taken from the red file name')
        form.addRow('Sample name', self.name)

        # options
        opt = QGroupBox('Options')
        of = QFormLayout(opt)
        self.mito_method = QComboBox()
        self.mito_method.addItem('Split objects: adaptive threshold + watershed', 'split')
        self.mito_method.addItem('MiNA classic: Otsu per cell', 'otsu')
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

        self.run_btn = QPushButton('Run analysis'); self.run_btn.setDefault(True)
        self.run_btn.setMinimumHeight(36); self.run_btn.clicked.connect(self.start)
        self.open_btn = QPushButton('Open output folder'); self.open_btn.setEnabled(False)
        self.open_btn.clicked.connect(self.open_out)
        btns = QHBoxLayout(); btns.addWidget(self.run_btn, 2); btns.addWidget(self.open_btn, 1)

        self.logbox = QPlainTextEdit(); self.logbox.setReadOnly(True)
        self.logbox.setPlaceholderText('Progress appears here.')

        left = QWidget(); lv = QVBoxLayout(left)
        lv.addWidget(inp); lv.addWidget(opt); lv.addLayout(btns); lv.addWidget(QLabel('Log')); lv.addWidget(self.logbox, 1)

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
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1); split.setSizes([440, 960])
        self.setCentralWidget(split)
        self.statusBar().showMessage('Choose the red, green and blue images, then press Run.')

    @staticmethod
    def _spin(lo, hi, val, dec, suffix='', step=None, tip=None):
        s = QDoubleSpinBox(); s.setRange(lo, hi); s.setDecimals(dec); s.setValue(val); s.setSuffix(suffix)
        if step:
            s.setSingleStep(step)
        if tip:
            s.setToolTip(tip)
        return s

    # file pickers
    def browse(self, key):
        start = self.paths[key].text() or next((e.text() for e in self.paths.values() if e.text()), '')
        path, _ = QFileDialog.getOpenFileName(self, f'Choose {key} image', os.path.dirname(start), TIFF_FILTER)
        if not path:
            return
        self.paths[key].setText(path)
        for k, p in pipeline.find_siblings(path).items():  # fill the other channels if their names match
            if not self.paths[k].text():
                self.paths[k].setText(p)

    def browse_out(self):
        d = QFileDialog.getExistingDirectory(self, 'Choose output folder', self.outdir.text())
        if d:
            self.outdir.setText(d)

    def open_out(self):
        if self.result:
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.result['outdir']))

    # run
    def start(self):
        p = {k: e.text().strip() for k, e in self.paths.items()}
        missing = [label for key, label in CHANNELS if not os.path.isfile(p[key])]
        if missing:
            QMessageBox.warning(self, 'Missing input', 'Choose an existing file for:\n' + '\n'.join(missing))
            return
        name = self.name.text().strip() or pipeline.sample_name(p['red'])
        outdir = self.outdir.text().strip() or os.path.join(os.path.dirname(p['red']), 'results', name)
        params = pipeline.Params(min_area_um2=self.min_area.value(), binuc_tau=self.binuc_tau.value(),
                                 exclude_binucleate=self.excl_binuc.isChecked(),
                                 include_edge_cells=self.incl_edge.isChecked(), pixel_size_um=self.px.value(),
                                 puncta_sensitivity=self.sens.value(), min_mito_area_um2=self.min_mito.value(),
                                 mito_method=self.mito_method.currentData())
        self.logbox.clear()
        self.run_btn.setEnabled(False); self.open_btn.setEnabled(False)
        self.statusBar().showMessage('Running…')
        self.thread = QThread()
        self.worker = Worker(dict(red_path=p['red'], green_path=p['green'], blue_path=p['blue'],
                                  outdir=outdir, params=params, name=name))
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.log.connect(self.logbox.appendPlainText)
        self.worker.finished.connect(self.done)
        self.worker.failed.connect(self.fail)
        for sig in (self.worker.finished, self.worker.failed):
            sig.connect(self.thread.quit)
        self.thread.start()

    def done(self, res):
        self.result = res
        self.run_btn.setEnabled(True); self.open_btn.setEnabled(True)
        self.statusBar().showMessage(f"Analysed {len(res['rows'])} cells, {len(res['mito_rows'])} mito objects, "
                                     f"{len(res['puncta_rows'])} green puncta. Results saved to {res['outdir']}")
        self.cells_tab.set_rows(transpose(res['rows']))
        self.mito_tab.set_rows(res['mito_rows'])
        self.puncta_tab.set_rows(res['puncta_rows'])
        self.corr_tab.set_data(res['rows'], res['mito_rows'])
        for k, v in self.views.items():
            v.set_image(res[k])
        self.tabs.setCurrentIndex(0)

    def fail(self, msg):
        self.run_btn.setEnabled(True)
        self.statusBar().showMessage('Analysis failed')
        QMessageBox.critical(self, 'Analysis failed', msg + '\n\nSee the log for details.')

    def closeEvent(self, e):
        if self.thread and self.thread.isRunning():
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
