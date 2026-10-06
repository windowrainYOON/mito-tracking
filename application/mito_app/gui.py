"""PySide6 desktop GUI for the per-cell MiNA pipeline."""
import csv, os, sys, traceback

from PySide6.QtCore import QObject, Qt, QThread, QUrl, Signal
from PySide6.QtGui import QDesktopServices, QPixmap
from PySide6.QtWidgets import (
    QApplication, QCheckBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout, QLabel,
    QLineEdit, QMainWindow, QMessageBox, QPlainTextEdit, QPushButton, QScrollArea, QSplitter, QTableWidget,
    QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget,
)

from . import __version__, pipeline

TIFF_FILTER = 'TIFF images (*.tif *.tiff);;All files (*)'
CHANNELS = (('red', 'Red (mitochondria)'), ('green', 'Green (protein of interest)'), ('blue', 'Blue (nuclei)'))


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


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f'Mito Analyzer {__version__}')
        self.resize(1280, 860)
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
        self.excl_binuc = QCheckBox('Exclude cells that share a weak border (look binucleate)')
        self.incl_edge = QCheckBox('Also analyse cells touching the image border')
        self.min_area = QDoubleSpinBox(); self.min_area.setRange(0, 1e5); self.min_area.setValue(250)
        self.min_area.setSuffix(' µm²')
        self.binuc_tau = QDoubleSpinBox(); self.binuc_tau.setRange(-1, 1); self.binuc_tau.setDecimals(3)
        self.binuc_tau.setSingleStep(0.01); self.binuc_tau.setValue(0.0)
        self.binuc_tau.setToolTip('Relative border contrast below which two neighbouring cells count as a weak-border pair')
        of.addRow(self.excl_binuc); of.addRow(self.incl_edge)
        of.addRow('Minimum cell area', self.min_area)
        of.addRow('Weak-border threshold', self.binuc_tau)

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
        self.table = QTableWidget()
        self.views = {k: ImageView() for k in ('overlay', 'cells', 'summary')}
        self.tabs.addTab(self.table, 'Per-cell results')
        self.tabs.addTab(self.views['overlay'], 'MiNA overlay')
        self.tabs.addTab(self.views['cells'], 'Cells (zoom)')
        self.tabs.addTab(self.views['summary'], 'Cell ROIs')

        split = QSplitter(); split.addWidget(left); split.addWidget(self.tabs)
        split.setStretchFactor(0, 0); split.setStretchFactor(1, 1); split.setSizes([440, 840])
        self.setCentralWidget(split)
        self.statusBar().showMessage('Choose the red, green and blue images, then press Run.')

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
                                 include_edge_cells=self.incl_edge.isChecked())
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
        self.statusBar().showMessage(f"Analysed {len(res['rows'])} cells. Results saved to {res['outdir']}")
        self.load_table(res['per_cell_csv'])
        for k, v in self.views.items():
            v.set_image(res[k])
        self.tabs.setCurrentIndex(0)

    def fail(self, msg):
        self.run_btn.setEnabled(True)
        self.statusBar().showMessage('Analysis failed')
        QMessageBox.critical(self, 'Analysis failed', msg + '\n\nSee the log for details.')

    def load_table(self, path):
        with open(path, newline='') as f:
            data = list(csv.reader(f))
        if not data:
            return
        head, body = data[0], data[1:]
        # cells as columns, metrics as rows: easier to read with ~20 metrics and a handful of cells
        self.table.clear()
        self.table.setRowCount(len(head) - 1); self.table.setColumnCount(len(body))
        self.table.setHorizontalHeaderLabels([r[0] for r in body])
        self.table.setVerticalHeaderLabels(head[1:])
        for j, r in enumerate(body):
            for i, v in enumerate(r[1:]):
                it = QTableWidgetItem(v); it.setFlags(it.flags() & ~Qt.ItemIsEditable)
                it.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.table.setItem(i, j, it)
        self.table.resizeColumnsToContents()

    def closeEvent(self, e):
        if self.thread and self.thread.isRunning():
            if QMessageBox.question(self, 'Quit', 'An analysis is still running. Quit anyway?') != QMessageBox.Yes:
                e.ignore(); return
        e.accept()


def main():
    app = QApplication(sys.argv)
    app.setApplicationName('Mito Analyzer')
    w = MainWindow(); w.show()
    return app.exec()
