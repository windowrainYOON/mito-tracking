"""Results browser: look at, filter, export and rename results that were already analysed.

Separate from the analysis (it never runs the pipeline): it reads the <sample>_per_cell.csv / _per_mito.csv /
_green_puncta.csv files under a results folder (pipeline.load_results) and shows one scope at a time (cells,
mito objects, green puncta, one row per image, one row per group) with only the parameters picked in the list.
The table as shown (picked parameters, filtered rows) is exported to Excel.

Renaming never touches the result files: group names go to groups.csv, condition / sample / cell names to
names.csv in the results folder, and every view that loads the results (this browser, the Analysis tab, the group
export) shows the new names. Clearing a name (or Reset names) goes back to the original one."""
import datetime, math, os

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSettings, QSortFilterProxyModel, Qt, QUrl
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QAbstractItemView, QApplication, QCheckBox, QComboBox, QFileDialog, QFormLayout, QGroupBox, QHBoxLayout,
    QHeaderView, QLabel, QLineEdit, QListWidget, QListWidgetItem, QMainWindow, QMessageBox, QPushButton, QSplitter,
    QTableView, QVBoxLayout, QWidget,
)

from . import pipeline, stats

SCOPES = (('cells', 'Cells (one row per cell)'), ('mito', 'Mito objects (one row per object)'),
          ('puncta', 'Green puncta (one row per punctum)'), ('images', 'Images (mean per image)'),
          ('groups', 'Groups (summary per group)'))
SOURCES = (('cells', 'cells'), ('mito', 'mito objects'), ('puncta', 'green puncta'))
# identifying columns per source, always shown (dataset is shown as "Condition")
ID_COLS = {'cells': ('group', 'dataset', 'preset', 'sample', 'cell', 'green_status', 'status'),
           'mito': ('group', 'dataset', 'preset', 'sample', 'cell', 'mito', 'green_status'),
           'puncta': ('group', 'dataset', 'preset', 'sample', 'cell', 'punctum', 'green_status')}
ID_LABELS = {'group': 'Group', 'dataset': 'Condition', 'preset': 'Preset', 'sample': 'Sample', 'cell': 'Cell',
             'mito': 'Mito object', 'punctum': 'Punctum', 'green_status': 'Green', 'status': 'QC',
             'n_images': 'Images', 'n_rows': 'n'}
# columns that are positions or links, not parameters
NOT_PARAMS = {'centroid_x', 'centroid_y', 'bbox_x', 'bbox_y', 'bbox_w', 'bbox_h', 'touches_border',
              'weak_border_with', 'seed', 'on_mito', 'mito', 'cell', 'punctum', 'n_rows'}
PUNCTA_LABELS = {'area_um2': 'Area (µm²)', 'green_mean': 'Green mean', 'green_max': 'Green max',
                 'green_integrated': 'Green integrated', 'overlap_fraction': 'Overlap with mito (fraction)',
                 'distance_to_mito_um': 'Distance to mito (µm)'}
STATS = (('mean', 'mean'), ('sd', 'SD'), ('sem', 'SEM'), ('median', 'median'), ('n', 'n'))
EDITABLE = ('group', 'dataset', 'sample', 'cell')


def labels_for(source):
    pairs = {'cells': stats.CELL_X + stats.CELL_Y, 'mito': stats.MITO_X + stats.MITO_Y}.get(source, [])
    out = dict(PUNCTA_LABELS) if source == 'puncta' else {}
    out.update(dict(pairs))
    return out


def is_number(v):
    return isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, (bool, np.bool_))


def param_keys(rows, source):
    """Numeric parameter columns of the rows, in file order."""
    skip = set(ID_COLS[source]) | NOT_PARAMS
    keys = []
    for r in rows:
        for k, v in r.items():
            if k not in skip and not k.startswith('_') and k not in keys and is_number(v):
                keys.append(k)
    order = {k: i for i, k in enumerate(rows[0])} if rows else {}
    return sorted(keys, key=lambda k: order.get(k, 1e9))


def load_puncta(samples, cells):
    """Rows of every <sample>_green_puncta.csv, with the same group / condition / sample / cell names as `cells`."""
    names = {(r['_folder'], r['_cell']): r['cell'] for r in cells}
    out = []
    for s in samples:
        p = os.path.join(s['folder'], f"{s['stem']}_green_puncta.csv")
        if not os.path.exists(p):
            continue
        meta = {k: s[k] for k in ('group', 'dataset', 'preset', 'sample')}
        nd = os.path.normpath(s['folder'])
        for r in pipeline.read_csv(p):
            c = str(r.get('cell', ''))
            r = dict(meta, **r)
            r['cell'] = names.get((nd, c), c)
            r['_folder'], r['_cell'] = nd, c
            out.append(r)
    return out


class Table(QAbstractTableModel):
    """Rows of dicts shown as `cols` with `headers`; `editable` columns call on_edit(row, key, text)."""

    def __init__(self):
        super().__init__()
        self.rows, self.cols, self.headers, self.editable, self.on_edit = [], [], [], set(), None

    def set_rows(self, rows, cols, headers, editable=()):
        self.beginResetModel()
        self.rows, self.cols, self.headers, self.editable = rows, list(cols), list(headers), set(editable)
        self.endResetModel()

    def rowCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()):
        return 0 if parent.isValid() else len(self.cols)

    def flags(self, idx):
        f = super().flags(idx)
        return f | Qt.ItemIsEditable if self.cols[idx.column()] in self.editable else f

    def data(self, idx, role=Qt.DisplayRole):
        v = self.rows[idx.row()].get(self.cols[idx.column()], '')
        if role in (Qt.DisplayRole, Qt.EditRole):
            if isinstance(v, float):
                return '' if math.isnan(v) else f'{v:.4g}'
            return str(v)
        if role == Qt.UserRole:
            if is_number(v):
                return -math.inf if isinstance(v, float) and math.isnan(v) else float(v)
            return str(v)
        if role == Qt.TextAlignmentRole and is_number(v):
            return int(Qt.AlignRight | Qt.AlignVCenter)
        if role == Qt.ToolTipRole and self.cols[idx.column()] in self.editable:
            return 'Double-click to rename (the result files are not changed; clear the name to go back)'

    def setData(self, idx, value, role=Qt.EditRole):
        if role != Qt.EditRole or self.on_edit is None:
            return False
        r, k = self.rows[idx.row()], self.cols[idx.column()]
        if str(value).strip() == str(r.get(k, '')):
            return False
        self.on_edit(r, k, str(value).strip())
        return True

    def headerData(self, s, orient, role=Qt.DisplayRole):
        if role == Qt.DisplayRole:
            return self.headers[s] if orient == Qt.Horizontal else str(s + 1)
        if role == Qt.ToolTipRole and orient == Qt.Horizontal:
            return self.cols[s]


class CheckList(QWidget):
    """A list of checkable items with All / None buttons and an optional text filter."""

    def __init__(self, title, search=False):
        super().__init__()
        self.list = QListWidget()
        self.search = QLineEdit(); self.search.setPlaceholderText('Search…'); self.search.setClearButtonEnabled(True)
        all_b, none_b = QPushButton('All'), QPushButton('None')
        for b in (all_b, none_b):
            b.setMaximumWidth(60)
        all_b.clicked.connect(lambda: self.set_all(True)); none_b.clicked.connect(lambda: self.set_all(False))
        head = QHBoxLayout(); head.addWidget(QLabel(f'<b>{title}</b>')); head.addStretch(1)
        head.addWidget(all_b); head.addWidget(none_b)
        lay = QVBoxLayout(self); lay.setContentsMargins(0, 0, 0, 0); lay.addLayout(head)
        if search:
            lay.addWidget(self.search)
            self.search.textChanged.connect(self.filter)
        else:
            self.search.hide()
        lay.addWidget(self.list, 1)

    def set_items(self, items, checked=None):
        """items: (key, label); checked: keys to tick (None = keep the current ticks, new items ticked)."""
        old = {self.list.item(i).data(Qt.UserRole): self.list.item(i).checkState() == Qt.Checked
               for i in range(self.list.count())}
        self.list.blockSignals(True); self.list.clear()
        for key, label in items:
            it = QListWidgetItem(label); it.setData(Qt.UserRole, key)
            it.setFlags(it.flags() | Qt.ItemIsUserCheckable)
            on = (key in checked) if checked is not None else old.get(key, True)
            it.setCheckState(Qt.Checked if on else Qt.Unchecked)
            self.list.addItem(it)
        self.list.blockSignals(False)
        self.filter(self.search.text())

    def checked(self):
        return [self.list.item(i).data(Qt.UserRole) for i in range(self.list.count())
                if self.list.item(i).checkState() == Qt.Checked]

    def set_all(self, on):
        self.list.blockSignals(True)
        for i in range(self.list.count()):
            if not self.list.item(i).isHidden():
                self.list.item(i).setCheckState(Qt.Checked if on else Qt.Unchecked)
        self.list.blockSignals(False)
        self.list.itemChanged.emit(self.list.item(0)) if self.list.count() else None

    def filter(self, text):
        t = text.strip().lower()
        for i in range(self.list.count()):
            it = self.list.item(i)
            it.setHidden(bool(t) and t not in f'{it.text()} {it.data(Qt.UserRole)}'.lower())


class ResultsBrowser(QMainWindow):
    def __init__(self, folder='', parent=None):
        super().__init__(parent)
        self.setWindowTitle('Mito Analyzer: results browser'); self.resize(1400, 860)
        self.root, self.data, self.samples = '', {}, []
        self.settings = QSettings('MitoAnalyzer', 'ResultsBrowser')
        self.picked = {}   # source -> parameter keys ticked

        self.folder = QLineEdit(); self.folder.setPlaceholderText('Results folder of earlier runs (searched recursively)')
        self.folder.returnPressed.connect(self.load)
        browse = QPushButton('Open folder…'); browse.clicked.connect(self.browse)
        reload_b = QPushButton('Reload'); reload_b.clicked.connect(self.load)
        self.reset_b = QPushButton('Reset names…'); self.reset_b.clicked.connect(self.reset_names)
        self.reset_b.setToolTip('Go back to the original condition, sample and cell names (group names are kept)')
        top = QHBoxLayout(); top.addWidget(QLabel('Results folder')); top.addWidget(self.folder, 1)
        for b in (browse, reload_b, self.reset_b):
            top.addWidget(b)

        # left: scope, filters, parameters
        self.scope = QComboBox()
        for k, l in SCOPES:
            self.scope.addItem(l, k)
        self.source = QComboBox()
        for k, l in SOURCES:
            self.source.addItem(f'Mean of {l}', k)
        self.unit = QComboBox(); self.unit.addItem('each cell / object is one value', 'rows')
        self.unit.addItem('each image mean is one value', 'images')
        self.stat_boxes = {}
        stat_row = QHBoxLayout()
        for k, l in STATS:
            cb = QCheckBox(l); cb.setChecked(k in ('mean', 'sd', 'n')); cb.toggled.connect(self.refresh)
            self.stat_boxes[k] = cb; stat_row.addWidget(cb)
        self.green = QComboBox()
        for l, k in (('All cells', ''), ('Green+ only', 'positive'), ('Green− only', 'negative')):
            self.green.addItem(l, k)
        self.qc = QCheckBox('Only cells that pass QC'); self.qc.setChecked(True)
        self.qc.setToolTip('Leave out cells touching the image border (and weak-border pairs, if excluded) — '
                           'the cells that were not analysed')
        form = QFormLayout()
        form.addRow('Scope', self.scope); form.addRow('Values', self.source); form.addRow('Group statistic', self.unit)
        form.addRow('', stat_row); form.addRow('Cells', self.green); form.addRow('', self.qc)
        self.f_group = CheckList('Groups'); self.f_cond = CheckList('Conditions'); self.f_sample = CheckList('Samples', True)
        self.params = CheckList('Parameters', True)
        filt = QSplitter(Qt.Vertical)
        for w in (self.f_group, self.f_cond, self.f_sample):
            filt.addWidget(w)
        left = QSplitter(Qt.Vertical)
        fbox = QWidget(); fl = QVBoxLayout(fbox); fl.setContentsMargins(0, 0, 0, 0); fl.addLayout(form); fl.addWidget(filt, 1)
        left.addWidget(fbox); left.addWidget(self.params); left.setSizes([480, 380])

        # right: table + export
        self.model = Table(); self.model.on_edit = self.rename
        self.proxy = QSortFilterProxyModel(); self.proxy.setSourceModel(self.model); self.proxy.setSortRole(Qt.UserRole)
        self.view = QTableView(); self.view.setModel(self.proxy); self.view.setSortingEnabled(True)
        self.view.setAlternatingRowColors(True)
        self.view.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed)
        self.view.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.count = QLabel()
        exp1 = QPushButton('Export this table to Excel…'); exp1.clicked.connect(self.export_table)
        exp2 = QPushButton('Export every scope to Excel…'); exp2.clicked.connect(self.export_all)
        exp2.setToolTip('One sheet per scope (cells, mito objects, puncta, images, groups) with the same filters '
                        'and the parameters ticked for each kind of value')
        hint = QLabel('Double-click a group, condition, sample or cell name to rename it.')
        hint.setStyleSheet('color: gray')
        bottom = QHBoxLayout(); bottom.addWidget(self.count); bottom.addWidget(hint); bottom.addStretch(1)
        bottom.addWidget(exp1); bottom.addWidget(exp2)
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(self.view, 1); rl.addLayout(bottom)

        split = QSplitter(); split.addWidget(left); split.addWidget(right); split.setSizes([360, 1040])
        self.info = QLabel(); self.info.setWordWrap(True)
        central = QWidget(); lay = QVBoxLayout(central); lay.addLayout(top); lay.addWidget(self.info); lay.addWidget(split, 1)
        self.setCentralWidget(central)

        self.scope.currentIndexChanged.connect(self.scope_changed)
        self.source.currentIndexChanged.connect(self.scope_changed)
        for w in (self.unit, self.green):
            w.currentIndexChanged.connect(self.refresh)
        self.qc.toggled.connect(self.refresh)
        for c in (self.f_group, self.f_cond, self.f_sample):
            c.list.itemChanged.connect(self.refresh)
        self.params.list.itemChanged.connect(self.params_changed)
        self.set_enabled(False)
        if folder:
            self.folder.setText(folder); self.load()

    # ---- loading ------------------------------------------------------------

    def set_enabled(self, on):
        for w in (self.scope, self.source, self.unit, self.green, self.qc, self.reset_b):
            w.setEnabled(on)

    def browse(self):
        d = QFileDialog.getExistingDirectory(self, 'Results folder', self.folder.text() or os.path.expanduser('~'))
        if d:
            self.folder.setText(d); self.load()

    def load(self, *_, keep=False):
        d = os.path.expanduser(self.folder.text().strip())
        if not os.path.isdir(d):
            QMessageBox.warning(self, 'Not a folder', f'{d or "(empty)"} is not a folder.'); return
        QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            cells, mito, samples = pipeline.load_results(d, keys=True)
            puncta = load_puncta(samples, cells)
        finally:
            QApplication.restoreOverrideCursor()
        if not samples:
            QMessageBox.warning(self, 'Nothing found', f'No <sample>_per_cell.csv found under {d}.'); return
        self.root, self.samples = d, samples
        self.data = {'cells': cells, 'mito': mito, 'puncta': puncta}
        status = {(r['_folder'], r['_cell']): r.get('status', 'ok') for r in cells}
        for r in mito + puncta:
            r['_status'] = status.get((r['_folder'], r['_cell']), 'ok')
        for r in cells:
            r['_status'] = r.get('status', 'ok')
        groups = sorted({s['group'] for s in samples}); conds = sorted({s['dataset'] for s in samples})
        self.f_group.set_items([(g, g) for g in groups], None if keep else set(groups))
        self.f_cond.set_items([(c, c) for c in conds], None if keep else set(conds))
        self.f_sample.set_items([(os.path.normpath(s['folder']), f"{s['dataset']} / {s['sample']}") for s in samples],
                                None if keep else {os.path.normpath(s['folder']) for s in samples})
        self.info.setText(f'{len(samples)} images in {len(groups)} group(s), {len(conds)} condition(s): '
                          f'{len(cells)} cells, {len(mito)} mito objects, {len(puncta)} green puncta.  {d}')
        self.set_enabled(True)
        self.scope_changed()

    # ---- what is shown ------------------------------------------------------

    def current_source(self):
        s = self.scope.currentData()
        return s if s in ID_COLS else self.source.currentData()

    def scope_changed(self, *_):
        sc = self.scope.currentData()
        self.source.setEnabled(sc in ('images', 'groups')); self.unit.setEnabled(sc == 'groups')
        for cb in self.stat_boxes.values():
            cb.setEnabled(sc == 'groups')
        src = self.current_source()
        lab = labels_for(src)
        keys = param_keys(self.data.get(src, []), src)
        self.params.set_items([(k, lab.get(k, k)) for k in keys], self.picked_keys(src))
        self.picked[src] = set(self.params.checked())
        self.refresh()

    def picked_keys(self, src):
        """Parameters ticked for a kind of value: this session's choice, else the last one saved, else all."""
        if src not in self.picked:
            try:
                stored = self.settings.value(f'params/{src}')
            except Exception:
                stored = None
            stored = [stored] if isinstance(stored, str) and stored else stored
            self.picked[src] = set(stored) if isinstance(stored, list) else \
                set(param_keys(self.data.get(src, []), src))
        return self.picked[src]

    def params_changed(self, *_):
        src = self.current_source()
        self.picked[src] = set(self.params.checked())
        try:
            self.settings.setValue(f'params/{src}', sorted(self.picked[src]))
        except Exception:
            pass
        self.refresh()

    def filtered(self, src):
        groups, conds = set(self.f_group.checked()), set(self.f_cond.checked())
        folders, g = set(self.f_sample.checked()), self.green.currentData()
        return [r for r in self.data.get(src, [])
                if r['group'] in groups and r['dataset'] in conds and r['_folder'] in folders
                and (not g or r.get('green_status') == g) and (not self.qc.isChecked() or r['_status'] == 'ok')]

    def table(self, scope, src=None, keys=None):
        """(rows, cols, headers, editable) of one scope with the current filters."""
        src = src or (scope if scope in ID_COLS else self.source.currentData())
        rows = self.filtered(src)
        lab = labels_for(src)
        if keys is None:
            keys = [k for k in param_keys(self.data.get(src, []), src) if k in self.picked_keys(src)]
        multi_preset = len({s['preset'] for s in self.samples}) > 1
        if scope in ID_COLS:
            ids = [c for c in ID_COLS[scope] if c != 'preset' or multi_preset]
            ids = [c for c in ids if any(c in r for r in rows[:1])] if rows else ids
            cols = ids + keys
            return rows, cols, [ID_LABELS.get(c, c) for c in ids] + [lab.get(k, k) for k in keys], \
                [c for c in EDITABLE if c in ids]
        by = ('group', 'dataset', 'preset', 'sample', '_folder')
        img = stats.sample_means(rows, keys, by=by)
        if scope == 'images':
            ids = ['group', 'dataset'] + (['preset'] if multi_preset else []) + ['sample', 'n_rows']
            heads = [ID_LABELS[c] for c in ids[:-1]] + [f'n {dict(SOURCES)[src]}']
            return img, ids + keys, heads + [lab.get(k, k) for k in keys], ['group', 'dataset', 'sample']
        # groups
        vals = img if self.unit.currentData() == 'images' else rows
        picked = [k for k, _ in STATS if self.stat_boxes[k].isChecked()] or ['mean']
        out = []
        for gname in sorted({r['group'] for r in rows}):
            gr = [r for r in vals if r['group'] == gname]
            d = dict(group=gname, n_images=len({r['_folder'] for r in rows if r['group'] == gname}),
                     n_rows=sum(r['group'] == gname for r in rows))
            for k in keys:
                desc = stats.describe(stats.numeric(gr, k))
                for s in picked:
                    d[f'{k}|{s}'] = desc[s]
            out.append(d)
        cols = ['group', 'n_images', 'n_rows'] + [f'{k}|{s}' for k in keys for s in picked]
        heads = ['Group', 'Images', f'n {dict(SOURCES)[src]}'] + [f'{lab.get(k, k)} {dict(STATS)[s]}'
                                                                   for k in keys for s in picked]
        return out, cols, heads, ['group']

    def refresh(self, *_):
        if not self.data:
            return
        rows, cols, heads, editable = self.table(self.scope.currentData())
        self.model.set_rows(rows, cols, heads, editable)
        self.view.sortByColumn(-1, Qt.AscendingOrder); self.proxy.sort(-1)
        self.view.resizeColumnsToContents()
        for i in range(len(cols)):
            if self.view.columnWidth(i) > 220:
                self.view.setColumnWidth(i, 220)
        n_par = len(self.picked.get(self.current_source(), ()))
        self.count.setText(f'{len(rows)} rows, {n_par} parameter(s)')

    # ---- renaming -----------------------------------------------------------

    def rename(self, row, key, new):
        """Rename a group / condition / sample / cell (empty = back to the original name)."""
        if key == 'group' and 'n_images' in row:   # groups scope: the whole group
            folders = [s['folder'] for s in self.samples if s['group'] == row['group']]
        elif key in ('group', 'dataset'):
            same = [s['folder'] for s in self.samples if s[key] == row[key]]
            folders = [row['_folder']]
            if len(same) > 1:
                what = ID_LABELS[key].lower()
                box = QMessageBox(QMessageBox.Question, f'Rename {what}',
                                  f'{len(same)} images have the {what} "{row[key]}". Rename it for all of them, '
                                  f'or move only "{row["sample"]}" to "{new or "(original)"}"?', parent=self)
                all_b = box.addButton(f'All {len(same)} images', QMessageBox.AcceptRole)
                one_b = box.addButton('Only this image', QMessageBox.ActionRole)
                box.addButton(QMessageBox.Cancel)
                box.exec()
                if box.clickedButton() == all_b:
                    folders = same
                elif box.clickedButton() != one_b:
                    return
        else:
            folders = [row['_folder']]
        by_folder = {os.path.normpath(s['folder']): s for s in self.samples}
        try:
            if key == 'group':
                pipeline.write_groups(self.root, [dict(outdir=f, dataset=by_folder[os.path.normpath(f)]['dataset'],
                                                       sample=by_folder[os.path.normpath(f)]['sample'], group=new)
                                                  for f in folders])
            elif key == 'cell':
                pipeline.write_names(self.root, [dict(folder=row['_folder'], cell=row['_cell'], field='cell',
                                                      name='' if new == row['_cell'] else new)])
            else:
                field = 'condition' if key == 'dataset' else 'sample'
                if field == 'sample':
                    s0 = by_folder[os.path.normpath(row['_folder'])]
                    clash = [s for s in self.samples if s is not s0 and s['dataset'] == s0['dataset']
                             and s['group'] == s0['group'] and s['sample'] == new]
                    if clash:
                        QMessageBox.warning(self, 'Rename sample', f'Another image in this group and condition is '
                                            f'already called "{new}"; their cells would be pooled as one image.')
                        return
                original = lambda f: (os.path.basename(os.path.dirname(os.path.dirname(f))) if field == 'condition'
                                      else by_folder[os.path.normpath(f)]['stem'])
                pipeline.write_names(self.root, [dict(folder=f, cell='', field=field,
                                                      name='' if new == original(f) else new) for f in folders])
        except OSError as e:
            QMessageBox.warning(self, 'Rename', f'Cannot write to {self.root}:\n{e}'); return
        self.load(keep=True)

    def reset_names(self):
        if not self.root:
            return
        if QMessageBox.question(self, 'Reset names', 'Go back to the original condition, sample and cell names '
                                'for every image under this folder? (Groups are kept; change them by renaming.)'
                                ) != QMessageBox.Yes:
            return
        try:
            pipeline.clear_names(self.root)
        except OSError as e:
            QMessageBox.warning(self, 'Reset names', str(e)); return
        self.load(keep=True)

    # ---- export -------------------------------------------------------------

    @staticmethod
    def sheet(rows, cols, heads):
        out = []
        for r in rows:
            out.append({h: r.get(c) for c, h in zip(cols, heads)})
        return out or [dict.fromkeys(heads)]

    def info_sheet(self):
        return [dict(item=k, value=v) for k, v in (
            ('results folder', self.root), ('exported', datetime.datetime.now().strftime('%Y-%m-%d %H:%M')),
            ('groups', ', '.join(self.f_group.checked())), ('conditions', ', '.join(self.f_cond.checked())),
            ('images', f'{len(self.f_sample.checked())} of {len(self.samples)}'),
            ('cells', self.green.currentText() + (', QC passed only' if self.qc.isChecked() else '')),
            ('group statistic', self.unit.currentText()))]

    def ask_path(self, stem):
        start = os.path.join(self.root, f'{stem}.xlsx')
        p, _ = QFileDialog.getSaveFileName(self, 'Export to Excel', start, 'Excel (*.xlsx)')
        return p if not p or p.lower().endswith('.xlsx') else p + '.xlsx'

    def write(self, path, sheets):
        try:
            pipeline.write_xlsx(path, sheets + [('info', self.info_sheet())])
        except OSError as e:
            QMessageBox.warning(self, 'Export', f'Cannot write {path}:\n{e}'); return
        QMessageBox.information(self, 'Export', f'Saved {path}')
        QDesktopServices.openUrl(QUrl.fromLocalFile(os.path.dirname(path)))

    def export_table(self):
        if not self.data:
            return
        sc = self.scope.currentData()
        p = self.ask_path(f'results_{sc}')
        if not p:
            return
        # rows in the order shown (after sorting by a header)
        order = [self.proxy.mapToSource(self.proxy.index(i, 0)).row() for i in range(self.proxy.rowCount())]
        rows = [self.model.rows[i] for i in order]
        self.write(p, [(sc, self.sheet(rows, self.model.cols, self.model.headers))])

    def export_all(self):
        if not self.data:
            return
        p = self.ask_path('results_all_scopes')
        if not p:
            return
        sheets = []
        for sc in ('cells', 'mito', 'puncta'):
            if self.data.get(sc):
                sheets.append((sc, self.sheet(*self.table(sc)[:3])))
        for sc in ('images', 'groups'):
            for src, l in SOURCES:
                if not self.data.get(src):
                    continue
                sheets.append((f'{sc}_{src}'[:31], self.sheet(*self.table(sc, src)[:3])))
        self.write(p, sheets)


def open_browser(folder='', parent=None):
    w = ResultsBrowser(folder, parent)
    w.setAttribute(Qt.WA_DeleteOnClose)
    w.show()
    return w


def main(folder=''):
    import sys
    from PySide6.QtWidgets import QApplication as App
    app = App.instance() or App(sys.argv)
    app.setApplicationName('Mito Analyzer')
    w = ResultsBrowser(folder); w.show()
    return app.exec()
