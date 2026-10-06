# PyInstaller spec: builds "Mito Analyzer.app" on macOS (a plain folder app elsewhere).
# Usage: pyinstaller --noconfirm --clean MitoAnalyzer.spec   (build_mac.sh does this for you)
from PyInstaller.utils.hooks import collect_data_files, collect_submodules

hidden = collect_submodules('skimage') + collect_submodules('tifffile') + ['roifile']
datas = collect_data_files('skimage')
excludes = ['tkinter', 'PySide6.QtWebEngineCore', 'PySide6.QtWebEngineWidgets', 'PySide6.Qt3DCore',
            'PySide6.QtQuick', 'PySide6.QtQml', 'PySide6.QtMultimedia', 'PySide6.QtCharts', 'PySide6.QtPdf',
            'IPython', 'pytest', 'dask', 'pandas']

a = Analysis(['run_app.py'], pathex=['.'], hiddenimports=hidden, datas=datas, excludes=excludes)
pyz = PYZ(a.pure)
exe = EXE(pyz, a.scripts, [], exclude_binaries=True, name='MitoAnalyzer', console=False,
          argv_emulation=False, target_arch=None)
coll = COLLECT(exe, a.binaries, a.datas, name='MitoAnalyzer')
app = BUNDLE(coll, name='Mito Analyzer.app', bundle_identifier='com.mitotracking.mitoanalyzer',
             info_plist={'NSHighResolutionCapable': True, 'CFBundleShortVersionString': '0.1.0',
                         'NSRequiresAquaSystemAppearance': False})
