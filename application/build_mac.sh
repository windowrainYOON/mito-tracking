#!/usr/bin/env bash
# Build "Mito Analyzer.app" into this folder (application/).
# Usage: ./build_mac.sh [--add-to-dock] [--locked] [--clean] [--dmg]
#   --add-to-dock  pin the app to the Dock (only if it is not there yet); rebuilds keep the same path
#   --locked       install the exact versions in requirements-lock.txt instead of the newest releases
#   --clean        recreate the build environment (.venv) from scratch
#   --dmg          also pack the app into "Mito Analyzer-<version>.dmg" to give to others (not tracked by git)
# Needs Python 3.11 or newer. If none is found and Homebrew is installed, python@3.12 is installed with it.
set -euo pipefail
cd "$(dirname "$0")"

DOCK=0; REQ=requirements.txt; CLEAN=0; DMG=0
for a in "$@"; do
  case "$a" in
    --add-to-dock) DOCK=1 ;;
    --locked) REQ=requirements-lock.txt ;;
    --clean) CLEAN=1 ;;
    --dmg) DMG=1 ;;
    -h|--help) sed -n '2,8p' "$0"; exit 0 ;;
    *) echo "Unknown option: $a (see --help)" >&2; exit 2 ;;
  esac
done

find_python() {
  local uvpy=""
  command -v uv >/dev/null 2>&1 && uvpy="$(uv python find '>=3.11' 2>/dev/null || true)"
  [ -z "$uvpy" ] && [ -x "$HOME/.local/bin/uv" ] && uvpy="$("$HOME/.local/bin/uv" python find '>=3.11' 2>/dev/null || true)"
  for c in ${PYTHON:-} python3.13 python3.12 python3.11 python3 $uvpy \
           "$HOME"/.local/bin/python3.1[1-9] /opt/homebrew/bin/python3 /usr/local/bin/python3 \
           /Library/Frameworks/Python.framework/Versions/3.1[1-9]/bin/python3; do
    if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
      echo "$c"; return 0
    fi
  done
  return 1
}

step() { echo; echo "==> $*"; }

PY="$(find_python || true)"
if [ -z "$PY" ]; then
  if command -v brew >/dev/null 2>&1; then
    step "Python 3.11+ not found; installing python@3.12 with Homebrew"
    brew install python@3.12
    PY="$(find_python || true)"
  fi
  if [ -z "$PY" ]; then
    cat >&2 <<'MSG'
Python 3.11 or newer is required. Either
  - install it from https://www.python.org/downloads/macos/ , or
  - without admin rights:  curl -LsSf https://astral.sh/uv/install.sh | sh && ~/.local/bin/uv python install 3.12
then run this script again.
MSG
    exit 1
  fi
fi
step "Using $("$PY" --version) at $(command -v "$PY")"

# rebuild the environment if asked, if it is missing, or if it was made with a Python that is gone or too old
if [ "$CLEAN" = 1 ] || ! .venv/bin/python -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
  rm -rf .venv
  "$PY" -m venv .venv
fi

LOG="$(pwd)/build.log"; : > "$LOG"
fail() { echo "$1 (last lines of $LOG below)" >&2; tail -40 "$LOG" >&2; exit 1; }
step "Installing packages ($REQ)"
.venv/bin/python -m pip install --upgrade pip >>"$LOG" 2>&1 || fail "pip upgrade failed"
.venv/bin/python -m pip install --upgrade -r "$REQ" pyinstaller >>"$LOG" 2>&1 || fail "Package install failed"
.venv/bin/python - <<'EOF'
import importlib
for m in ('numpy', 'scipy', 'skimage', 'tifffile', 'roifile', 'matplotlib', 'PySide6', 'openpyxl', 'PyInstaller'):
    print(f'  {m:11s} {importlib.import_module(m).__version__}')
EOF

step "Checking the app starts"
.venv/bin/python -c 'from mito_app import pipeline, gui' || {
  echo "The app does not import with these package versions. Try: ./build_mac.sh --locked" >&2; exit 1; }

step "Building (takes a few minutes)"
.venv/bin/python -m PyInstaller --noconfirm --clean --log-level WARN MitoAnalyzer.spec >>"$LOG" 2>&1 \
  || fail "PyInstaller build failed"

rm -rf "Mito Analyzer.app"
cp -R "dist/Mito Analyzer.app" "Mito Analyzer.app"
rm -rf build dist
# unsigned local build: drop the quarantine flag so Gatekeeper does not block it
xattr -cr "Mito Analyzer.app" 2>/dev/null || true
APP="$(pwd)/Mito Analyzer.app"
if [ "$DOCK" = 1 ]; then
  if defaults read com.apple.dock persistent-apps 2>/dev/null | grep -q "Mito%20Analyzer.app"; then
    echo "Already in the Dock."
  else
    defaults write com.apple.dock persistent-apps -array-add \
      "<dict><key>tile-data</key><dict><key>file-data</key><dict><key>_CFURLString</key><string>file://${APP// /%20}/</string><key>_CFURLStringType</key><integer>15</integer></dict></dict></dict>"
    killall Dock
    echo "Added to the Dock."
  fi
fi
# refresh Finder/Dock icon cache for the rebuilt bundle
touch "$APP"
if [ "$DMG" = 1 ]; then
  VER="$(.venv/bin/python -c 'import mito_app; print(mito_app.__version__)')"
  DMGF="$(pwd)/Mito Analyzer-$VER.dmg"
  step "Packing $DMGF"
  STAGE="$(mktemp -d)"
  cp -R "$APP" "$STAGE/"
  ln -s /Applications "$STAGE/Applications"
  rm -f "$DMGF"
  hdiutil create -volname "Mito Analyzer $VER" -srcfolder "$STAGE" -fs HFS+ -format UDZO -ov "$DMGF" >>"$LOG" 2>&1 \
    || { rm -rf "$STAGE"; fail "hdiutil failed"; }
  rm -rf "$STAGE"
  echo "DMG: $DMGF ($(du -h "$DMGF" | cut -f1)). Open it and drag the app onto Applications."
  echo "Unsigned: on another Mac, the first launch needs System Settings > Privacy & Security > Open Anyway."
fi
step "Built: $APP"
echo "Open it with: open \"$APP\""
