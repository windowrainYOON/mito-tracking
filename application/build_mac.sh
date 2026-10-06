#!/usr/bin/env bash
# Build "Mito Analyzer.app" into this folder (application/).
# Needs Python 3.11+ (python.org installer or `brew install python@3.12`).
# Usage: ./build_mac.sh [--add-to-dock]   -> application/Mito Analyzer.app
#   --add-to-dock  pin the app to the Dock (only if it is not there yet); rebuilds keep the same path
set -euo pipefail
cd "$(dirname "$0")"

PY=""
for c in ${PYTHON:-} python3.13 python3.12 python3.11 python3; do
  if command -v "$c" >/dev/null 2>&1 && "$c" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; then
    PY="$c"; break
  fi
done
if [ -z "$PY" ]; then
  echo "Python 3.11 or newer is required. Install it with 'brew install python@3.12' or from python.org." >&2
  exit 1
fi
echo "Using $($PY --version) at $(command -v "$PY")"

"$PY" -m venv .venv
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt pyinstaller
.venv/bin/python -m PyInstaller --noconfirm --clean MitoAnalyzer.spec

rm -rf "Mito Analyzer.app"
cp -R "dist/Mito Analyzer.app" "Mito Analyzer.app"
# unsigned local build: drop the quarantine flag so Gatekeeper does not block it
xattr -cr "Mito Analyzer.app" 2>/dev/null || true
APP="$(pwd)/Mito Analyzer.app"
if [ "${1:-}" = "--add-to-dock" ]; then
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
echo
echo "Built: $(pwd)/Mito Analyzer.app"
echo "Open it with: open \"$(pwd)/Mito Analyzer.app\""
