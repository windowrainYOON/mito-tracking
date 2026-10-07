#!/usr/bin/env bash
# One-line installer / updater for Mito Analyzer (macOS):
#   curl -fsSL https://raw.githubusercontent.com/windowrainYOON/mito-tracking/main/install.sh | bash
# Downloads the code (no git or Xcode needed), finds Python 3.11+ or installs one in your home folder
# without admin rights (uv), builds "Mito Analyzer.app" and pins it to the Dock. Run it again to update.
# Settings (environment variables):
#   MITO_DIR=~/Applications/MitoAnalyzer   where the code and the app go
#   MITO_REF=main                          branch or tag to install
#   MITO_DOCK=1                            0 = do not pin to the Dock
set -euo pipefail

REPO="windowrainYOON/mito-tracking"
REF="${MITO_REF:-main}"
DEST="${MITO_DIR:-$HOME/Applications/MitoAnalyzer}"
DOCK="${MITO_DOCK:-1}"

say() { echo; echo "==> $*"; }
die() { echo "Error: $*" >&2; exit 1; }

if [ "$(uname -s)" != Darwin ] && [ -z "${MITO_SKIP_BUILD:-}" ]; then
  die "this installer builds the macOS app. On Windows / Linux run from source (see README)."
fi
command -v curl >/dev/null || die "curl is needed"
command -v tar >/dev/null || die "tar is needed"

# 1. code: download the branch as a tarball and swap it in, keeping the build environment for faster updates
say "Downloading $REPO ($REF)"
TMP="$(mktemp -d)"; trap 'rm -rf "$TMP"' EXIT
curl -fsSL "https://codeload.github.com/$REPO/tar.gz/$REF" -o "$TMP/src.tgz" \
  || die "could not download https://github.com/$REPO ($REF)"
mkdir -p "$TMP/src"
tar -xzf "$TMP/src.tgz" -C "$TMP/src" --strip-components 1
[ -f "$TMP/src/application/build_mac.sh" ] || die "the download does not contain application/build_mac.sh"
mkdir -p "$(dirname "$DEST")"
if [ -d "$DEST/application/.venv" ]; then
  mv "$DEST/application/.venv" "$TMP/src/application/.venv"
fi
rm -rf "$DEST"
mv "$TMP/src" "$DEST"
echo "Code in $DEST"

# 2. Python 3.11+: use one that is there, else install uv + Python 3.12 into ~/.local (no admin rights)
ok_py() { command -v "$1" >/dev/null 2>&1 && "$1" -c 'import sys; sys.exit(sys.version_info < (3, 11))' 2>/dev/null; }
PY=""
for c in ${PYTHON:-} python3.13 python3.12 python3.11 python3 "$HOME"/.local/bin/python3.1[1-9] \
         /opt/homebrew/bin/python3 /usr/local/bin/python3 /Library/Frameworks/Python.framework/Versions/3.1[1-9]/bin/python3; do
  if ok_py "$c"; then PY="$c"; break; fi
done
UV="$(command -v uv || true)"; [ -z "$UV" ] && [ -x "$HOME/.local/bin/uv" ] && UV="$HOME/.local/bin/uv"
if [ -z "$PY" ] && [ -n "$UV" ]; then
  PY="$("$UV" python find '>=3.11' 2>/dev/null || true)"
fi
if [ -z "$PY" ]; then
  say "No Python 3.11+ found; installing uv and Python 3.12 into your home folder (no admin rights needed)"
  if [ -z "$UV" ]; then
    curl -LsSf https://astral.sh/uv/install.sh | env UV_NO_MODIFY_PATH=1 sh
    UV="$HOME/.local/bin/uv"
  fi
  "$UV" python install 3.12
  PY="$("$UV" python find 3.12)"
fi
ok_py "$PY" || die "could not find or install Python 3.11+"
echo "Python: $("$PY" --version) ($PY)"

# 3. build (and pin to the Dock)
if [ -n "${MITO_SKIP_BUILD:-}" ]; then
  echo "MITO_SKIP_BUILD set: not building"; exit 0
fi
ARGS=()
[ "$DOCK" = 1 ] && ARGS+=(--add-to-dock)
PYTHON="$PY" "$DEST/application/build_mac.sh" ${ARGS[@]+"${ARGS[@]}"}

APP="$DEST/application/Mito Analyzer.app"
say "Installed: $APP"
echo "Update later by running the same command again."
open "$APP" 2>/dev/null || true
