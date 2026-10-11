#!/usr/bin/env bash
# Build the documentation site in CI: pinned Doxygen, a venv with
# doc/requirements.txt, then `make -C doc html` (warnings are errors).
#
# Shared by .github/workflows/docs.yml (PR check) and docs-deploy.yml (Pages).
# Uses only the runner's system python3 and curl, so the workflows need no
# marketplace actions beyond actions/checkout.

set -euo pipefail

DOXYGEN_VERSION="1.12.0"
DOXYGEN_SHA256="3c42c3f3fb206732b503862d9c9c11978920a8214f223a3950bbf2520be5f647"

DOC_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="${RUNNER_TEMP:-${TMPDIR:-/tmp}}/moris-docs-tools"
mkdir -p "$TOOLS_DIR"

# Doxygen, pinned: Breathe output (and so the -W build) depends on its version.
if [ ! -x "$TOOLS_DIR/doxygen/bin/doxygen" ]; then
    curl -sSLf -o "$TOOLS_DIR/doxygen.tgz" \
        "https://www.doxygen.nl/files/doxygen-${DOXYGEN_VERSION}.linux.bin.tar.gz"
    echo "${DOXYGEN_SHA256}  $TOOLS_DIR/doxygen.tgz" | sha256sum --check
    mkdir -p "$TOOLS_DIR/doxygen"
    tar xzf "$TOOLS_DIR/doxygen.tgz" -C "$TOOLS_DIR/doxygen" --strip-components=1
    rm "$TOOLS_DIR/doxygen.tgz"
fi
"$TOOLS_DIR/doxygen/bin/doxygen" --version

# Sphinx toolchain (+ pillow/pytest for the gallery tests)
python3 -m venv "$TOOLS_DIR/venv"
"$TOOLS_DIR/venv/bin/pip" install --quiet --upgrade pip
"$TOOLS_DIR/venv/bin/pip" install --quiet -r "$DOC_DIR/requirements.txt" pillow pytest

"$TOOLS_DIR/venv/bin/python" -m pytest -q "$DOC_DIR/tests"

# Example-gallery scenes: the gallery-assets release asset (doc/gallery/build_gallery.py --upload).
# In CI only a missing release is tolerated; any other failure must not deploy a gallery-less site.
# Locally the step is skipped unless MORIS_DOCS_GALLERY=download (an unpacked doc/_static/gallery is kept).
if [ -n "${GITHUB_ACTIONS:-}" ] || [ "${MORIS_DOCS_GALLERY:-}" = "download" ]; then
    rm -rf "$DOC_DIR/_static/gallery"
    mkdir -p "$DOC_DIR/_static/gallery"
    if view=$(gh release view gallery-assets 2>&1); then
        gh release download gallery-assets -p gallery-assets.zip -D "$TOOLS_DIR" --clobber
        unzip -q -o "$TOOLS_DIR/gallery-assets.zip" -d "$DOC_DIR/_static/gallery"
        echo "gallery: $(ls "$DOC_DIR/_static/gallery" | grep -c '\.vtksz$') scenes"
    elif echo "$view" | grep -q "release not found"; then
        echo "::warning::gallery-assets release not found; building without the example scenes"
    else
        echo "$view" >&2
        exit 1
    fi
fi

make -C "$DOC_DIR" html \
    DOXYGEN="$TOOLS_DIR/doxygen/bin/doxygen" \
    SPHINXBUILD="$TOOLS_DIR/venv/bin/sphinx-build"

echo "Doxygen warnings: $(wc -l < "$DOC_DIR/_build/doxygen/doxygen-warnings.log") lines (not fatal)"
