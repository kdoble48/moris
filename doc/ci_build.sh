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

# Sphinx toolchain
python3 -m venv "$TOOLS_DIR/venv"
"$TOOLS_DIR/venv/bin/pip" install --quiet --upgrade pip
"$TOOLS_DIR/venv/bin/pip" install --quiet -r "$DOC_DIR/requirements.txt"

make -C "$DOC_DIR" html \
    DOXYGEN="$TOOLS_DIR/doxygen/bin/doxygen" \
    SPHINXBUILD="$TOOLS_DIR/venv/bin/sphinx-build"

echo "Doxygen warnings: $(wc -l < "$DOC_DIR/_build/doxygen/doxygen-warnings.log") lines (not fatal)"
