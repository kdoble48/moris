# Sphinx configuration for the MORIS documentation site.
#
# Build with `make -C doc html` (see doc/Makefile). Prose is Markdown (MyST);
# the C++ API pages pull Doxygen XML through Breathe.

import os
import re
import sys
from pathlib import Path

DOC_DIR = Path(__file__).resolve().parent
SOURCE_DIR = Path(os.environ.get("MORIS_SOURCE_DIR", DOC_DIR.parent))
DOXYGEN_XML = Path(os.environ.get("DOXYGEN_OUTPUT_DIR", DOC_DIR / "_build" / "doxygen")) / "xml"
USE_DOXYGEN = not os.environ.get("MORIS_DOCS_NO_DOXYGEN") and (DOXYGEN_XML / "index.xml").exists()

sys.path.insert(0, str(DOC_DIR / "_ext"))


def _moris_version():
    """Read the version from the top-level CMakeLists.txt project() call."""
    text = (SOURCE_DIR / "CMakeLists.txt").read_text(errors="replace")
    match = re.search(r"project\s*\(\s*moris[^)]*VERSION\s+([0-9.]+)", text, re.IGNORECASE)
    return match.group(1) if match else ""


# -- Project -----------------------------------------------------------------

project = "MORIS"
author = "University of Colorado AMReC, MORIS developers"
copyright = "2022-2026, University of Colorado"
release = _moris_version()
version = release

# -- General -----------------------------------------------------------------

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinx.ext.githubpages",
    "sphinx_copybutton",
    "moris_sphinx",
]
if USE_DOXYGEN:
    extensions.append("breathe")

# moris_sphinx: render placeholder notes for Breathe directives without Doxygen XML
moris_doxygen_available = USE_DOXYGEN

source_suffix = {".md": "markdown", ".rst": "restructuredtext"}
root_doc = "index"
exclude_patterns = ["_build", "_ext", "Thumbs.db", ".DS_Store", "README.md"]

myst_enable_extensions = ["colon_fence", "deflist", "dollarmath", "fieldlist", "attrs_inline"]
myst_heading_anchors = 3

# The design records in doc/internal/ are working notes written for GitHub:
# their relative links point at repository files, not at site pages.
suppress_warnings = ["myst.xref_missing", "myst.header"]

# -- Breathe -----------------------------------------------------------------

breathe_projects = {"moris": str(DOXYGEN_XML)}
breathe_default_project = "moris"
breathe_default_members = ("members",)
breathe_domain_by_extension = {"hpp": "cpp", "h": "cpp"}
breathe_show_include = True
cpp_index_common_prefix = ["moris::"]

# -- HTML --------------------------------------------------------------------

html_theme = "furo"
html_title = f"MORIS {release}".strip()
html_static_path = ["_static"]
html_css_files = ["moris.css"]
html_baseurl = os.environ.get("MORIS_DOCS_BASEURL", "")
html_theme_options = {
    "source_repository": "https://github.com/kdoble48/moris",
    "source_branch": "develop",
    "source_directory": "doc/",
    "navigation_with_keys": True,
}
