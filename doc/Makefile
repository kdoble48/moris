# Build the MORIS documentation site.
#
#   make -C doc html        Doxygen XML/HTML, then Sphinx  ->  doc/_build/html/index.html
#   make -C doc html-nodox  Sphinx only (API pages render as stubs; fast prose iteration)
#   make -C doc clean
#
# Tools: pip install -r doc/requirements.txt; doxygen >= 1.9 on PATH (or DOXYGEN=...).

DOCDIR       := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
SRCDIR       := $(abspath $(DOCDIR)/..)
BUILDDIR     ?= $(DOCDIR)/_build
DOXYGEN      ?= doxygen
SPHINXBUILD  ?= sphinx-build
SPHINXOPTS   ?= -W --keep-going -j auto

export MORIS_SOURCE_DIR   := $(SRCDIR)
export DOXYGEN_OUTPUT_DIR := $(BUILDDIR)/doxygen

.PHONY: html html-nodox doxygen clean

html: doxygen
	$(SPHINXBUILD) -b html $(SPHINXOPTS) $(DOCDIR) $(BUILDDIR)/html
	rm -rf $(BUILDDIR)/html/doxygen
	cp -r $(DOXYGEN_OUTPUT_DIR)/html $(BUILDDIR)/html/doxygen

html-nodox:
	MORIS_DOCS_NO_DOXYGEN=1 $(SPHINXBUILD) -b html $(SPHINXOPTS) $(DOCDIR) $(BUILDDIR)/html

doxygen:
	mkdir -p $(DOXYGEN_OUTPUT_DIR)
	cd $(DOCDIR) && $(DOXYGEN) Doxyfile

clean:
	rm -rf $(BUILDDIR)
