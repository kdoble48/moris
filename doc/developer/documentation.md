# Writing documentation

This site is built with [Sphinx](https://www.sphinx-doc.org). The prose pages are Markdown
([MyST](https://myst-parser.readthedocs.io)); the C++ API pages are generated from Doxygen
comments in the headers and rendered through [Breathe](https://breathe.readthedocs.io).

## Building the site

```bash
python3 -m venv ~/.venvs/moris-docs
~/.venvs/moris-docs/bin/pip install -r doc/requirements.txt
make -C doc html SPHINXBUILD=~/.venvs/moris-docs/bin/sphinx-build
```

Open `doc/_build/html/index.html`. The build needs `doxygen` on the `PATH`, or pass
`DOXYGEN=/path/to/doxygen`. CI pins Doxygen 1.12.0; other versions generally work, but because
warnings are errors, a different version can surface warnings CI does not.

| Target | What it does |
|---|---|
| `make -C doc html` | Doxygen (XML and HTML), then Sphinx; the Doxygen HTML is copied to `_build/html/doxygen/` |
| `make -C doc html-nodox` | Sphinx only — fast, for editing prose; API directives render as placeholders |
| `make -C doc clean` | remove `doc/_build` |

The Sphinx step runs with warnings as errors, as CI does. From a configured CMake tree,
`cmake --build build_opt --target docs` runs the same build when Sphinx and Doxygen are found.

## Publishing

Two workflows share the build script `doc/ci_build.sh`, which installs the pinned Doxygen and the
Python requirements and runs `make -C doc html`:

- `.github/workflows/docs.yml` builds the site as a check on every pull request that touches
  `doc/` or a header. It uses no action besides `actions/checkout`.
- `.github/workflows/docs-deploy.yml` builds and publishes the site to GitHub Pages on every
  push to `develop`. It needs Pages enabled with source "GitHub Actions" and the Pages actions
  on the repository's Actions allow-list.

## Adding a page

1. Write `doc/<section>/<page>.md`, starting with a single `#` title.
2. Add it to the `toctree` of the section's `index.md` (or, for a new section, to
   `doc/index.md`).
3. Link to other pages with `` {doc}`../user-guide/interlinks` `` and to headings with
   `[text](../user-guide/interlinks.md#heading-anchor)` (anchors are generated for `#` to
   `###` headings).

### Example pages

Every example under `projects/EXA` has a page at `doc/examples/exa/<category>/<Example>.md`, where
`<category>` is the first directory level (`optimization`, `structure`, `thermal`, `fluid`,
`mesh_generation`) and `<Example>` the example's directory name. A page opens with a summary
paragraph and the generated facts block,

````md
```{exa-summary} projects/EXA/structure/linear/Plane_Strain
```
````

which lists the decks, the ctest names, and the features detected in the deck source. It
continues with **Problem**, **Setup in MORIS** (one paragraph per configured module), **What the
test checks**, and **Run it**; `doc/examples/exa/optimization/LevelSet_Boxbeam_NewIO.md` is the
reference. Every value on the page comes from the deck, its `example_test_case.cpp`, or its
`CMakeLists.txt`. A new example without a page still appears in the catalog, linked to its
source directory.

### Regenerating the example gallery

Each example page carries an interactive 3D scene of the example's own output, and the
{doc}`../examples/example-gallery` shows them as a card grid. The scenes are static vtk.js files
(`.vtksz`) shown in an iframe by trame-vtk's standalone `viewer.html`; no server runs anywhere. They
are generated locally — they need ParaView and a build tree in which the examples have been run — and
stored outside git, as the `gallery-assets` release asset that the docs workflows download.

```bash
# 1. run the examples in a build configured with MORIS_USE_EXAMPLES=ON
ctest --test-dir build_opt -R '<pattern>' --output-on-failure

# 2. export one scene per example (ParaView 6 pvpython; vendors trame-vtk into doc/gallery/_vendor)
python3 doc/gallery/build_gallery.py --build build_opt --out $MORIS_RUNS_DIR/docs_gallery

# 3. check the thumbnails in the output directory, then publish
python3 doc/gallery/build_gallery.py --build build_opt --out $MORIS_RUNS_DIR/docs_gallery --upload
```

| File | Role |
|---|---|
| `doc/gallery/gallery_select.py` | picks the file and field per example: a VIS result field, else the XTK integration mesh by phase, else a GEN field, else the bare mesh |
| `doc/gallery/export_scenes.py` | pvpython: renders the choice and writes `<id>.png`, `<id>-thumb-raw.png`, `<id>.vtksz`, `<id>.json`; for a series of three or more steps (an `.e-s.NNNN` optimization series, or a transient run) also one frame per step, up to 60, under one color range |
| `doc/gallery/scene_export.py` | serializes a render view with trame-vtk (partitioned-dataset and categorical-color workarounds) |
| `doc/gallery/build_gallery.py` | one pvpython process per example, thumbnails, the animated WebP from the frames, scene validation, `manifest.json`, zip, upload |
| `doc/gallery/overrides.json` | optional per-example choices: `{"<Example>": {"file": …, "array": …, "skip": true}}` |
| `doc/_ext/moris_gallery.py` | Sphinx: scenes on the example pages, the gallery cards, `viewer.html` |

`build_gallery.py` runs Windows ParaView 6 through WSL when it is installed, otherwise the `pvpython`
on the `PATH` (`--pvpython` overrides). To preview without uploading, unzip the archive into
`doc/_static/gallery/` and build the site; CI downloads the release instead
(`MORIS_DOCS_GALLERY=download bash doc/ci_build.sh` does the same locally). Each upload changes the
manifest's timestamp, which the pages append to the asset URLs, so browsers and the GitHub Pages
cache load the new scenes without renamed files.

Design records — RFCs, reviews, plans — go in `doc/internal/` and are listed on the
{doc}`../design/index` page.

## Documenting C++ code

Doxygen reads every header under `projects/*/src`. Comments in Javadoc style appear in the API
reference; the first sentence is the brief description:

```cpp
/**
 * Computes the residual of the linear elastic bulk weak form.
 *
 * Assembles the internal-force contribution of the current integration cell into the
 * leader residual, using the constitutive model registered as "ElastLinIso".
 *
 * @param aWStar  integration weight times the Jacobian determinant
 */
void compute_residual( real aWStar ) override;
```

Use `//!` for one-line member comments. Formulas go in `\f$ … \f$` (inline) or `\f[ … \f]`
(display) and render with MathJax.

To feature a class on a curated API page, add a Breathe directive to the page under `doc/api/`:

````md
```{doxygenclass} moris::fem::IWG
:members: compute_residual, compute_jacobian
```
````

The fully qualified name must match Doxygen's; search
`doc/_build/doxygen/xml/index.xml` when unsure.
