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
