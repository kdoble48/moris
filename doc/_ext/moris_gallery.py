"""Sphinx extension: interactive 3D scenes on the example pages, and the example gallery.

On ``builder-inited`` it reads ``_static/gallery/manifest.json`` (downloaded from the
``gallery-assets`` release by CI, or written by ``doc/gallery/build_gallery.py``), writes the card
grid ``examples/_generated/example_gallery.md``, and copies trame-vtk's standalone viewer to
``_static/viewer.html``. The ``{exa-summary}`` directive (exa_catalog.py) calls :func:`scene_markdown`
to put the Interactive / Static tabs on each example page. Without a manifest the pages render
without scenes and the gallery page says how to generate them.

The scenes are static vtk.js files shown by the viewer in an iframe: no trame server runs anywhere.
The page functions don't import Sphinx, so ``doc/tests/test_gallery_pages.py`` calls them directly.
Adapted from kdoble48/nastran-pv (docs/_ext/npv_gallery.py).
"""
from __future__ import annotations

import html
import json
import re
import shutil
from pathlib import Path
from urllib.parse import quote

ASSETS = "_static/gallery"
GENERATED = Path("examples") / "_generated" / "example_gallery.md"
COLLAPSE_OVER = 12
CATEGORY_TITLES = {
    "optimization": "Optimization and sensitivity analysis",
    "structure": "Structural mechanics",
    "thermal": "Thermal",
    "fluid": "Fluid dynamics",
    "mesh_generation": "Mesh generation",
}

_STATE = {"manifest": None, "by_directory": {}, "version": ""}


def hexrgb(rgb) -> str:
    return "#" + "".join(f"{max(0, min(255, round(float(c) * 255))):02x}" for c in rgb)


def _num(x) -> str:
    return f"{float(x):.4g}"


def legend_html(leg) -> str:
    title = f'<div class="moris-legend-title">{html.escape(leg["array"])}</div>'
    if leg["kind"] == "continuous":
        lo, hi = (float(x) for x in leg["range"])
        span = (hi - lo) or 1.0
        grad = ", ".join(f"{hexrgb(s[1:4])} {100 * (float(s[0]) - lo) / span:.1f}%" for s in leg["stops"])
        return (f'<div class="moris-legend">{title}'
                f'<div class="moris-bar" style="background: linear-gradient(90deg, {grad})"></div>'
                f'<div class="moris-bar-range"><span>{_num(lo)}</span><span>{_num(hi)}</span></div></div>')
    items = "".join(f'<li><span class="moris-swatch" style="background:{hexrgb(e["rgb"])}"></span>'
                    f'{html.escape(str(e["label"]))}</li>' for e in leg["entries"])
    body = f'<ul class="moris-cats">{items}</ul>'
    if len(leg["entries"]) > COLLAPSE_OVER:
        body = f'<details><summary>{len(leg["entries"])} values</summary>{body}</details>'
    return f'<div class="moris-legend">{title}{body}</div>'


def _tag(version) -> str:
    return f"?v={version}" if version else ""


def root_prefix(docname: str) -> str:
    """Relative path from a page's HTML file to the site root (``examples/exa/fluid/X`` -> ``../../../``)."""
    return "../" * docname.count("/")


def _size(nbytes) -> str:
    return f"{nbytes / 1e6:.1f} MB" if nbytes >= 1_000_000 else f"{max(1, round(nbytes / 1e3))} KB"


def scene_markdown(meta, docname, version="", build="") -> str:
    """Interactive / Static tabs, legends and provenance for one example page (MyST markdown)."""
    prefix = root_prefix(docname)
    # the viewer fetches fileURL relative to itself; the per-build tag defeats GitHub Pages' 10-minute cache
    url = quote(f"gallery/{meta['scene']}{_tag(version)}", safe="/")
    iframe = (f'<iframe class="moris-scene" title="3D scene: {html.escape(meta["title"])}" loading="lazy" '
              f'src="{prefix}_static/viewer.html?fileURL={url}"></iframe>')
    image = (f'<img class="moris-static" src="{prefix}{ASSETS}/{meta["png"]}{_tag(version)}" '
             f'alt="{html.escape(meta["caption"])}: {html.escape(meta["title"])}">')
    hint = f'Drag to rotate, scroll to zoom, right-drag to pan. Scene: {_size(meta["scene_bytes"])}.'
    legends = "".join(legend_html(lg) for lg in meta.get("legends", []))
    summary = "\n".join(t for t in meta.get("text", []) if t).strip()
    side = legends + (f'<pre class="moris-summary">{html.escape(summary)}</pre>' if summary else "")
    source = meta.get("source", {})
    run = {"showcase": " of a full-length showcase run (`MORIS_EXA_MODE=showcase`)",
           "regression": " of the regression test run"}.get(source.get("run"), "")
    provenance = (f'{meta["caption"]}. Rendered from `{source.get("file", "?")}`{run}'
                  f'{", written " + source["modified"] if source.get("modified") else ""}'
                  f'{" by " + build if build else ""}.')
    animation = ""
    anim = meta.get("animation")
    if anim:
        frames = len(anim["frames"])
        if anim.get("iterations"):
            what = "optimization iterations"
            per = anim["steps"] // anim["iterations"]
            sampled = (f'{frames} frames over {anim["iterations"]} {what}'
                       + (f' ({per} time steps each)' if per > 1 else ""))
        else:
            what = "time steps"
            sampled = f'{frames} of {anim["steps"]} {what}' if frames < anim["steps"] else f'all {anim["steps"]} {what}'
        anim_legend = (f'```{{raw}} html\n<div class="moris-side">{legend_html(anim["legend"])}</div>\n```\n\n'
                       if anim.get("legend") else "")
        animation = (f'````{{tab-item}} Animation\n'
                     f'<img class="moris-static" src="{prefix}{ASSETS}/{anim["file"]}{_tag(version)}" '
                     f'alt="Animation over the {what}: {html.escape(meta["title"])}">\n\n'
                     f'{anim_legend}{sampled}, one color range across all frames.\n````\n')
    paraview = paraview_markdown(meta, _STATE.get("manifest") or {})
    note = ""
    if meta.get("kind") == "xtk":
        note = ("```{note}\nThis is the XTK **integration mesh**, not the background (HMR) mesh: each cell the "
                "interface crosses is subdivided into triangles that conform to the interface, and uncut cells "
                "are kept whole.\n```\n\n")
    return (f'`````{{tab-set}}\n'
            f'````{{tab-item}} Interactive\n{iframe}\n\n{hint}\n````\n'
            f'{animation}'
            f'````{{tab-item}} Static\n{image}\n````\n'
            f'`````\n\n'
            + note
            + (f'```{{raw}} html\n<div class="moris-side">{side}</div>\n```\n\n' if side else "")
            + f'{provenance}{paraview}\n')


def paraview_markdown(meta, manifest) -> str:
    """The "Open in ParaView" sentence: a link to the example's bundle (a release asset), or why there is none."""
    bundle = meta.get("paraview") or {}
    base = manifest.get("paraview_url", "")
    if bundle.get("file") and base:
        size = _size(bundle.get("bytes", 0))
        return (f' **[Open in ParaView]({base}{bundle["file"]})** ({size}): the output file, the ParaView state and '
                f'an `open.py`; run `pvpython open.py` in the unpacked folder, or load `state.pvsm` in ParaView and '
                f'point it at the folder.')
    if bundle.get("skipped"):
        return f' No ParaView bundle for this example: {bundle["skipped"]}.'
    return ""


def gallery_markdown(manifest) -> str:
    """Card grid of every exported example, grouped by category (included by examples/example-gallery.md)."""
    if not manifest:
        return ("The scenes have not been generated for this build. Maintainers: see "
                "{doc}`../developer/documentation` (Regenerating the example gallery).\n")
    out = [f'{len(manifest["views"])} examples, each rendered from the output its ctest run wrote. '
           "Open one to rotate and zoom its scene in the browser.\n"]
    categories = []
    for view in manifest["views"]:
        if view["category"] not in categories:
            categories.append(view["category"])
    for category in categories:
        views = [v for v in manifest["views"] if v["category"] == category]
        out.append(f'## {CATEGORY_TITLES.get(category, category)}\n')
        cards = [f'```{{grid-item-card}} {v["title"]}\n:link: /examples/exa/{v["category"]}/{v["id"]}\n'
                 f':link-type: doc\n:img-top: ../{ASSETS}/{v["thumb"]}\n\n{v["caption"]}\n```' for v in views]
        out.append("`````{grid} 1 2 3 3\n:gutter: 3\n\n" + "\n".join(cards) + "\n`````\n")
    out.append(f'Generated {manifest["generated"][:10]} from docs commit `{manifest["commit"]}`'
               f'{", outputs from " + manifest["build"] if manifest.get("build") else ""}.\n')
    return "\n".join(out)


def version_tag(manifest) -> str:
    """One tag per gallery build: 2026-10-11T14:02:05+00:00 -> 20261011T140205."""
    return re.sub(r"[^0-9T]", "", (manifest or {}).get("generated", ""))[:15]


def scene_for(directory: str, docname: str) -> str:
    """Scene markdown for an example directory, or "" when the manifest has none."""
    meta = _STATE["by_directory"].get(directory)
    if not meta:
        return ""
    return scene_markdown(meta, docname, _STATE["version"], (_STATE["manifest"] or {}).get("build", ""))


def load(srcdir) -> dict | None:
    path = Path(srcdir) / ASSETS / "manifest.json"
    manifest = json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
    _STATE.update(manifest=manifest, version=version_tag(manifest),
                  by_directory={v["directory"]: v for v in (manifest or {}).get("views", [])})
    return manifest


def _builder_inited(app):
    from sphinx.util import logging

    src = Path(app.srcdir)
    manifest = load(src)
    target = src / GENERATED
    target.parent.mkdir(parents=True, exist_ok=True)
    text = gallery_markdown(manifest)
    if not target.exists() or target.read_text(encoding="utf-8") != text:
        target.write_text(text, encoding="utf-8")
    try:
        from trame_vtk.tools.vtksz2html import HTML_VIEWER_PATH
    except ImportError:
        if manifest:
            logging.getLogger(__name__).warning("moris_gallery: trame-vtk not installed; scenes will not load")
        return
    shutil.copy(HTML_VIEWER_PATH, src / "_static" / "viewer.html")


def setup(app):
    app.connect("builder-inited", _builder_inited)
    return {"version": "1", "parallel_read_safe": True, "parallel_write_safe": True}
