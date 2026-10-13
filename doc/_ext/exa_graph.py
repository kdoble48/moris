"""Relational graph of the EXA examples: which examples are alike, and why.

Two examples are related when their decks share detected features (the tags ``exa_catalog``
finds in the sources: physics, optimizer, geometry representation, ghost, Nitsche, dimension,
heat method, adaptive refinement, …). The similarity of a pair is a weighted Jaccard index of
the two tag sets plus a small bonus for sharing a category or a group directory, so the graph
cannot drift from the code any more than the catalog can.

Two renderings, both from :func:`graph_data`:

- ``related_markdown`` — the "Related examples" block on each example page: the example's
  nearest neighbours as an ego graph (``<div class="exa-graph">``, drawn by
  ``_static/exa_graph.js`` with D3) and as a plain list with the shared tags, for readers
  without JavaScript and for the search index;
- ``{exa-graph}`` — the whole graph on one page (``examples/graph.md``), colored by category,
  filterable by tag.

Pure functions (``similarity``, ``neighbours``, ``graph_data``) import nothing from Sphinx:
``doc/tests/test_exa_graph.py`` covers them.
"""

import html
import json
import posixpath

D3_URL = "https://cdnjs.cloudflare.com/ajax/libs/d3/7.9.0/d3.min.js"
# tag weights: dimension and deck style say little about what a deck does; the physics says most
TAG_WEIGHTS = {"2D": 0.5, "3D": 0.5, "Deck API v2": 0.3, "Nitsche": 0.7, "transient": 0.8}
PHYSICS_WEIGHT = 1.5
PHYSICS_TAGS = {"contact", "nonlinear elasticity", "linear elasticity", "diffusion", "advection",
                "incompressible flow", "compressible flow", "turbulence", "Hamilton-Jacobi", "Helmholtz filter",
                "fluid-structure"}
SAME_CATEGORY = 0.15
SAME_GROUP = 0.10
MAX_NEIGHBOURS = 5
MAP_NEIGHBOURS = 3    # the full map stays readable with each example's strongest links only
MIN_SIMILARITY = 0.15


def tag_weight(tag):
    if tag in PHYSICS_TAGS:
        return PHYSICS_WEIGHT
    return TAG_WEIGHTS.get(tag, 1.0)


def similarity(a, b):
    """Weighted Jaccard index of the two examples' tags, plus the category and group bonuses."""
    tags_a, tags_b = set(a.get("tags", [])), set(b.get("tags", []))
    union = tags_a | tags_b
    if not union:
        score = 0.0
    else:
        shared = sum(tag_weight(t) for t in tags_a & tags_b)
        score = shared / sum(tag_weight(t) for t in union)
    if a.get("category") == b.get("category"):
        score += SAME_CATEGORY
        if a.get("group") and a.get("group") == b.get("group"):
            score += SAME_GROUP
    return round(score, 4)


def shared_tags(a, b):
    """The tags both examples carry, in the first example's order."""
    other = set(b.get("tags", []))
    return [t for t in a.get("tags", []) if t in other]


def neighbours(example, examples, limit=MAX_NEIGHBOURS, minimum=MIN_SIMILARITY):
    """The most similar other examples, best first: ``[(example, score), ...]``."""
    scored = [(other, similarity(example, other)) for other in examples if other["leaf"] != example["leaf"]]
    scored = [(other, score) for other, score in scored if score >= minimum]
    scored.sort(key=lambda item: (-item[1], item[0]["title"].lower()))
    return scored[:limit]


def graph_data(examples, urls=None, thumbs=None, centre=None, limit=MAX_NEIGHBOURS):
    """Nodes and edges for D3. Every node keeps its ``limit`` best neighbours; an edge exists when
    either end keeps the other. With ``centre`` (a leaf name) only the centre and its neighbours
    are returned (the ego graph). ``urls``/``thumbs`` map leaf -> page URL / thumbnail URL."""
    urls, thumbs = urls or {}, thumbs or {}
    by_leaf = {e["leaf"]: e for e in examples}
    if centre:
        kept = [by_leaf[centre]] + [other for other, _ in neighbours(by_leaf[centre], examples, limit)]
    else:
        kept = list(examples)
    kept_leaves = {e["leaf"] for e in kept}
    edges = {}
    for example in kept:
        for other, score in neighbours(example, examples, limit):
            if other["leaf"] not in kept_leaves:
                continue
            if centre and centre not in (example["leaf"], other["leaf"]):
                pass  # neighbour-neighbour edges stay: they show which neighbours form a family
            key = tuple(sorted((example["leaf"], other["leaf"])))
            edges.setdefault(key, {"source": key[0], "target": key[1], "weight": score,
                                   "shared": shared_tags(by_leaf[key[0]], by_leaf[key[1]])})
    nodes = [{"id": e["leaf"], "title": e["title"], "category": e["category"], "group": e.get("group", ""),
              "tags": e.get("tags", []), "url": urls.get(e["leaf"], ""), "thumb": thumbs.get(e["leaf"], ""),
              "centre": e["leaf"] == centre}
             for e in sorted(kept, key=lambda e: (e["category"], e["title"].lower()))]
    return {"nodes": nodes, "edges": sorted(edges.values(), key=lambda e: (-e["weight"], e["source"]))}


# --- Sphinx side -------------------------------------------------------------------------------

def _page_docname(example):
    return f"examples/exa/{example['category']}/{example['leaf']}"


def _page_urls(app, examples, docname):
    """Relative URL from ``docname``'s page to each example page that exists (leaf -> URL)."""
    from pathlib import Path

    here = posixpath.dirname(docname)
    urls = {}
    for example in examples:
        target = _page_docname(example)
        if (Path(app.srcdir) / f"{target}.md").exists():
            urls[example["leaf"]] = posixpath.relpath(f"{target}.html", here or ".")
    return urls


def _thumb_urls(examples, docname):
    """Gallery thumbnails for the tooltips, when a gallery manifest is loaded."""
    try:
        import moris_gallery
    except ImportError:
        return {}
    prefix = moris_gallery.root_prefix(docname)
    version = moris_gallery._tag(moris_gallery._STATE.get("version"))
    thumbs = {}
    for example in examples:
        meta = moris_gallery._STATE.get("by_directory", {}).get(example["directory"])
        if meta and meta.get("thumb"):
            thumbs[example["leaf"]] = f"{prefix}{moris_gallery.ASSETS}/{meta['thumb']}{version}"
    return thumbs


def _graph_html(data, mode, height):
    payload = html.escape(json.dumps(data, separators=(",", ":")), quote=False)
    return (f'<div class="exa-graph exa-graph-{mode}" data-mode="{mode}" style="height:{height}px">'
            f'<script type="application/json">{payload}</script>'
            '<noscript>The graph needs JavaScript; the list below holds the same links.</noscript></div>')


def related_markdown(app, example, examples, docname):
    """The "Related examples" block of one example page (MyST markdown with raw HTML)."""
    near = neighbours(example, examples)
    if not near:
        return ""
    urls = _page_urls(app, examples, docname)
    data = graph_data(examples, urls, _thumb_urls(examples, docname), centre=example["leaf"])
    map_url = posixpath.relpath("examples/graph.html", posixpath.dirname(docname) or ".")
    # a sphinx-design dropdown (open), not a heading: the block must not enter the page's own outline
    lines = [":::{dropdown} Related examples", ":open:", ":icon: git-branch", "",
             "Nearest examples by the features detected in the decks; an edge is wider the more "
             "features two examples share. Click a node to open its page; the whole map is on "
             f"[the example map]({map_url}).",
             "", _graph_html(data, "ego", 320), ""]
    for other, score in near:
        title = f"[{other['title']}]({urls[other['leaf']]})" if other["leaf"] in urls else other["title"]
        shared = ", ".join(shared_tags(example, other)) or "same category"
        lines.append(f"- {title} — shares *{shared}*")
    lines.append(":::")
    return "\n".join(lines)


def full_graph_markdown(app, examples, docname):
    urls = _page_urls(app, examples, docname)
    data = graph_data(examples, urls, _thumb_urls(examples, docname), limit=MAP_NEIGHBOURS)
    return _graph_html(data, "map", 720)


def _directive():
    from sphinx.util.docutils import SphinxDirective

    class ExaGraphDirective(SphinxDirective):
        """```{exa-graph}``` — the relational graph of every example on one page."""

        has_content = False

        def run(self):
            from pathlib import Path

            import exa_catalog

            examples = exa_catalog.collect_examples(Path(self.config.moris_source_dir))
            return self.parse_text_to_nodes(full_graph_markdown(self.env.app, examples, self.env.docname))

    return ExaGraphDirective


def setup(app):
    app.add_directive("exa-graph", _directive())
    app.add_js_file(D3_URL)
    app.add_js_file("exa_graph.js")
