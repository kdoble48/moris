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


TAG_SHARE = 0.65          # the tag overlap carries most of the score ...
ANATOMY_SHARE = 0.35      # ... the deck anatomy (DOFs, solvers, size, block structure) the rest
BLOCK_KINDS = ("properties", "constitutive models", "stabilization parameters", "IWGs", "IQIs", "geometries")


def tag_similarity(a, b):
    """Weighted Jaccard index of the two examples' tags (0..1)."""
    tags_a, tags_b = set(a.get("tags", [])), set(b.get("tags", []))
    union = tags_a | tags_b
    if not union:
        return 0.0
    return sum(tag_weight(t) for t in tags_a & tags_b) / sum(tag_weight(t) for t in union)


def _jaccard(x, y):
    x, y = set(x), set(y)
    return len(x & y) / len(x | y) if x | y else None


def _stack(anatomy):
    solvers = anatomy.get("solvers", {})
    return [f"{s} (linear)" for s in solvers.get("linear", [])] + \
           [f"{s} (nonlinear)" for s in solvers.get("nonlinear", [])]


def anatomy_similarity(a, b):
    """``(score 0..1, what they share)`` from the deck anatomies: the DOF set, the solver stack, the
    mesh size (log scale) and the block structure (counts per kind, cosine). A part counts only
    when both decks define it; ``(None, [])`` when nothing is comparable."""
    import math

    aa, ab = a.get("anatomy") or {}, b.get("anatomy") or {}
    parts, shared = [], []
    dofs = _jaccard(aa.get("dofs", []), ab.get("dofs", []))
    if dofs is not None:
        parts.append((0.35, dofs))
        if dofs == 1.0:
            shared.append("same DOFs: " + ", ".join(aa["dofs"]))
        elif dofs > 0:
            shared.append("shared DOFs: " + ", ".join(d for d in aa["dofs"] if d in ab["dofs"]))
    solvers = _jaccard(_stack(aa), _stack(ab))
    if solvers is not None:
        parts.append((0.25, solvers))
        if solvers == 1.0:
            shared.append("same solver stack: " + ", ".join(s.split(" ")[0] for s in _stack(aa)))
    cells_a, cells_b = aa.get("mesh", {}).get("cells"), ab.get("mesh", {}).get("cells")
    if cells_a and cells_b:
        proximity = max(0.0, 1.0 - abs(math.log10(cells_a) - math.log10(cells_b)) / 1.5)
        parts.append((0.2, proximity))
        if proximity > 0.7:
            shared.append(f"similar mesh size (~{cells_a:,} and ~{cells_b:,} background cells)")
    counts_a = [aa.get("blocks", {}).get(k, 0) for k in BLOCK_KINDS]
    counts_b = [ab.get("blocks", {}).get(k, 0) for k in BLOCK_KINDS]
    norm = math.sqrt(sum(x * x for x in counts_a)) * math.sqrt(sum(x * x for x in counts_b))
    if norm:
        cosine = sum(x * y for x, y in zip(counts_a, counts_b)) / norm
        parts.append((0.2, cosine))
        if cosine > 0.95 and abs(sum(counts_a) - sum(counts_b)) <= 0.3 * max(sum(counts_a), sum(counts_b)):
            shared.append("similar deck structure (block counts)")
    if not parts:
        return None, []
    return sum(w * s for w, s in parts) / sum(w for w, _ in parts), shared


def similarity(a, b):
    """Tag overlap and deck anatomy, weighted, plus the category and group bonuses."""
    anatomy, _ = anatomy_similarity(a, b)
    score = tag_similarity(a, b) if anatomy is None else TAG_SHARE * tag_similarity(a, b) + ANATOMY_SHARE * anatomy
    if a.get("category") == b.get("category"):
        score += SAME_CATEGORY
        if a.get("group") and a.get("group") == b.get("group"):
            score += SAME_GROUP
    return round(score, 4)


def shared_tags(a, b):
    """The tags both examples carry, in the first example's order."""
    other = set(b.get("tags", []))
    return [t for t in a.get("tags", []) if t in other]


def shared_features(a, b):
    """What an edge stands for: the shared tags, then the anatomy matches."""
    return shared_tags(a, b) + anatomy_similarity(a, b)[1]


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
                                   "shared": shared_features(by_leaf[key[0]], by_leaf[key[1]])})
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
        shared = "; ".join(shared_features(example, other)) or "same category"
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
