"""Sphinx extension that catalogs the example decks under projects/EXA.

At builder-inited it writes doc/examples/_generated/exa_catalog.md (gitignored),
which doc/examples/catalog.md pulls in with ``{include}``. It scans every example
directory registered with ``moris_add_example`` in its CMakeLists.txt and renders
one table per category (optimization, structure, thermal, fluid, mesh
generation). Each row is derived from the sources at build time — nothing is
hand-maintained, so the catalog cannot drift from the code:

- decks and ctest cases from the ``moris_add_example`` call (INPUTS, PROCS, CASES);
- the title and feature bullets from example_documentation.dox, when filled in;
- feature tags detected in the deck source (physics family from the IWG types,
  optimizer, geometry fields, ghost stabilization, spatial dimension, deck style).
"""

import re
from pathlib import Path

from sphinx.util import logging
from sphinx.util.docutils import SphinxDirective

LOGGER = logging.getLogger(__name__)

CATEGORY_TITLES = {
    "optimization": "Optimization and sensitivity analysis",
    "structure": "Structural mechanics",
    "thermal": "Thermal",
    "fluid": "Fluid dynamics",
    "mesh_generation": "Mesh generation",
}

# IWG type prefix -> physics tag (first match wins; order matters)
PHYSICS_TAGS = [
    ("STRUC_NONLINEAR_CONTACT", "contact"),
    ("STRUC_LINEAR_CONTACT", "contact"),
    ("STRUC_NON_LINEAR", "nonlinear elasticity"),
    ("STRUC_LINEAR", "linear elasticity"),
    ("SPATIALDIFF", "diffusion"),
    ("ADVECTION", "advection"),
    ("INCOMPRESSIBLE_NS", "incompressible flow"),
    ("COMPRESSIBLE_NS", "compressible flow"),
    ("SPALART_ALLMARAS", "turbulence"),
    ("TIME_CONTINUITY", "transient"),
    ("HJTEST", "Hamilton-Jacobi"),
    ("HELMHOLTZ", "Helmholtz filter"),
    ("FS_STRUC_INTERFACE", "fluid-structure"),
]

FIELD_TAGS = {
    "USER_DEFINED": "user-defined geometry",
    "SIGNED_DISTANCE_OBJECT": "surface mesh (OBJ)",
    "SIGNED_DISTANCE_IMAGE": "image geometry",
    "NODAL_FROM_FILE": "field from file",
    "VOXEL": "voxel geometry",
}


def _strip_cmake_comments(text):
    return "\n".join(line.split("#", 1)[0] for line in text.splitlines())


def _parse_add_example(cmake_text):
    """Return a list of dicts for each moris_add_example(...) call."""
    calls = []
    body_text = _strip_cmake_comments(cmake_text)
    # first non-empty set(VAR values...) per variable: the ungated default
    variables = {}
    for name, values in re.findall(r"set\s*\(\s*(\w+)\s+([^)]*)\)", body_text):
        if values.split() and name not in variables:
            variables[name] = values.split()
    for match in re.finditer(r"moris_add_example\s*\((.*?)\)", body_text, re.DOTALL):
        tokens = []
        for token in match.group(1).split():
            reference = re.fullmatch(r"\$\{(\w+)\}", token)
            tokens += variables.get(reference.group(1), []) if reference else [token]
        keywords = {"NAME", "TEST_BASE", "PROCS", "NO_PROCS_SUFFIX", "SOURCES", "INPUTS",
                    "EXTRA_SO_INCLUDES", "CASES"}
        args, current = {}, None
        for token in tokens:
            if token in keywords:
                current = token
                args.setdefault(current, [])
            elif current:
                args[current].append(token)
        calls.append(args)
    return calls


def _dox_info(dox_path):
    """Title and real (non-TODO) feature bullets from example_documentation.dox."""
    if not dox_path.exists():
        return None, []
    text = dox_path.read_text(errors="replace")
    title_match = re.search(r"\\page\s+\S+\s+(.+)", text)
    title = title_match.group(1).strip().rstrip(".") if title_match else None
    features = []
    description = re.search(r"# Problem Description(.*?)(?:\n#\s|\*/)", text, re.DOTALL)
    if description:
        for line in description.group(1).splitlines():
            line = line.strip()
            if line.startswith("- ") and "TODO" not in line:
                features.append(line[2:].strip())
    return title, features


def _deck_tags(deck_text):
    tags = []
    iwg_types = set(re.findall(r"IWG_Type::(\w+)", deck_text))
    iwg_types |= set(re.findall(r'"IWG_type"\s*,\s*"?(\w+)', deck_text))
    for prefix, tag in PHYSICS_TAGS:
        if any(iwg.upper().startswith(prefix) for iwg in iwg_types) and tag not in tags:
            tags.append(tag)

    algorithms = set(re.findall(r"Optimization_Algorithm_Type::(\w+)", deck_text))
    algorithms |= {name.upper() for name in re.findall(r"create_(gcmma|sweep|lbfgs|sqp|mma|rol)_parameter_list", deck_text)}
    for algorithm in sorted(algorithms):
        tags.append("sensitivity sweep" if algorithm == "SWEEP" else algorithm)

    for field, tag in FIELD_TAGS.items():
        if re.search(r"(Field_Type|Geometry_Type)::" + field + r"\b", deck_text):
            tags.append(tag)
    if "create_field_array_parameter_list" in deck_text:
        tags.append("field array")
    if re.search(r'"ghost_stab"\s*,\s*(true|tUseGhost|sUseGhost|1)', deck_text):
        tags.append("ghost")
    if "NITSCHE" in deck_text:
        tags.append("Nitsche")
    if "STK_XTK" in deck_text or "STK_FEM" in deck_text:
        tags.append("STK mesh")
    if "MORIS_DECK" in deck_text:
        tags.append("Deck API v2")

    dimension = _dimension(deck_text)
    if dimension:
        tags.insert(0, dimension)
    return tags


def _dimension(deck_text):
    """2D/3D from the HMR element counts, when written literally."""
    match = re.search(r'"number_of_elements_per_dimension"\s*,\s*([^)]*)\)', deck_text)
    if not match:
        return None
    value = match.group(1)
    if "{" in value or "tNum" in value:
        return None
    count = len([part for part in re.split(r"[,\s\"]+", value) if part])
    return {2: "2D", 3: "3D"}.get(count)


def _github_url(app, relative_path, directory=False):
    repository = app.config.html_theme_options.get("source_repository", "").rstrip("/")
    branch = app.config.html_theme_options.get("source_branch", "develop")
    kind = "tree" if directory else "blob"
    return f"{repository}/{kind}/{branch}/{relative_path}" if repository else None


def _built_directories(exa_dir):
    """Directories reachable from projects/EXA/CMakeLists.txt through uncommented add_subdirectory()."""
    reached, pending = set(), [exa_dir]
    while pending:
        directory = pending.pop()
        if directory in reached:
            continue
        reached.add(directory)
        cmake_file = directory / "CMakeLists.txt"
        if not cmake_file.exists():
            continue
        text = _strip_cmake_comments(cmake_file.read_text(errors="replace"))
        for name in re.findall(r"add_subdirectory\s*\(\s*([^\s)]+)", text):
            pending.append((directory / name).resolve())
    return reached


def collect_examples(source_dir):
    exa_dir = source_dir / "projects" / "EXA"
    built = _built_directories(exa_dir.resolve())
    examples = []
    for cmake_file in sorted(exa_dir.rglob("CMakeLists.txt")):
        calls = _parse_add_example(cmake_file.read_text(errors="replace"))
        if not calls:
            continue
        directory = cmake_file.parent
        relative = directory.relative_to(exa_dir)
        title, features = _dox_info(directory / "example_documentation.dox")
        for call in calls:
            test_base = (call.get("TEST_BASE") or call.get("NAME") or [directory.name])[0]
            inputs = call.get("INPUTS") or [test_base]
            decks = [directory / f"{name}.cpp" for name in inputs]
            deck_text = "\n".join(deck.read_text(errors="replace") for deck in decks if deck.exists())
            examples.append({
                "category": relative.parts[0],
                "group": "/".join(relative.parts[1:-1]),
                "name": (call.get("NAME") or [directory.name])[0],
                "title": title if len(calls) == 1 and title else directory.name.replace("_", " "),
                "directory": directory.relative_to(source_dir).as_posix(),
                "decks": [deck.relative_to(source_dir).as_posix() for deck in decks if deck.exists()],
                "leaf": directory.name,
                "built": directory.resolve() in built,
                "test_base": test_base,
                "no_procs_suffix": "NO_PROCS_SUFFIX" in call,
                "procs": call.get("PROCS", []),
                "cases": call.get("CASES", []),
                "features": features,
                "tags": _deck_tags(deck_text),
            })
    return examples


def _cell(text):
    return text.replace("|", "\\|")


def page_docname(example):
    """Docname of the example's own page: examples/exa/<category>/<leaf directory>."""
    return f"examples/exa/{example['category']}/{example['leaf']}"


def ctest_names(example):
    """The ctest names moris_add_example registers for this example."""
    bases = ([f"{example['test_base']}_{case}" for case in example["cases"]]
             if example["cases"] else [example["test_base"]])
    if example["no_procs_suffix"]:
        return bases
    return [f"{base}-{procs}-procs" for base in bases for procs in example["procs"]]


def summary_markdown(app, example):
    """The generated facts block shown at the top of each example page."""
    tags = " ".join(f'<span class="exa-tag">{tag}</span>' for tag in example["tags"])
    decks = ", ".join(f"[{Path(deck).name}]({_github_url(app, deck)})" for deck in example["decks"])
    tests = ", ".join(f"`{name}`" for name in ctest_names(example))
    if not example["built"]:
        tests = ("**not built** — its directory is commented out of the parent `CMakeLists.txt`, "
                 f"so these tests are not registered: {tests}")
    source = _github_url(app, example["directory"], directory=True)
    rows = [
        ("Source", f"[`{example['directory']}`]({source})"),
        ("Decks", decks or "—"),
        ("ctest", tests or "—"),
        ("Detected features", tags or "—"),
    ]
    lines = ["| | |", "|---|---|"] + [f"| **{key}** | {value} |" for key, value in rows]
    return "\n".join(lines)


def render_markdown(app, examples):
    lines = [
        f"{len(examples)} example decks, generated from `projects/EXA` at build time. "
        "Tags are detected in the deck source; each example is also a ctest integration "
        "test (run with `MORIS_USE_EXAMPLES=ON`).",
        "",
    ]
    categories = sorted({example["category"] for example in examples},
                        key=lambda name: list(CATEGORY_TITLES).index(name) if name in CATEGORY_TITLES else 99)
    for category in categories:
        members = [example for example in examples if example["category"] == category]
        lines += [f"## {CATEGORY_TITLES.get(category, category)}", "",
                  "| Example | Features | Decks | Tests |", "|---|---|---|---|"]
        for example in sorted(members, key=lambda item: (item["group"], item["title"].lower())):
            directory_url = _github_url(app, example["directory"], directory=True)
            if (Path(app.srcdir) / f"{page_docname(example)}.md").exists():
                title = f"{{doc}}`{_cell(example['title'])} </{page_docname(example)}>`"
            elif directory_url:
                title = f"[{_cell(example['title'])}]({directory_url})"
            else:
                title = _cell(example["title"])
            if example["group"]:
                title += f"<br><small>{example['group']}</small>"
            features = " ".join(f'<span class="exa-tag">{tag}</span>' for tag in example["tags"])
            if example["features"]:
                features += "<br><small>" + "; ".join(_cell(feature) for feature in example["features"]) + "</small>"
            decks = "<br>".join(
                f"[{Path(deck).name}]({_github_url(app, deck)})" for deck in example["decks"])
            tests = ", ".join(f"{procs} proc{'s' if procs != '1' else ''}" for procs in example["procs"])
            if example["cases"]:
                tests += f"<br>cases {', '.join(example['cases'])}"
            if not example["built"]:
                tests = "not built"
            lines.append(f"| {title} | {features or '—'} | {decks or '—'} | {tests or '—'} |")
        lines.append("")
    return "\n".join(lines)


GENERATED_FILE = Path("examples") / "_generated" / "exa_catalog.md"


def write_catalog(app):
    """builder-inited: regenerate doc/examples/_generated/exa_catalog.md."""
    source_dir = Path(app.config.moris_source_dir)
    examples = collect_examples(source_dir)
    if not examples:
        LOGGER.warning("exa_catalog: no examples found under %s/projects/EXA", source_dir)
    target = Path(app.srcdir) / GENERATED_FILE
    target.parent.mkdir(parents=True, exist_ok=True)
    text = render_markdown(app, examples)
    # rewrite only on change, so incremental builds stay incremental
    if not target.exists() or target.read_text() != text:
        target.write_text(text)


class ExaSummaryDirective(SphinxDirective):
    """```{exa-summary} projects/EXA/<category>/.../<Example>``` — generated facts table."""

    required_arguments = 1
    has_content = False

    def run(self):
        source_dir = Path(self.config.moris_source_dir)
        directory = self.arguments[0].strip().rstrip("/")
        matches = [example for example in collect_examples(source_dir)
                   if example["directory"] == directory]
        if not matches:
            LOGGER.warning("exa-summary: no moris_add_example registered in %s", directory,
                           location=self.get_location())
            return []
        self.env.note_dependency(str(source_dir / directory / "CMakeLists.txt"))
        import moris_gallery  # scenes from the gallery manifest, when one was downloaded or built

        scene = moris_gallery.scene_for(directory, self.env.docname)
        return self.parse_text_to_nodes(scene + "\n" + summary_markdown(self.env.app, matches[0]))


def setup(app):
    app.add_config_value("moris_source_dir", "", "env")
    app.connect("builder-inited", write_catalog)
    app.add_directive("exa-summary", ExaSummaryDirective)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
