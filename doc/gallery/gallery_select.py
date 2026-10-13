"""Choose what to render for an example from the Exodus files its ctest run left behind.

Pure Python (no ParaView import), so the rules are unit-tested: ``doc/tests/test_gallery_select.py``.
``export_scenes.py`` gathers a :class:`FileInfo` per candidate file under pvpython and calls
:func:`choose`.

Preference, highest first:

1. ``result`` — a VIS output with a physics field (displacement, temperature, velocity, …): the
   last file of an ``.e-s.NNNN`` series (last optimization iteration), last time step.
2. ``xtk`` — the XTK integration mesh (``xtk_temp*.exo``), colored by element block (phase): the
   cut geometry, which says more about a level-set design than a raw GEN field.
3. ``gen`` — a GEN output with level-set / geometry fields (``GEN_*`` or ``*_GEN*`` names).
4. ``mesh`` — any other mesh with no fields (input or mesh-generation meshes), drawn plain.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import PurePath

# Physics fields in order of preference; each entry: (array names that match, label, vector?).
PHYSICS_FIELDS = [
    (("U",), "|U|", True),
    (("TEMP", "Temp"), "TEMP", False),
    (("V",), "|V|", True),
    (("P",), "P", False),
    (("LVLSET", "LEVELSET", "LS", "Levelset", "LevelsetField", "Level_Set_Field"), None, False),
    (("VOLUME", "VOL"), None, False),
    (("STRESS",), None, False),
    (("HEAT",), None, False),
]
# Arrays that are bookkeeping, never worth coloring by.
IGNORED_ARRAYS = {"ids", "Node_Id", "Node_Index", "cell_id", "object_id", "ClusterVolume", "PHID", "THETA",
                  # added by ParaView's legacy Exodus reader
                  "GlobalNodeId", "PedigreeNodeId", "GlobalElementId", "PedigreeElementId", "ObjectId"}
SERIES = re.compile(r"^(?P<stem>.+)\.e-s\.(?P<index>\d+)$")
KIND_RANK = {"result": 0, "xtk": 1, "gen": 2, "mesh": 3}


@dataclass
class FileInfo:
    path: str
    point_arrays: list[str] = field(default_factory=list)
    points: int = 0
    times: int = 0


@dataclass
class Choice:
    path: str
    kind: str              # result | gen | xtk | mesh
    array: str | None      # array to color by (None: blocks for xtk, plain for mesh)
    label: str | None      # legend label (e.g. "|U|")
    vector: bool = False


def kind_of(info: FileInfo) -> str:
    name = PurePath(info.path.replace("\\", "/")).name
    if name.startswith("xtk_temp"):
        return "xtk"
    useful = [a for a in info.point_arrays if a not in IGNORED_ARRAYS]
    if not useful:
        return "mesh"
    if name.startswith("GEN_") or "_GEN" in name or "_Gen" in name or name.startswith("gen_"):
        return "gen"
    return "result"


def pick_array(arrays: list[str]) -> tuple[str | None, str | None, bool]:
    """``(array, label, vector)`` for the preferred field in ``arrays``."""
    for names, label, vector in PHYSICS_FIELDS:
        for name in names:
            if name in arrays:
                return name, label or name, vector
    useful = [a for a in arrays if a not in IGNORED_ARRAYS]
    return (useful[0], useful[0], False) if useful else (None, None, False)


def last_of_series(paths: list[str]) -> list[str]:
    """Keep only the highest ``.e-s.NNNN`` file of each series (the last iteration)."""
    best: dict[str, tuple[int, str]] = {}
    plain = []
    for path in paths:
        match = SERIES.match(PurePath(path.replace("\\", "/")).name)
        if not match:
            plain.append(path)
            continue
        key = str(PurePath(path.replace("\\", "/")).parent / match["stem"])
        index = int(match["index"])
        if key not in best or index > best[key][0]:
            best[key] = (index, path)
    return plain + [path for _, path in best.values()]


def choose(infos: list[FileInfo], override: dict | None = None) -> Choice | None:
    """The file and array to render, or None when no candidate has any geometry."""
    override = override or {}
    candidates = [i for i in infos if i.points > 0]
    if override.get("file"):
        # a file name or the stem of an .e-s series (e.g. "Deck.exo" selects Deck.exo.e-s.0101)
        candidates = [i for i in candidates if PurePath(i.path.replace("\\", "/")).name == override["file"]
                      or PurePath(i.path.replace("\\", "/")).name.startswith(override["file"] + ".e-s.")]
    keep = set(last_of_series([i.path for i in candidates]))
    candidates = [i for i in candidates if i.path in keep]
    if not candidates:
        return None

    def rank(info: FileInfo):
        kind = kind_of(info)
        array = pick_array(info.point_arrays)[0] if kind in ("result", "gen") else None
        field_rank = next((n for n, (names, _, _) in enumerate(PHYSICS_FIELDS) if array in names), 99)
        return (KIND_RANK[kind], field_rank, -info.points, info.path)

    best = min(candidates, key=rank)
    kind = kind_of(best)
    if kind in ("result", "gen"):
        array, label, vector = pick_array(best.point_arrays)
        if override.get("array") in best.point_arrays:
            array = override["array"]
            label = override.get("label", array)
            vector = bool(override.get("vector", False))
        return Choice(best.path, kind, array, label, vector)
    return Choice(best.path, kind, None, None)


LEVEL_SET_ARRAYS = ("LVLSET", "LEVELSET", "LS", "Levelset", "LevelsetField", "Level_Set_Field")


def interface_array(arrays: list[str]) -> str | None:
    """The level-set array whose threshold contour is the material interface, if the file has one."""
    return next((name for name in LEVEL_SET_ARRAYS if name in arrays), None)


ISOCONTOUR = re.compile(r'"isocontour_threshold"\s*,\s*([-+0-9.eE]+)')


def isocontour_threshold(deck_text: str) -> float:
    """GEN's ``isocontour_threshold`` set in a deck (the level-set value of the interface); 0 by default."""
    match = ISOCONTOUR.search(deck_text)
    try:
        return float(match.group(1)) if match else 0.0
    except ValueError:  # a named constant: keep the default rather than guess
        return 0.0


MAX_FRAMES = 60


def frame_unit(path: str) -> str:
    """An ``.e-s.NNNN`` series holds one optimization iteration per file; otherwise steps are time steps."""
    return "iteration" if SERIES.match(PurePath(path.replace("\\", "/")).name) else "time step"


def series_numbers(path: str, all_paths: list[str]) -> list[int]:
    """Sorted ``NNNN`` indices of the ``.e-s.NNNN`` series ``path`` belongs to (empty if not a series)."""
    match = SERIES.match(PurePath(path.replace("\\", "/")).name)
    if not match:
        return []
    numbers = set()
    for other in all_paths:
        m = SERIES.match(PurePath(other.replace("\\", "/")).name)
        if m and m["stem"] == match["stem"]:
            numbers.add(int(m["index"]))
    return sorted(numbers)


def step_labels(path: str, all_paths: list[str], steps: int) -> tuple[list[str], int | None]:
    """A label per reader time step, and the number of optimization iterations (None for a transient run).

    The IOSS reader opens a whole ``.e-s.NNNN`` series as one time series: one file per optimization
    iteration, each holding one or more time steps. Labels use the real file numbers, so a restarted
    run reads ``iteration 11`` rather than ``1``."""
    numbers = series_numbers(path, all_paths)
    if not numbers or steps % len(numbers):
        return [f"time step {i + 1} of {steps}" for i in range(steps)], None
    per = steps // len(numbers)
    labels = []
    for i in range(steps):
        label = f"iteration {numbers[i // per]} of {numbers[-1]}"
        if per > 1:
            label += f", time step {i % per + 1} of {per}"
        labels.append(label)
    return labels, len(numbers)


def sample_indices(count: int, limit: int = MAX_FRAMES) -> list[int]:
    """Up to ``limit`` evenly spaced indices into ``count`` steps, always including the first and last."""
    if count <= limit:
        return list(range(count))
    return sorted({round(i * (count - 1) / (limit - 1)) for i in range(limit)})


def series_paths(path: str, all_paths: list[str]) -> list[str]:
    """The files of the ``.e-s.NNNN`` series ``path`` belongs to, in iteration order ([path] if not a series).

    Each file is read on its own: XTK renumbers element blocks between iterations (a phase that
    appears later takes over another phase's block id), so a reader that opens the series as one
    time series takes the first file's block layout and mislabels the blocks of later files."""
    match = SERIES.match(PurePath(path.replace("\\", "/")).name)
    if not match:
        return [path]
    parent = PurePath(path.replace("\\", "/")).parent
    found = {}
    for other in all_paths:
        p = PurePath(other.replace("\\", "/"))
        m = SERIES.match(p.name)
        if m and m["stem"] == match["stem"] and p.parent == parent:
            found[int(m["index"])] = other
    return [found[k] for k in sorted(found)] or [path]


PHASE = re.compile(r"_(?P<cut>[cn])_p(?P<phase>\d+)$", re.IGNORECASE)
SET_PHASE = re.compile(r"HMR_dummy_[cn]_p(\d+)", re.IGNORECASE)
STRING_DECL = re.compile(r"std::string\s+(\w+)\s*=\s*([^;]+);")
STRING_LITERAL = re.compile(r'^(?:std::string\s*\(\s*)?"([^"]*)"\s*\)?$')
PARAMETER = re.compile(r'"(IWG_name|IWG_type|dof_residual|mesh_set_names)"\s*,\s*(.+?)\s*\)\s*;', re.DOTALL)


def _string_value(expr: str, variables: dict[str, str]) -> str:
    """Value of a C++ string expression built from literals, ``std::string( "..." )`` and string variables
    joined with ``+``; unknown parts (function calls, other types) evaluate to "". A conditional
    ``c ? a : b`` gives both branches, comma-joined (e.g. the 2D and the 3D DOF list)."""
    ternary = re.fullmatch(r'[^?"]*\?\s*(.+?)\s*:\s*([^:]+)', expr.strip())
    if ternary:
        return _string_value(ternary.group(1), variables) + "," + _string_value(ternary.group(2), variables)
    parts, depth, current, quoted = [], 0, "", False
    for ch in expr:
        if ch == '"':
            quoted = not quoted
        elif not quoted and ch in "()":
            depth += 1 if ch == "(" else -1
        if ch == "+" and not quoted and depth == 0:
            parts.append(current)
            current = ""
        else:
            current += ch
    parts.append(current)
    value = ""
    for part in (p.strip() for p in parts):
        literal = STRING_LITERAL.match(part)
        if literal:
            value += literal.group(1)
        elif re.fullmatch(r"\w+", part):
            value += variables.get(part, "")
    return value


def bulk_phases(deck_texts: list[str], array: str | None, vector: bool) -> set[int] | None:
    """Phases with physics for ``array``: those whose sets a bulk IWG solving for that field integrates over.

    Void phases often carry auxiliary bulk IWGs (heat-method or level-set regularization fields on
    THETA/PHID) or none at all, and the solution field is nonzero there too (ersatz material), so the
    residual DOF is the test: |U| needs a bulk IWG with residual UX/UY(/UZ), TEMP one with residual TEMP.
    None when the decks give no answer (e.g. Deck API v2 helpers), meaning: show every phase."""
    if not array:
        return None

    def matches(dof: str) -> bool:
        return dof == array or (vector and len(dof) == len(array) + 1 and dof.startswith(array) and dof[-1] in "XYZ")

    phases: set[int] = set()
    for text in deck_texts:
        variables: dict[str, str] = {}
        for name, expr in STRING_DECL.findall(text):
            variables[name] = _string_value(expr, variables)
        # one IWG per "IWG_name": its type, residual and sets follow before the next one
        for chunk in re.split(r'(?="IWG_name")', text)[1:]:
            chunk = re.split(r'"IQI_name"|add_parameter_list', chunk)[0]
            values = {key: expr for key, expr in PARAMETER.findall(chunk)}
            if "BULK" not in values.get("IWG_type", "") or "mesh_set_names" not in values:
                continue
            residual = _string_value(values.get("dof_residual", ""), variables)
            if not any(matches(d.strip()) for d in re.split(r"[,;]", residual)):
                continue
            sets = _string_value(values["mesh_set_names"], variables)
            phases.update(int(p) for p in SET_PHASE.findall(sets))
    return phases or None


def block_phase(name: str) -> tuple[int | None, str]:
    """Phase index and a readable label for an XTK element block name (``HMR_dummy_c_p2``)."""
    match = PHASE.search(name)
    if not match:
        return None, name
    state = "cut" if match["cut"].lower() == "c" else "uncut"
    return int(match["phase"]), f"phase {match['phase']} ({state})"
