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
        candidates = [i for i in candidates if PurePath(i.path.replace("\\", "/")).name == override["file"]]
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


PHASE = re.compile(r"_(?P<cut>[cn])_p(?P<phase>\d+)$", re.IGNORECASE)


def block_phase(name: str) -> tuple[int | None, str]:
    """Phase index and a readable label for an XTK element block name (``HMR_dummy_c_p2``)."""
    match = PHASE.search(name)
    if not match:
        return None, name
    state = "cut" if match["cut"].lower() == "c" else "uncut"
    return int(match["phase"]), f"phase {match['phase']} ({state})"
