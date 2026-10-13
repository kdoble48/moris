"""The relational graph of the examples (doc/_ext/exa_graph.py): similarity, neighbours, graph data."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_ext"))

import exa_graph as eg  # noqa: E402


def ex(leaf, tags, category="optimization", group=""):
    return {"leaf": leaf, "title": leaf.replace("_", " "), "category": category, "group": group, "tags": tags}


BOXBEAM = ex("LevelSet_Boxbeam", ["2D", "linear elasticity", "GCMMA", "field array", "ghost", "Nitsche",
                                  "level-set design"])
BOXBEAM_AR = ex("LevelSet_Boxbeam_Adaptive_Refinement", ["2D", "linear elasticity", "GCMMA", "field array",
                                                         "ghost", "Nitsche", "level-set design", "adaptive refinement"])
ZZ = ex("Zienkiewicz_Zhu", ["2D", "linear elasticity", "diffusion", "GCMMA", "ghost", "Nitsche", "level-set design",
                            "heat method", "adaptive refinement", "von Mises stress"])
HEAT = ex("HeatConduction_2D", ["2D", "diffusion", "Nitsche", "ghost"], category="thermal", group="diffusion")
LAPLACE = ex("Laplace_2D", ["2D", "diffusion", "Nitsche", "eigenvalue"], category="thermal", group="diffusion")
MESH = ex("obj_sphere", ["3D", "surface mesh (OBJ)"], category="mesh_generation")
ALL = [BOXBEAM, BOXBEAM_AR, ZZ, HEAT, LAPLACE, MESH]


def test_similarity_orders_by_shared_features():
    assert eg.similarity(BOXBEAM, BOXBEAM_AR) > eg.similarity(BOXBEAM, ZZ) > eg.similarity(BOXBEAM, HEAT)
    assert eg.similarity(BOXBEAM, MESH) < eg.MIN_SIMILARITY            # nothing in common but a tag or two
    assert eg.similarity(HEAT, LAPLACE) > eg.similarity(HEAT, ZZ)        # same category + group bonus
    assert eg.similarity(BOXBEAM, BOXBEAM_AR) == eg.similarity(BOXBEAM_AR, BOXBEAM)
    assert eg.shared_tags(BOXBEAM, ZZ) == ["2D", "linear elasticity", "GCMMA", "ghost", "Nitsche", "level-set design"]


def test_neighbours_best_first_and_bounded():
    near = eg.neighbours(BOXBEAM, ALL, limit=2)
    assert [e["leaf"] for e, _ in near] == ["LevelSet_Boxbeam_Adaptive_Refinement", "Zienkiewicz_Zhu"]
    assert all(score >= eg.MIN_SIMILARITY for _, score in near)
    assert BOXBEAM not in [e for e, _ in eg.neighbours(BOXBEAM, ALL)]    # never itself
    assert eg.neighbours(MESH, ALL) == []                                 # an isolated example has no block


def test_graph_data_ego_and_map():
    ego = eg.graph_data(ALL, urls={"Zienkiewicz_Zhu": "Zienkiewicz_Zhu.html"}, centre="LevelSet_Boxbeam")
    ids = {n["id"] for n in ego["nodes"]}
    assert "LevelSet_Boxbeam" in ids and "obj_sphere" not in ids
    assert sum(n["centre"] for n in ego["nodes"]) == 1
    assert next(n for n in ego["nodes"] if n["id"] == "Zienkiewicz_Zhu")["url"] == "Zienkiewicz_Zhu.html"
    edge = next(e for e in ego["edges"] if {e["source"], e["target"]} == {"LevelSet_Boxbeam", "Zienkiewicz_Zhu"})
    assert "GCMMA" in edge["shared"] and 0 < edge["weight"] <= 1.3
    full = eg.graph_data(ALL)
    assert {n["id"] for n in full["nodes"]} == {e["leaf"] for e in ALL}
    keys = [(e["source"], e["target"]) for e in full["edges"]]
    assert len(keys) == len(set(keys))                                    # each pair once
    assert not any(n["centre"] for n in full["nodes"])
