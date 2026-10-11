"""Rules that pick the file and field each example's scene shows (doc/gallery/gallery_select.py)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gallery"))

import gallery_select as gs  # noqa: E402

BIN = "/b/projects/EXA/optimization/X/bin/"


def info(name, arrays=(), points=100):
    return gs.FileInfo(BIN + name, list(arrays), points, 1)


def test_result_beats_xtk_and_gen():
    choice = gs.choose([info("xtk_temp.exo"), info("GEN_X.exo", ["Field 0"]),
                        info("X.exo", ["Node_Id", "U", "TEMP"])])
    assert (Path(choice.path).name, choice.kind, choice.array, choice.label, choice.vector) == \
        ("X.exo", "result", "U", "|U|", True)


def test_xtk_beats_gen():
    choice = gs.choose([info("GEN_X.exo", ["Field 0"]), info("xtk_temp.exo", ["ids"])])
    assert choice.kind == "xtk" and choice.array is None


def test_field_preference_order():
    assert gs.pick_array(["P", "TEMP", "V"])[0] == "TEMP"
    assert gs.pick_array(["Node_Id", "Node_Index", "Field 3"]) == ("Field 3", "Field 3", False)
    assert gs.pick_array(["Node_Id"]) == (None, None, False)


def test_last_file_of_series_wins():
    paths = [BIN + f"X.exo.e-s.{i:04d}" for i in (1, 12, 3)] + [BIN + "xtk_temp.exo"]
    assert sorted(Path(p).name for p in gs.last_of_series(paths)) == ["X.exo.e-s.0012", "xtk_temp.exo"]
    choice = gs.choose([info(Path(p).name, ["U"]) for p in paths[:3]])
    assert Path(choice.path).name == "X.exo.e-s.0012"


def test_kinds():
    assert gs.kind_of(info("stk_fem_mesh.exo", ["GlobalNodeId", "PedigreeNodeId"])) == "mesh"
    assert gs.kind_of(info("stk_fem_mesh.exo", [])) == "mesh"
    assert gs.kind_of(info("other.exo", ["Damage"])) == "result"  # an unknown field still counts
    assert gs.kind_of(info("Laplace_GEN.exo.e-s.0000", ["Field 2"])) == "gen"
    assert gs.kind_of(info("xtk_temp_ip.exo")) == "xtk"


def test_empty_files_and_override():
    assert gs.choose([info("X.exo", ["U"], points=0)]) is None
    choice = gs.choose([info("X.exo", ["U", "TEMP"]), info("xtk_temp.exo")],
                       {"file": "X.exo", "array": "TEMP"})
    assert (choice.array, choice.label, choice.vector) == ("TEMP", "TEMP", False)
    choice = gs.choose([info("X.exo", ["U"]), info("xtk_temp.exo")], {"file": "xtk_temp.exo"})
    assert choice.kind == "xtk"


def test_block_phase_labels():
    assert gs.block_phase("hmr_dummy_c_p2") == (2, "phase 2 (cut)")
    assert gs.block_phase("HMR_dummy_n_p10") == (10, "phase 10 (uncut)")
    assert gs.block_phase("block_main") == (None, "block_main")
