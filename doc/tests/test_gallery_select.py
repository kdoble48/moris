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
    series = [info(f"A.exo.e-s.{i:04d}", ["U"], points=500) for i in (1, 2, 101)] \
        + [info("A_Quadratic.exo.e-s.0020", ["U"], points=900)]
    choice = gs.choose(series, {"file": "A.exo"})                   # series stem selects the series
    assert Path(choice.path).name == "A.exo.e-s.0101"


def test_interface_array_and_threshold():
    assert gs.interface_array(["U", "LVLSET", "PHID"]) == "LVLSET"
    assert gs.interface_array(["U", "TEMP"]) is None
    assert gs.isocontour_threshold('aParameterLists.set( "isocontour_threshold", 0.5 );') == 0.5
    assert gs.isocontour_threshold('set( "isocontour_threshold", tThreshold );') == 0.0
    assert gs.isocontour_threshold("no GEN threshold here") == 0.0


def test_step_labels():
    restart = [BIN + f"R.exo.e-s.{i:04d}" for i in range(11, 17)] + [BIN + "GEN_R.exo.e-s.0011"]
    labels, iterations = gs.step_labels(restart[0], restart, 6)
    assert iterations == 6 and labels[0] == "iteration 11 of 16" and labels[-1] == "iteration 16 of 16"
    mach = [BIN + "M.exo.e-s.0001", BIN + "M.exo.e-s.0002"]
    labels, iterations = gs.step_labels(mach[0], mach, 6)
    assert iterations == 2
    assert labels[:4] == ["iteration 1 of 2, time step 1 of 3", "iteration 1 of 2, time step 2 of 3",
                          "iteration 1 of 2, time step 3 of 3", "iteration 2 of 2, time step 1 of 3"]
    labels, iterations = gs.step_labels(BIN + "T.exo", [], 20)
    assert iterations is None and labels[-1] == "time step 20 of 20"


def test_animation_frames():
    assert gs.frame_unit(BIN + "X.exo.e-s.0007") == "iteration"
    assert gs.frame_unit(BIN + "X.exo") == "time step"
    assert gs.sample_indices(5) == [0, 1, 2, 3, 4]
    picked = gs.sample_indices(151, 60)
    assert len(picked) == 60 and picked[0] == 0 and picked[-1] == 150


def test_block_phase_labels():
    assert gs.block_phase("hmr_dummy_c_p2") == (2, "phase 2 (cut)")
    assert gs.block_phase("HMR_dummy_n_p10") == (10, "phase 10 (uncut)")
    assert gs.block_phase("block_main") == (None, "block_main")


ZZ_LIKE_DECK = """
    std::string tBulkSets  = "HMR_dummy_n_p9,HMR_dummy_c_p9";
    std::string tVoidSets  = "HMR_dummy_n_p8,HMR_dummy_c_p8";
    std::string tFrameSets = "HMR_dummy_n_p11,HMR_dummy_c_p11";
    std::string tTotalDomain = tBulkSets + "," + tVoidSets;
    std::string tDofStrg = tIs3D ? "UX,UY,UZ" : "UX,UY";
        aParameterLists.set( "IWG_name", "IWGDiffusionThetaBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "THETA" );
        aParameterLists.set( "mesh_set_names", tTotalDomain );
        aParameterLists.set( "IWG_name", "IWGBulkU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", tDofStrg );
        aParameterLists.set( "mesh_set_names", tBulkSets );
        aParameterLists.set( "IWG_name", std::string( "IWGBulkUFrame" ) );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", std::string( "UX,UY" ) );
        aParameterLists.set( "mesh_set_names", tFrameSets );
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "mesh_set_names", "iside_b0_11_b1_15" );
"""


def test_bulk_phases_from_deck():
    # the void phase p8 only carries an auxiliary (THETA) bulk IWG: not shown for |U|
    assert gs.bulk_phases([ZZ_LIKE_DECK], "U", True) == {9, 11}
    # for the auxiliary field itself, the void counts
    assert gs.bulk_phases([ZZ_LIKE_DECK], "THETA", False) == {8, 9}
    # nothing to go on: None (every phase is shown)
    assert gs.bulk_phases([ZZ_LIKE_DECK], "TEMP", False) is None
    assert gs.bulk_phases(["tFem.set( \"dof_residual\", deck::Dofs::Displacement2D );"], "U", True) is None
    assert gs.bulk_phases([ZZ_LIKE_DECK], None, False) is None


def test_series_paths_in_iteration_order():
    files = [BIN + "D.exo.e-s.0010", BIN + "D.exo.e-s.0002", BIN + "xtk_temp.exo", BIN + "D.exo.e-s.0001",
             "/other/bin/D.exo.e-s.0003"]
    assert gs.series_paths(BIN + "D.exo.e-s.0010", files) == [BIN + "D.exo.e-s.0001", BIN + "D.exo.e-s.0002",
                                                               BIN + "D.exo.e-s.0010"]
    assert gs.series_paths(BIN + "xtk_temp.exo", files) == [BIN + "xtk_temp.exo"]
