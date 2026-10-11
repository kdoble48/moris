"""Packaging helpers of doc/gallery/build_gallery.py: file discovery, thumbnails, scene checks, manifest."""
import json
import sys
import zipfile
from pathlib import Path

import pytest
from PIL import Image

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gallery"))

import build_gallery as bg  # noqa: E402


def test_exodus_files_skip_parallel_pieces(tmp_path):
    for name in ["xtk_temp.exo", "xtk_temp.exo.4.0", "X.exo.e-s.0003", "X.exo.e-s.0003.4.1", "mesh.e",
                 "receipt.xml", "X.so"]:
        (tmp_path / name).write_text("")
    assert [p.name for p in bg.exodus_files(tmp_path)] == ["X.exo.e-s.0003", "mesh.e", "xtk_temp.exo"]


def test_run_inputs_and_showcase_job(tmp_path):
    source, build = tmp_path / "src", tmp_path / "bin"
    source.mkdir()
    build.mkdir()
    for name in ["Deck.cpp", "mesh.e", "Voxel_file.hdf5", "CMakeLists.txt"]:
        (source / name).write_text("")
    for name in ["Deck.so", "mesh.e", "Voxel_file.hdf5", "xtk_temp.exo", "Deck.exo.e-s.0003",
                 "ADV_Alg_0_Iter_3.hdf5", "Parameter_Receipt.xml", "stray.txt"]:
        (build / name).write_text("")
    assert [p.name for p in bg.run_inputs(build, source)] == ["Deck.so", "Voxel_file.hdf5", "mesh.e"]

    example = {"leaf": "Deck_Example", "name": "Deck_Example", "directory": "projects/EXA/x/Deck_Example"}
    job = bg.showcase_job(example, tmp_path / "build_opt", tmp_path / "EXA-test.exe", tmp_path / "runs")
    assert job.startswith("# Docs gallery: full showcase run of Deck_Example")
    assert "# mem: 2\n" in job
    assert "MORIS_EXA_MODE=showcase OMP_NUM_THREADS=1 nice -n 19" in job
    assert "'[EXA_Deck_Example]'" in job and job.rstrip().endswith("date > DONE")


def test_crop_thumbnail_centers_content(tmp_path):
    im = Image.new("RGB", (1600, 1000), "white")
    im.paste(Image.new("RGB", (200, 100), "black"), (100, 100))
    im.save(tmp_path / "raw.png")
    bg.crop_thumbnail(tmp_path / "raw.png", tmp_path / "thumb.png")
    thumb = Image.open(tmp_path / "thumb.png")
    assert thumb.size == (800, 500)
    assert thumb.getpixel((400, 250)) == (0, 0, 0)        # content centered
    assert thumb.getpixel((2, 2)) == (255, 255, 255)      # white margin


def test_make_animation_crops_and_labels(tmp_path):
    paths = []
    for k in range(3):
        im = Image.new("RGB", (1600, 1000), "white")
        im.paste(Image.new("RGB", (400, 200), (0, 0, 80 * k)), (600, 400))
        paths.append(tmp_path / f"f{k}.png")
        im.save(paths[-1])
    bg.make_animation(paths, tmp_path / "a.webp", labels=["it 1", "it 2", "it 3"])
    anim = Image.open(tmp_path / "a.webp")
    assert anim.n_frames == 3
    assert anim.width < 600 and anim.height < 300      # cropped to the content (+ label band)


def _scene(tmp_path, scene, hashes=True):
    path = tmp_path / "s.vtksz"
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("index.json", json.dumps({"scene": scene, "hashes": {"h": {}} if hashes else {}}))
    return path


def test_validate_scene(tmp_path):
    good = {"type": "vtkRenderWindow", "dependencies": [
        {"type": "vtkOpenGLActor", "dependencies": [{"type": "vtkPolyData"}]}]}
    assert "vtkPolyData" in bg.validate_scene(_scene(tmp_path, good))
    with pytest.raises(ValueError, match="no actors"):
        bg.validate_scene(_scene(tmp_path, {"type": "vtkRenderWindow"}))
    with pytest.raises(ValueError, match="indexed lookup"):
        bad = dict(good, dependencies=good["dependencies"] + [{"type": "x", "properties": {"indexedLookup": 1}}])
        bg.validate_scene(_scene(tmp_path, bad))
    with pytest.raises(ValueError, match="no geometry"):
        bg.validate_scene(_scene(tmp_path, good, hashes=False))


def test_manifest_order_and_pack(tmp_path):
    views = [{"id": "B", "category": "thermal", "title": "b", "png": "B.png", "thumb": "B-thumb.png",
              "scene": "B.vtksz"},
             {"id": "A", "category": "optimization", "title": "a", "png": "A.png", "thumb": "A-thumb.png",
              "scene": "A.vtksz"}]
    for v in views:
        (tmp_path / f"{v['id']}.json").write_text(json.dumps(v))
        for key in ("png", "thumb", "scene"):
            (tmp_path / v[key]).write_text("x")
    (tmp_path / "examples.json").write_text("[]")
    manifest = bg.collect_manifest(tmp_path, "abc", "build_opt @ def")
    assert [v["id"] for v in manifest["views"]] == ["A", "B"]
    names = bg.pack(tmp_path, manifest, tmp_path / "out.zip")
    assert names[-1] == "manifest.json" and "A.vtksz" in names and "examples.json" not in names
    assert set(zipfile.ZipFile(tmp_path / "out.zip").namelist()) == set(names)
