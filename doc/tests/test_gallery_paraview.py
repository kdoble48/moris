"""Per-example ParaView bundles: zipped as their own release assets, linked from the page."""
import json
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gallery"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_ext"))

import build_gallery as bg  # noqa: E402
import moris_gallery as mg  # noqa: E402


def test_bundle_paraview_zips_each_folder(tmp_path):
    folder = tmp_path / "Deck_Example-paraview"
    folder.mkdir()
    (folder / "Deck.exo").write_bytes(b"x" * 1000)
    (folder / "state.pvsm").write_text("<ParaView/>")
    (folder / "open.py").write_text("pass\n")
    (tmp_path / "Deck_Example.png").write_bytes(b"png")           # not a bundle folder: untouched
    bundles = bg.bundle_paraview(tmp_path)
    assert list(bundles) == ["Deck_Example"]
    assert bundles["Deck_Example"]["file"] == "Deck_Example-paraview.zip"
    assert not folder.exists() and (tmp_path / "Deck_Example.png").exists()
    with zipfile.ZipFile(tmp_path / "Deck_Example-paraview.zip") as z:
        assert sorted(z.namelist()) == ["Deck_Example-paraview/Deck.exo", "Deck_Example-paraview/open.py",
                                        "Deck_Example-paraview/state.pvsm"]
    assert bundles["Deck_Example"]["bytes"] == (tmp_path / "Deck_Example-paraview.zip").stat().st_size


def test_bundles_stay_out_of_the_gallery_archive(tmp_path):
    for name in ("A.png", "A-thumb.png", "A.vtksz", "A-paraview.zip"):
        (tmp_path / name).write_bytes(b"x")
    manifest = {"views": [{"id": "A", "png": "A.png", "thumb": "A-thumb.png", "scene": "A.vtksz",
                           "paraview": {"file": "A-paraview.zip", "bytes": 1}}]}
    names = bg.pack(tmp_path, manifest, tmp_path / "gallery.zip")
    assert "A-paraview.zip" not in names and "A.vtksz" in names
    assert json.loads((tmp_path / "manifest.json").read_text())["views"][0]["paraview"]["file"] == "A-paraview.zip"


def test_page_links_the_bundle_or_says_why_not():
    manifest = {"paraview_url": "https://github.com/x/y/releases/download/gallery-assets/"}
    meta = {"paraview": {"file": "A-paraview.zip", "bytes": 3_100_000}}
    md = mg.paraview_markdown(meta, manifest)
    assert "[Open in ParaView](https://github.com/x/y/releases/download/gallery-assets/A-paraview.zip)" in md
    assert "(3.1 MB)" in md and "pvpython open.py" in md
    assert mg.paraview_markdown({"paraview": {"skipped": "the output file is 51 MB"}}, manifest) == \
        " No ParaView bundle for this example: the output file is 51 MB."
    assert mg.paraview_markdown({}, manifest) == ""
    assert mg.paraview_markdown(meta, {}) == ""                      # no release URL known: no dead link
