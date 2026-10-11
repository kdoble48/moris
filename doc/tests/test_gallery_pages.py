"""Page markup for the example scenes (doc/_ext/moris_gallery.py); no Sphinx needed."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_ext"))

import moris_gallery as mg  # noqa: E402

META = {
    "id": "Laplace_2D", "directory": "projects/EXA/thermal/diffusion/Laplace_2D", "category": "thermal",
    "title": "Laplace 2D", "kind": "result", "array": "TEMP", "caption": "Solution field TEMP from the VIS output",
    "png": "Laplace_2D.png", "thumb": "Laplace_2D-thumb.png", "scene": "Laplace_2D.vtksz", "scene_bytes": 250_000,
    "text": ["Laplace.exo.e-s.0000", "100 cells, 2D"],
    "legends": [{"array": "TEMP", "kind": "continuous", "range": [0.0, 2.0],
                 "stops": [[0.0, 0.0, 0.0, 1.0], [2.0, 1.0, 0.0, 0.0]]}],
    "source": {"file": "Laplace.exo.e-s.0000", "modified": "2026-07-01"},
}
MANIFEST = {"generated": "2026-10-11T14:02:05+00:00", "commit": "abc1234", "build": "build_opt @ 3bc5a1c",
            "views": [META]}


def test_root_prefix():
    assert mg.root_prefix("index") == ""
    assert mg.root_prefix("examples/exa/thermal/Laplace_2D") == "../../../"


def test_version_tag():
    assert mg.version_tag(MANIFEST) == "20261011T140205"
    assert mg.version_tag(None) == ""


def test_scene_markdown_paths_and_cache_tag():
    md = mg.scene_markdown(META, "examples/exa/thermal/Laplace_2D", "20261011T140205", "build_opt @ 3bc5a1c")
    assert "````{tab-item} Interactive" in md and "````{tab-item} Static" in md
    assert 'src="../../../_static/viewer.html?fileURL=gallery/Laplace_2D.vtksz%3Fv%3D20261011T140205"' in md
    assert 'src="../../../_static/gallery/Laplace_2D.png?v=20261011T140205"' in md
    assert "Scene: 250 KB." in md
    assert "Rendered from `Laplace.exo.e-s.0000`, written 2026-07-01 by build_opt @ 3bc5a1c." in md
    assert "linear-gradient(90deg, #0000ff 0.0%, #ff0000 100.0%)" in md


def test_categorical_legend_collapses_long_lists():
    entries = [{"value": str(i), "label": f"phase {i}", "rgb": [0.5, 0.5, 0.5]} for i in range(13)]
    out = mg.legend_html({"array": "object_id", "kind": "categorical", "entries": entries})
    assert "<details><summary>13 values</summary>" in out
    short = mg.legend_html({"array": "object_id", "kind": "categorical", "entries": entries[:3]})
    assert "<details>" not in short and short.count("moris-swatch") == 3


def test_gallery_markdown():
    md = mg.gallery_markdown(MANIFEST)
    assert "## Thermal" in md
    assert ":link: /examples/exa/thermal/Laplace_2D" in md
    assert ":img-top: ../_static/gallery/Laplace_2D-thumb.png" in md
    assert "outputs from build_opt @ 3bc5a1c" in md
    assert "have not been generated" in mg.gallery_markdown(None)


def test_scene_for_uses_loaded_manifest(tmp_path):
    (tmp_path / "_static" / "gallery").mkdir(parents=True)
    import json

    (tmp_path / "_static" / "gallery" / "manifest.json").write_text(json.dumps(MANIFEST))
    mg.load(tmp_path)
    assert "viewer.html" in mg.scene_for(META["directory"], "examples/exa/thermal/Laplace_2D")
    assert mg.scene_for("projects/EXA/none", "x") == ""
    mg.load(tmp_path / "missing")
    assert mg.scene_for(META["directory"], "x") == ""
