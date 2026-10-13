"""Attached views: where a view script is looked up, and the overrides.json key check."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "gallery"))

import build_gallery as bg  # noqa: E402


def test_view_script_lookup(tmp_path):
    example = {"leaf": "Deck_Example", "directory": "projects/EXA/optimization/Deck_Example"}
    views, source = tmp_path / "views", tmp_path / "src"
    assert bg.view_script_for(example, source, views) is None
    (source / example["directory"]).mkdir(parents=True)
    (source / example["directory"] / "gallery_view.py").write_text("def view(scene): pass\n")
    assert bg.view_script_for(example, source, views).name == "gallery_view.py"
    views.mkdir()
    (views / "Deck_Example.py").write_text("def view(scene): pass\n")
    assert bg.view_script_for(example, source, views) == views / "Deck_Example.py"   # the docs copy wins


def test_override_keys_are_checked():
    problems = bg.check_overrides({"X": {"camera": {"zoom": 1.2}, "_why": "a note", "cammera": {}},
                                  "Y": {"clip": {"normal": [0, 1, 0]}}})
    assert len(problems) == 1 and problems[0].startswith("X: unknown override 'cammera'")


def test_shipped_overrides_and_views_are_valid():
    overrides = bg.json.loads((bg.HERE / "overrides.json").read_text())
    assert bg.check_overrides(overrides) == []
    for script in (bg.HERE / "views").glob("*.py"):
        assert "def view(" in script.read_text()
