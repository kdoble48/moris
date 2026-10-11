"""Generate the example-gallery scenes and (optionally) upload them as a GitHub release asset.

    python doc/gallery/build_gallery.py --build <MORIS build tree with examples run>
        [--out DIR] [--only ID ...] [--pvpython PATH] [--vendor] [--upload]

Steps: list every projects/EXA example (the catalog's own scan); collect the Exodus files its
ctest run left in ``<build>/projects/EXA/...``; export one scene per example with pvpython
(``export_scenes.py``); crop thumbnails; validate scenes; write manifest.json; zip the listed files
to ``<out>.zip``; with --upload, replace the asset on the ``gallery-assets`` release. Slow
(pvpython): run it in the background. Run the examples first (``ctest -L``/``-R`` in a build with
MORIS_USE_EXAMPLES=ON); examples without output are skipped.

Adapted from kdoble48/nastran-pv (docs/build_gallery.py).
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw, ImageFont

HERE = Path(__file__).resolve().parent
DOC = HERE.parent
ROOT = DOC.parent
sys.path.insert(0, str(DOC / "_ext"))

RELEASE = "gallery-assets"
REPO = "kdoble48/moris"
SCENE_WARN_BYTES = 3_000_000
CATEGORY_ORDER = ["optimization", "structure", "thermal", "fluid", "mesh_generation"]
WINDOWS_PVPYTHON = Path("/mnt/c/Program Files/ParaView 6.0.0/bin/pvpython.exe")
VENDOR_PACKAGES = ["trame-vtk==2.11.17", "trame-client==4.1.4"]
# a decomposed parallel output (name.exo.4.0, name.exo.e-s.0001.4.0) is not a standalone file
PARALLEL_PIECE = re.compile(r"\.\d+\.\d+$")


def exodus_files(directory: Path) -> list[Path]:
    """Single-file Exodus outputs in an example's build directory (and its bin/)."""
    files = []
    for path in sorted(directory.rglob("*")):
        name = path.name
        if not path.is_file() or PARALLEL_PIECE.search(name):
            continue
        if name.endswith((".exo", ".e")) or ".exo.e-s." in name or ".e-s." in name:
            files.append(path)
    return files


def default_pvpython() -> str:
    if WINDOWS_PVPYTHON.exists():
        return str(WINDOWS_PVPYTHON)
    return shutil.which("pvpython") or "pvpython"


def to_tool_path(path: Path | str, pvpython: str) -> str:
    """The path as the pvpython executable sees it (Windows path for a Windows pvpython under WSL)."""
    if pvpython.lower().endswith(".exe"):
        return subprocess.run(["wslpath", "-w", str(path)], capture_output=True, text=True, check=True).stdout.strip()
    return str(path)


def page_title(example: dict) -> str:
    """The example page's ``# `` heading, so gallery cards and pages agree; the catalog title otherwise."""
    page = DOC / "examples" / "exa" / example["category"] / f"{example['leaf']}.md"
    if page.exists():
        for line in page.read_text(encoding="utf-8").splitlines():
            if line.startswith("# "):
                return line[2:].strip()
    return example["title"]


def build_examples(build: Path, only: list[str] | None, overrides: dict) -> list[dict]:
    from exa_catalog import collect_examples, ctest_names

    examples = []
    for example in collect_examples(ROOT):
        if only and example["leaf"] not in only:
            continue
        files = exodus_files(build / example["directory"])
        if not files:
            continue
        tests = ctest_names(example)
        examples.append({
            "id": example["leaf"], "directory": example["directory"], "category": example["category"],
            "title": page_title(example), "files": [str(f) for f in files],
            "override": overrides.get(example["leaf"], {}),
            "reproduce": [f"ctest -R '^{tests[0]}$' --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON"
                          if tests else "# not registered as a ctest",
                          f"python3 doc/gallery/build_gallery.py --build <build tree> --only {example['leaf']}"],
        })
    return examples


def crop_thumbnail(src, dst, size=(800, 500), margin=0.06, threshold=8) -> None:
    """Crop to the non-white pixels, pad to the ``size`` aspect ratio with white, resize."""
    im = Image.open(src).convert("RGB")
    diff = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).convert("L")
    x0, y0, x1, y1 = diff.point(lambda v: 255 if v > threshold else 0).getbbox() or (0, 0, *im.size)
    w, h = (x1 - x0) * (1 + 2 * margin), (y1 - y0) * (1 + 2 * margin)
    aspect = size[0] / size[1]
    w, h = (h * aspect, h) if w / h < aspect else (w, w / aspect)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    box = [round(cx - w / 2), round(cy - h / 2), round(cx + w / 2), round(cy + h / 2)]
    canvas = Image.new("RGB", (box[2] - box[0], box[3] - box[1]), (255, 255, 255))
    part = im.crop((max(box[0], 0), max(box[1], 0), min(box[2], im.width), min(box[3], im.height)))
    canvas.paste(part, (max(-box[0], 0), max(-box[1], 0)))
    canvas.resize(size, Image.LANCZOS).save(dst)


def validate_scene(path) -> set[str]:
    """Class names in a .vtksz; ValueError unless it has an actor, polydata and data arrays."""
    with zipfile.ZipFile(path) as z:
        state = json.loads(z.read("index.json"))
    classes, indexed, stack = set(), False, [state["scene"]]
    while stack:
        node = stack.pop()
        classes.add(node.get("type"))
        indexed |= bool(node.get("properties", {}).get("indexedLookup"))
        stack.extend(node.get("dependencies", []))
    if not any(c and c.endswith("Actor") for c in classes):
        raise ValueError(f"{Path(path).name}: no actors in scene")
    if indexed:
        raise ValueError(f"{Path(path).name}: indexed lookup in scene (vtk.js cannot render it)")
    if "vtkPolyData" not in classes or not state.get("hashes"):
        raise ValueError(f"{Path(path).name}: no geometry data in scene")
    return classes


def collect_manifest(out_dir, commit, build_label) -> dict:
    views = [json.loads(p.read_text(encoding="utf-8")) for p in Path(out_dir).glob("*.json")
             if p.name not in ("manifest.json", "examples.json")]
    rank = {c: i for i, c in enumerate(CATEGORY_ORDER)}
    views.sort(key=lambda v: (rank.get(v["category"], 99), v["title"].lower()))
    now = dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")
    return {"generated": now, "commit": commit, "build": build_label, "views": views}


def listed_files(manifest) -> list[str]:
    names = set()
    for v in manifest["views"]:
        names.update((v["png"], v["thumb"], v["scene"]))
        if v.get("animation"):
            names.add(v["animation"]["file"])
    return sorted(names)


def make_animation(frame_paths, dst, labels=None, width=960, frame_ms=250, hold_ms=1500, margin=0.04,
                   threshold=8) -> None:
    """Animated WebP from PNG frames, cropped to the union of their non-white content; ``labels`` are
    stamped in the top-left corner of each frame (above the content: the crop leaves a label band)."""
    frames = [Image.open(p).convert("RGB") for p in frame_paths]
    box = None
    for im in frames:
        diff = ImageChops.difference(im, Image.new("RGB", im.size, (255, 255, 255))).convert("L")
        b = diff.point(lambda v: 255 if v > threshold else 0).getbbox()
        if b:
            box = b if box is None else (min(box[0], b[0]), min(box[1], b[1]), max(box[2], b[2]), max(box[3], b[3]))
    box = box or (0, 0, *frames[0].size)
    pad_x, pad_y = round((box[2] - box[0]) * margin), round((box[3] - box[1]) * margin)
    box = (max(box[0] - pad_x, 0), max(box[1] - pad_y, 0),
           min(box[2] + pad_x, frames[0].width), min(box[3] + pad_y, frames[0].height))
    scale = min(1.0, width / (box[2] - box[0]))
    size = (round((box[2] - box[0]) * scale), round((box[3] - box[1]) * scale))
    frames = [im.crop(box).resize(size, Image.LANCZOS) for im in frames]
    if labels:
        band = 34
        font = ImageFont.load_default(size=22)
        stamped = []
        for im, text in zip(frames, labels):
            canvas = Image.new("RGB", (im.width, im.height + band), (255, 255, 255))
            canvas.paste(im, (0, band))
            ImageDraw.Draw(canvas).text((6, 5), text, fill=(40, 40, 40), font=font)
            stamped.append(canvas)
        frames = stamped
    durations = [frame_ms] * (len(frames) - 1) + [hold_ms]
    frames[0].save(dst, save_all=True, append_images=frames[1:], duration=durations, loop=0, quality=80, method=4)


def pack(out_dir, manifest, zip_path) -> list[str]:
    out_dir = Path(out_dir)
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=1), encoding="utf-8")
    names = listed_files(manifest) + ["manifest.json"]
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for name in names:
            z.write(out_dir / name, name)
    return names


def _git_commit(path) -> str:
    r = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=path, capture_output=True, text=True)
    return r.stdout.strip() or "unknown"


def vendor() -> None:
    """Install trame-vtk (pure Python) where pvpython's scene_export.py finds it."""
    subprocess.run([sys.executable, "-m", "pip", "install", "--quiet", "--no-deps", "--upgrade",
                    "--target", str(HERE / "_vendor"), *VENDOR_PACKAGES], check=True)


def upload(zip_path) -> None:
    if subprocess.run(["gh", "release", "view", RELEASE, "--repo", REPO], capture_output=True).returncode:
        subprocess.run(["gh", "release", "create", RELEASE, "--repo", REPO, "--title", "Gallery assets",
                        "--latest=false", "--notes",
                        "Generated by doc/gallery/build_gallery.py; downloaded by the docs workflows."], check=True)
    subprocess.run(["gh", "release", "upload", RELEASE, str(zip_path), "--repo", REPO, "--clobber"], check=True)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build", type=Path, required=True, help="MORIS build tree whose examples were run")
    runs = os.environ.get("MORIS_RUNS_DIR")
    ap.add_argument("--out", type=Path,
                    default=Path(runs) / "docs_gallery" if runs else DOC / "_build" / "gallery")
    ap.add_argument("--only", nargs="*", help="example directory names to export (default: all)")
    ap.add_argument("--pvpython", default=default_pvpython())
    ap.add_argument("--vendor", action="store_true", help="(re)install trame-vtk into doc/gallery/_vendor first")
    ap.add_argument("--upload", action="store_true", help=f"replace the asset on the '{RELEASE}' release")
    args = ap.parse_args(argv)

    if args.vendor or not (HERE / "_vendor" / "trame_vtk").exists():
        vendor()
    out = args.out.resolve()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    overrides_file = HERE / "overrides.json"
    overrides = json.loads(overrides_file.read_text()) if overrides_file.exists() else {}
    examples = build_examples(args.build.resolve(), args.only, overrides)
    for example in examples:
        example["files"] = [to_tool_path(f, args.pvpython) for f in example["files"]]
    examples_json = out / "examples.json"
    examples_json.write_text(json.dumps(examples, indent=1), encoding="utf-8")
    print(f"{len(examples)} examples with output under {args.build}")
    # one pvpython process per example: a reader crash (IOSS aborts on some meshes) loses one scene, not all
    command = [args.pvpython, to_tool_path(HERE / "export_scenes.py", args.pvpython),
               to_tool_path(examples_json, args.pvpython), to_tool_path(out, args.pvpython)]
    log = out / "export.log"
    failed = []
    with open(log, "w", encoding="utf-8") as log_file:
        for example in examples:
            result = subprocess.run(command + [example["id"]], stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                    text=True, errors="replace")
            log_file.write(f"=== {example['id']} (exit {result.returncode})\n{result.stdout}\n")
            status = [line for line in result.stdout.splitlines() if line.startswith(("exported", "skip"))]
            print(status[-1] if status else f"failed {example['id']} (exit {result.returncode})", flush=True)
            if not status:
                failed.append(example["id"])
    if failed:
        print(f"{len(failed)} examples crashed pvpython (see {log}): {' '.join(failed)}")
    if not any(out.glob("*.vtksz")):
        print("nothing exported")
        return 1
    build_label = f"{args.build.resolve().name} @ {_git_commit(args.build.resolve().parent)}"
    manifest = collect_manifest(out, _git_commit(ROOT), build_label)
    for v in manifest["views"]:
        crop_thumbnail(out / v["thumb_raw"], out / v["thumb"])
        if v.get("animation"):
            make_animation([out / "frames" / f for f in v["animation"]["frames"]], out / v["animation"]["file"],
                           v["animation"].get("labels"))
        validate_scene(out / v["scene"])
        size = (out / v["scene"]).stat().st_size
        if size > SCENE_WARN_BYTES:
            print(f"warning: {v['scene']} is {size / 1e6:.1f} MB (budget {SCENE_WARN_BYTES / 1e6:.0f} MB)")
    zip_path = out.with_suffix(".zip")
    names = pack(out, manifest, zip_path)
    print(f"{len(manifest['views'])} scenes, {len(names)} files -> {zip_path} ({zip_path.stat().st_size / 1e6:.1f} MB)")
    if args.upload:
        upload(zip_path)
        print(f"uploaded to release '{RELEASE}' on {REPO}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
