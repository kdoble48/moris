"""Export one interactive scene per MORIS example (run with pvpython, ParaView 6).

    pvpython doc/gallery/export_scenes.py EXAMPLES.json OUT_DIR

EXAMPLES.json (written by ``doc/gallery/build_gallery.py``) is a list of
{"id", "directory", "category", "title", "files": [Exodus paths], "override": {...}, "reproduce": [...]}.
For each example the file and array to show are chosen by ``gallery_select.choose`` from what the
example's ctest run left behind; the last time step is rendered. Writes <id>.png (with the scalar bar),
<id>-thumb-raw.png (overlays hidden; cropped later), <id>.vtksz (vtk.js scene) and <id>.json
(metadata). Examples with nothing to show are skipped. Exit status 1 only when nothing was exported.

Adapted from kdoble48/nastran-pv (paraview/export_site.py).
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from paraview import simple as pv  # noqa: E402

import gallery_select as gs  # noqa: E402
import scene_export as se  # noqa: E402

SIZE = [1600, 1000]
EDGE_CELL_LIMIT = 4000           # draw cell edges only where they stay legible
XTK_EDGE_CELL_LIMIT = 20000      # on the XTK integration mesh the edges show the immersed cut
PHASE_COLORS = [                 # one color per phase index; phase 0 light gray
    [0.86, 0.86, 0.86], [0.17, 0.42, 0.69], [0.84, 0.38, 0.30], [0.30, 0.62, 0.36],
    [0.58, 0.40, 0.74], [0.93, 0.69, 0.13], [0.35, 0.70, 0.80], [0.55, 0.34, 0.29],
]
KIND_CAPTION = {
    "result": "Solution field {label} from the VIS output",
    "gen": "Geometry field {label} from the GEN output",
    "xtk": "XTK integration mesh, colored by element block (phase)",
    "mesh": "Mesh",
}


def new_view():
    view = pv.CreateRenderView()
    view.ViewSize = SIZE
    view.Background = [1, 1, 1]
    try:
        view.UseColorPaletteForBackground = 0
    except AttributeError:
        pass
    view.OrientationAxesVisibility = 1
    return view


def open_reader(path):
    """IOSS reader (ParaView 6 default); the legacy Exodus reader when IOSS cannot read the file."""
    try:
        reader = pv.IOSSReader(FileName=[path])
        reader.UpdatePipelineInformation()
        if reader.GetDataInformation() is not None:
            return reader
    except Exception:  # IOSS rejects some MORIS meshes (e.g. duplicate side-set names)
        pass
    reader = pv.ExodusIIReader(FileName=[path])
    reader.UpdatePipelineInformation()
    return reader


def file_info(path):
    """A gallery_select.FileInfo for ``path`` (arrays and point count at the last time step)."""
    reader = open_reader(path)
    times = list(getattr(reader, "TimestepValues", []) or [])
    reader.UpdatePipeline(times[-1] if times else 0.0)
    info = gs.FileInfo(path, list(reader.PointData.keys()), reader.GetDataInformation().GetNumberOfPoints(),
                       len(times))
    pv.Delete(reader)
    return info


def is_3d(bounds):
    return (bounds[5] - bounds[4]) > 1e-9 * max(1.0, bounds[1] - bounds[0], bounds[3] - bounds[2])


def set_camera(view, three_d, face_cut=False):
    """2D: straight down z, parallel projection. 3D: isometric, from +y when a y-cutaway face must show."""
    view.ResetCamera()
    cam = view.GetActiveCamera()
    focal = cam.GetFocalPoint()
    dist = cam.GetDistance()
    if three_d:
        view.CameraParallelProjection = 0
        side = 1.0 if face_cut else -1.0
        cam.SetPosition(focal[0] + dist / math.sqrt(3), focal[1] + side * dist / math.sqrt(3),
                        focal[2] + dist / math.sqrt(3))
        cam.SetViewUp(0, 0, 1)
    else:
        view.CameraParallelProjection = 1
        cam.SetPosition(focal[0], focal[1], focal[2] + dist)
        cam.SetViewUp(0, 1, 0)
    view.ResetCamera()
    cam.Zoom(1.08)


def color_blocks(rep, reader, view):
    """Categorical coloring by element block, one color per phase, block names as labels."""
    blocks = list(getattr(reader, "ElementBlocks", []) or [])
    pv.ColorBy(rep, ("CELLS", "object_id"))
    lut = pv.GetColorTransferFunction("object_id")
    lut.InterpretValuesAsCategories = 1
    annotations, colors = [], []
    for index, name in enumerate(blocks, start=1):
        phase, label = gs.block_phase(name)
        annotations += [str(index), label]
        colors += PHASE_COLORS[(phase or 0) % len(PHASE_COLORS)]
    lut.Annotations = annotations
    lut.IndexedColors = colors
    rep.SetScalarBarVisibility(view, True)


def color_field(rep, source, choice, view):
    lut_name = choice.label if choice.vector else choice.array
    pv.ColorBy(rep, ("POINTS", lut_name))
    lut = pv.GetColorTransferFunction(lut_name)
    data_range = source.PointData[lut_name].GetRange(-1)
    signed = data_range[0] < 0 < data_range[1]
    if signed:
        lut.ApplyPreset("Cool to Warm", True)
        bound = max(abs(data_range[0]), abs(data_range[1]))
        lut.RescaleTransferFunction(-bound, bound)
    else:
        lut.ApplyPreset("Viridis (matplotlib)", True)
        lut.RescaleTransferFunction(*data_range)
    rep.SetScalarBarVisibility(view, True)


def render(example, choice, out):
    """Build the view for ``choice``; return (view, caption text, step info)."""
    reader = open_reader(choice.path)
    times = list(getattr(reader, "TimestepValues", []) or [])
    view = new_view()
    if times:
        view.ViewTime = times[-1]
    reader.UpdatePipeline(times[-1] if times else 0.0)
    bounds = reader.GetDataInformation().GetBounds()
    three_d = is_3d(bounds)
    blocks = len(list(getattr(reader, "ElementBlocks", []) or [])) > 1 and "object_id" in reader.CellData.keys()
    by_block = choice.kind == "xtk" or (choice.kind == "mesh" and blocks)
    source = pv.MergeBlocks(Input=reader)
    if choice.vector:
        source = pv.Calculator(Input=source, ResultArrayName=choice.label, Function=f"mag({choice.array})")
    cutaway = None
    if three_d and by_block:
        # the outer phase hides the inner ones: keep the half y <= mid-plane
        cutaway = 0.5 * (bounds[2] + bounds[3])
        source = pv.Clip(Input=source, ClipType="Plane", Crinkleclip=0)
        source.ClipType.Origin = [0.5 * (bounds[0] + bounds[1]), cutaway, 0.5 * (bounds[4] + bounds[5])]
        source.ClipType.Normal = [0.0, 1.0, 0.0]
        source.Invert = 1
    if three_d:
        source = pv.ExtractSurface(Input=source)
    source.UpdatePipeline(times[-1] if times else 0.0)
    rep = pv.Show(source, view)
    cells = source.GetDataInformation().GetNumberOfCells()
    if cells <= (XTK_EDGE_CELL_LIMIT if choice.kind == "xtk" else EDGE_CELL_LIMIT):
        rep.SetRepresentationType("Surface With Edges")
        rep.EdgeColor = [0.55, 0.55, 0.55]
    if by_block:
        color_blocks(rep, reader, view)
    elif choice.kind == "mesh":
        rep.DiffuseColor = [0.75, 0.80, 0.86]
        rep.AmbientColor = [0.75, 0.80, 0.86]
    else:
        color_field(rep, source, choice, view)
    set_camera(view, three_d, face_cut=cutaway is not None)
    step = f"time step {len(times)} of {len(times)}" if len(times) > 1 else None
    if cutaway is not None:
        step = (step + ", " if step else "") + f"cut away at y = {cutaway:.4g}"
    return view, three_d, cells, step, by_block


def export_example(example, out):
    override = example.get("override", {})
    if override.get("skip"):
        return None, "skipped by override"
    infos = [file_info(path) for path in example["files"]]
    choice = gs.choose(infos, override)
    if choice is None:
        return None, "no file with geometry"
    view, three_d, cells, step, by_block = render(example, choice, out)
    base = example["id"]
    file_name = os.path.basename(choice.path.replace("\\", "/"))
    caption = KIND_CAPTION[choice.kind].format(label=choice.label)
    if choice.kind == "mesh" and by_block:
        caption = "Mesh, colored by element block"
    summary = [f"{file_name}" + (f", {step}" if step else ""), f"{cells} cells, {'3D' if three_d else '2D'}"]
    over = se.overlays(view)
    pv.SaveScreenshot(os.path.join(out, f"{base}.png"), view, ImageResolution=SIZE)
    se.write_vtksz(view, os.path.join(out, f"{base}.vtksz"))
    meta = {
        "id": base, "directory": example["directory"], "category": example["category"],
        "title": example["title"], "kind": choice.kind, "array": choice.label, "caption": caption,
        "png": f"{base}.png", "thumb": f"{base}-thumb.png", "thumb_raw": f"{base}-thumb-raw.png",
        "scene": f"{base}.vtksz", "scene_bytes": os.path.getsize(os.path.join(out, f"{base}.vtksz")),
        "text": summary, "legends": over["legends"],
        "source": {"file": file_name, "modified": time.strftime(
            "%Y-%m-%d", time.localtime(os.path.getmtime(choice.path)))},
        "reproduce": example.get("reproduce", []),
    }
    se.hide_overlays(view)
    pv.SaveScreenshot(os.path.join(out, meta["thumb_raw"]), view, ImageResolution=SIZE)
    pv.Delete(view)
    for proxy in list(pv.GetSources().values()):
        pv.Delete(proxy)
    return meta, None


def main(argv=None):
    """``export_scenes.py EXAMPLES.json OUT_DIR [ID ...]`` — only the listed examples when IDs are given."""
    argv = sys.argv[1:] if argv is None else argv
    examples_json, out, only = argv[0], argv[1], set(argv[2:])
    os.makedirs(out, exist_ok=True)
    with open(examples_json, encoding="utf-8") as f:
        examples = [e for e in json.load(f) if not only or e["id"] in only]
    written = 0
    for example in examples:
        try:
            meta, reason = export_example(example, out)
        except Exception as exc:  # one broken output must not stop the gallery
            meta, reason = None, f"error: {exc}"
        if meta is None:
            print(f"skip {example['id']}: {reason}", flush=True)
            continue
        with open(os.path.join(out, f"{example['id']}.json"), "w", encoding="utf-8") as f:
            json.dump(meta, f, indent=1, allow_nan=False)
        written += 1
        print(f"exported {example['id']} ({meta['kind']}, {meta['array']}, {meta['scene_bytes'] / 1e6:.2f} MB)",
              flush=True)
    if not written:
        print("nothing exported")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
