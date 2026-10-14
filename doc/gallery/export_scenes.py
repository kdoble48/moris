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
import traceback

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
    "xtk": "XTK integration mesh, colored by phase: cut cells subdivided into triangles along the interface, "
           "uncut cells whole",
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
    """IOSS reader (ParaView 6 default) on this one file; the legacy Exodus reader when IOSS cannot read it.

    Related files are not scanned: an ``.e-s.NNNN`` series is read file by file (see
    ``gallery_select.series_paths``), since XTK renumbers element blocks between iterations."""
    try:
        reader = pv.IOSSReader(FileName=[path])
        reader.ScanForRelatedFiles = 0
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


# Keeps the boundary edges of one phase's cells that have the phase on one side only: on an adaptively
# refined mesh, the plain boundary edges also include the seams at hanging nodes (fine cells meeting
# one coarse cell), which lie inside the phase. 2D only.
OUTLINE_SCRIPT = r"""
import vtk
inp = self.GetInput()
out = self.GetOutput()
surf = vtk.vtkDataSetSurfaceFilter()
surf.SetInputData(inp)
edges_filter = vtk.vtkFeatureEdges()
edges_filter.SetInputConnection(surf.GetOutputPort())
edges_filter.BoundaryEdgesOn()
edges_filter.FeatureEdgesOff()
edges_filter.NonManifoldEdgesOff()
edges_filter.ManifoldEdgesOff()
edges_filter.ColoringOff()
edges_filter.Update()
edges = edges_filter.GetOutput()
out.SetPoints(edges.GetPoints())
lines = vtk.vtkCellArray()
if inp.GetNumberOfCells() > 0:
    locator = vtk.vtkCellLocator()
    locator.SetDataSet(inp)
    locator.BuildLocator()
    pts = edges.GetPoints()
    for i in range(edges.GetNumberOfCells()):
        ids = edges.GetCell(i).GetPointIds()
        a, b = pts.GetPoint(ids.GetId(0)), pts.GetPoint(ids.GetId(1))
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = (dx * dx + dy * dy) ** 0.5
        if length == 0.0:
            continue
        nx, ny = -dy / length * 0.05 * length, dx / length * 0.05 * length
        mx, my, mz = 0.5 * (a[0] + b[0]), 0.5 * (a[1] + b[1]), 0.5 * (a[2] + b[2])
        if (locator.FindCell([mx + nx, my + ny, mz]) >= 0) != (locator.FindCell([mx - nx, my - ny, mz]) >= 0):
            lines.InsertNextCell(2)
            lines.InsertCellPoint(ids.GetId(0))
            lines.InsertCellPoint(ids.GetId(1))
out.SetLines(lines)
"""


def available_blocks(reader):
    try:
        return list(reader.ElementBlocks.Available)
    except AttributeError:
        return list(getattr(reader, "ElementBlocks", []) or [])


def select_phases(reader, phases):
    """Read only the element blocks of ``phases`` (None: all). True when the selection took effect."""
    if not phases:
        return False
    keep = [n for n in available_blocks(reader) if gs.block_phase(n)[0] in phases]
    if not keep:
        return False
    reader.ElementBlocks = keep
    return True


def set_file(readers, path, view):
    """Point every reader at ``path`` (one file of a series), redo its block selection, follow its time."""
    for reader, phases in readers:
        reader.FileName = [path]
        reader.UpdatePipelineInformation()
        select_phases(reader, phases)
    times = list(getattr(readers[0][0], "TimestepValues", []) or [])
    view.ViewTime = times[-1] if times else 0.0
    return view.ViewTime


def add_outline(view, path, phase, bounds):
    """Black outline of one bulk phase: its interfaces with void and with other phases, and its outer edge."""
    reader = open_reader(path)
    select_phases(reader, {phase})
    merged = pv.MergeBlocks(Input=reader)
    merged.MergePoints = 1
    outline = pv.ProgrammableFilter(Input=merged, OutputDataSetType="vtkPolyData", Script=OUTLINE_SCRIPT)
    lifted = pv.Transform(Input=outline)
    lifted.Transform.Translate = [0.0, 0.0, 1e-3 * max(bounds[1] - bounds[0], bounds[3] - bounds[2])]
    rep = pv.Show(lifted, view)
    rep.ColorArrayName = ["POINTS", ""]
    rep.DiffuseColor = [0.0, 0.0, 0.0]
    rep.AmbientColor = [0.0, 0.0, 0.0]
    rep.LineWidth = 3.0
    return reader


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


def background_edges(source, view):
    """The edges of the quad/hex cells of ``source`` as a gray wireframe over the colored surface."""
    # the quad/hex family in every order: MORIS writes 9-node quads / 27-node hexes for quadratic meshes
    wanted = {"Quad", "Pixel", "Quadratic Quadrilateral", "Bi-Quadratic Quadrilateral", "Lagrange Quadrilateral",
              "Hexahedron", "Voxel", "Quadratic Hexahedron", "Tri-Quadratic Hexahedron",
              "Bi-Quadratic Quadratic Hexahedron", "Lagrange Hexahedron"}
    quads = pv.ExtractCellsByType(Input=source)
    quads.CellTypes = [t for t in quads.GetProperty("CellTypes").Available if t in wanted]
    rep = pv.Show(quads, view)
    rep.SetRepresentationType("Wireframe")
    rep.ColorArrayName = ["POINTS", ""]
    rep.DiffuseColor = [0.55, 0.55, 0.55]
    rep.AmbientColor = [0.55, 0.55, 0.55]
    rep.LineWidth = 1.0
    return rep


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
    override = example.get("override", {})
    reader = open_reader(choice.path)
    # bulk phases only: void phases (often present in the VIS output, with an ersatz solution) are dropped
    phases, hidden = None, []
    if choice.kind in ("result", "gen"):
        phases = set(override["phases"]) if override.get("phases") else \
            gs.bulk_phases(example.get("deck_texts", []), choice.array, choice.vector)
        present = {gs.block_phase(n)[0] for n in available_blocks(reader)} - {None}
        if not select_phases(reader, phases):
            phases = None
        else:
            hidden = sorted(present - phases)
    readers = [(reader, phases)]
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
    volume = source            # the 3D dataset before its surface is extracted: what clips act on
    if three_d:
        source = pv.ExtractSurface(Input=source)
    source.UpdatePipeline(times[-1] if times else 0.0)
    rep = pv.Show(source, view)
    cells = source.GetDataInformation().GetNumberOfCells()
    # override "edges": true draws the cells whatever their number (adaptive refinement is the point there)
    edge_rep = None
    if example.get("override", {}).get("edges") or \
            cells <= (XTK_EDGE_CELL_LIMIT if choice.kind == "xtk" else EDGE_CELL_LIMIT):
        if choice.kind in ("result", "gen"):
            # background cells only: the cut cells' triangles/tets would draw as a dense dark band along
            # every interface, which the black outline already marks
            edge_rep = background_edges(source, view)
        else:
            rep.SetRepresentationType("Surface With Edges")
            rep.EdgeColor = [0.55, 0.55, 0.55]
    if by_block:
        color_blocks(rep, reader, view)
    if choice.kind == "xtk":
        # name the mesh in the image itself: the triangles are the cut, not the background discretization
        tag = pv.Text(Text="Integration mesh (XTK)")
        tag_rep = pv.Show(tag, view)
        tag_rep.WindowLocation = "Upper Left Corner"
        tag_rep.FontSize = 28
        tag_rep.Color = [0.15, 0.15, 0.15]
    elif choice.kind == "mesh" and not by_block:
        rep.DiffuseColor = [0.75, 0.80, 0.86]
        rep.AmbientColor = [0.75, 0.80, 0.86]
    elif choice.kind in ("result", "gen"):
        color_field(rep, source, choice, view)
    interface = None
    if phases and not three_d and override.get("interface", True) is not False:
        # the bulk phases' own boundaries: material/void interfaces, interfaces between phases, outer edges
        for phase in sorted(phases):
            readers.append((add_outline(view, choice.path, phase, bounds), {phase}))
        interface = "outline"
    else:
        level_set = gs.interface_array(list(reader.PointData.keys()))
        if level_set and choice.kind in ("result", "gen") and not three_d and override.get("interface", True) is not False:
            value = override["interface"] if isinstance(override.get("interface"), (int, float)) \
                and not isinstance(override.get("interface"), bool) else example.get("interface_value", 0.0)
            interface = (level_set, value)
            add_interface(view, reader, level_set, value, bounds)    # built on ``reader``: follows its file
    set_camera(view, three_d, face_cut=cutaway is not None)
    # the attached view: declarative keys in overrides.json, then the example's view script
    scene = Scene(example, choice, view, reader, source, rep, bounds, three_d, times, readers, phases, volume)
    scene.edge_rep = edge_rep
    apply_view_overrides(scene, override)
    script = example.get("view_script")
    if script:
        load_view_script(script).view(scene)
    source = scene.source
    step = None
    series = gs.series_paths(choice.path, example.get("files", []))
    if len(series) > 1:
        step = f"iteration {series_number(choice.path)}, the last of {len(series)} files in the .e-s series"
    elif len(times) > 1:
        step = gs.step_labels(choice.path, [], len(times))[0][-1]
    if cutaway is not None:
        step = (step + ", " if step else "") + f"cut away at y = {cutaway:.4g}"
    if scene.notes:
        step = (step + ", " if step else "") + ", ".join(scene.notes)
    shown = (sorted(phases), hidden) if phases else None
    return view, three_d, cells, step, by_block, source, times, interface, readers, shown, scene.caption_extra


class Scene:
    """What an attached view can change after the default pipeline has been built.

    ``source`` is the displayed dataset (its representation is ``rep``); a helper that inserts a
    filter re-shows the result with the same coloring. The camera is already set to the default
    (straight down z in 2D, isometric in 3D); ``camera()`` adjusts it from there."""

    def __init__(self, example, choice, view, reader, source, rep, bounds, three_d, times, readers, phases,
                 volume=None):
        self.example, self.choice, self.view, self.reader = example, choice, view, reader
        self.source, self.rep, self.bounds, self.three_d, self.times = source, rep, bounds, three_d, times
        self.readers, self.phases = readers, phases
        self.volume = volume if volume is not None else source    # 3D: the cells; 2D: the same as source
        self.notes, self.caption_extra = [], ""

    def _cut(self, filter_source, face_cut=False):
        """Show a filter applied to the volume: in 3D its surface is extracted, and the camera is
        reset to the new bounds (from the +y side when a y-cut face is to be seen)."""
        self.volume = filter_source
        shown = pv.ExtractSurface(Input=filter_source) if self.three_d else filter_source
        self.replace(shown)
        self.bounds = self.source.GetDataInformation().GetBounds()
        set_camera(self.view, self.three_d, face_cut=face_cut)
        return self.source

    @property
    def time(self):
        return self.times[-1] if self.times else 0.0

    def replace(self, new_source):
        """Show ``new_source`` (a filter built on ``self.source``) in place of the current source."""
        coloring = self.rep.ColorArrayName
        representation = self.rep.Representation
        pv.Hide(self.source, self.view)
        self.source = new_source
        self.source.UpdatePipeline(self.time)
        self.rep = pv.Show(self.source, self.view)
        self.rep.Representation = representation
        if coloring and coloring[1]:
            pv.ColorBy(self.rep, (coloring[0], coloring[1]))
            self.rep.SetScalarBarVisibility(self.view, True)
        return self.source

    def clip(self, normal, origin=None, invert=False, crinkle=False, note=True):
        """Keep the half-space on the side ``normal`` points to (``invert`` keeps the other)."""
        b = self.bounds
        origin = origin or [0.5 * (b[0] + b[1]), 0.5 * (b[2] + b[3]), 0.5 * (b[4] + b[5])]
        clip = pv.Clip(Input=self.volume, ClipType="Plane", Crinkleclip=int(crinkle))
        clip.ClipType.Origin = list(origin)
        clip.ClipType.Normal = list(normal)
        clip.Invert = int(invert)
        axis = "xyz"[max(range(3), key=lambda i: abs(normal[i]))]
        if note:
            index = "xyz".index(axis)
            where = origin[index]
            if abs(where) < 1e-9 * max(1.0, b[2 * index + 1] - b[2 * index]):
                where = 0.0    # a mid-plane of a symmetric domain: -2.2e-16 reads as 0
            self.notes.append(f"clipped at {axis} = {where:.4g}")
        # keeping y <= origin (invert with +y normal) exposes a face toward +y: look at it from there
        return self._cut(clip, face_cut=(axis == "y" and bool(invert) == (normal[1] > 0)))

    def threshold(self, array, lower, upper, association="POINTS"):
        """Keep the cells whose ``array`` lies in [lower, upper]."""
        threshold = pv.Threshold(Input=self.volume, Scalars=[association, array])
        threshold.LowerThreshold, threshold.UpperThreshold = lower, upper
        return self._cut(threshold)

    def contour(self, array, values, color=(0.0, 0.0, 0.0), width=3.0):
        """Draw isolines/isosurfaces of a point array over the scene."""
        contour = pv.Contour(Input=self.source, ContourBy=["POINTS", array], Isosurfaces=list(values))
        rep = pv.Show(contour, self.view)
        rep.ColorArrayName = ["POINTS", ""]
        rep.DiffuseColor = rep.AmbientColor = list(color)
        rep.LineWidth = width
        return contour

    def _hide_scalar_bar(self):
        name = self.rep.ColorArrayName[1] if self.rep.ColorArrayName else ""
        if name:
            pv.HideScalarBarIfNotNeeded(pv.GetColorTransferFunction(name), self.view)
            self.rep.SetScalarBarVisibility(self.view, False)

    def color_by(self, array, association="POINTS", preset=None):
        self._hide_scalar_bar()
        pv.ColorBy(self.rep, (association, array))
        lut = pv.GetColorTransferFunction(array)
        if preset:
            lut.ApplyPreset(preset, True)
        self.rep.SetScalarBarVisibility(self.view, True)
        return lut

    def color_blocks(self):
        """Color by element block (one color per phase, block names in the legend) instead of a field."""
        self._hide_scalar_bar()
        color_blocks(self.rep, self.reader, self.view)

    def color_range(self, lower, upper):
        name = self.rep.ColorArrayName[1]
        lut = pv.GetColorTransferFunction(name)
        lut.RescaleTransferFunction(lower, upper)
        lut.AutomaticRescaleRangeMode = "Never"
        return lut

    def edges(self, on=True, color=(0.55, 0.55, 0.55)):
        """Background-cell edges (quads/hexes) on or off; the cut cells stay smooth."""
        if self.choice.kind not in ("result", "gen"):
            self.rep.SetRepresentationType("Surface With Edges" if on else "Surface")
            if on:
                self.rep.EdgeColor = list(color)
            return
        if getattr(self, "edge_rep", None) is not None:
            pv.Hide(self.edge_rep.Input, self.view)
            self.edge_rep = None
        if on:
            self.edge_rep = background_edges(self.source, self.view)
            self.edge_rep.DiffuseColor = self.edge_rep.AmbientColor = list(color)

    def camera(self, azimuth=0.0, elevation=0.0, zoom=1.0, position=None, focal=None, up=None, parallel=None):
        """Adjust the camera from the default: rotate (degrees), zoom (> 1 closer), or set it outright."""
        cam = self.view.GetActiveCamera()
        if position:
            cam.SetPosition(*position)
        if focal:
            cam.SetFocalPoint(*focal)
        if up:
            cam.SetViewUp(*up)
        if parallel is not None:
            self.view.CameraParallelProjection = int(parallel)
        if azimuth:
            cam.Azimuth(azimuth)
        if elevation:
            cam.Elevation(elevation)
            cam.OrthogonalizeViewUp()
        if zoom and zoom != 1.0:
            cam.Zoom(zoom)

    def reset_camera(self):
        set_camera(self.view, self.three_d)

    def caption(self, text):
        """Appended to the page caption."""
        self.caption_extra += (" " if self.caption_extra else "") + text

    def note(self, text):
        """Appended to the one-line summary under the legend."""
        self.notes.append(text)


VIEW_KEYS = {"camera", "clip", "contour", "color_range", "caption", "edges", "threshold"}


def apply_view_overrides(scene, override):
    """The declarative part of an attached view (overrides.json), applied in a fixed order."""
    if "clip" in override:
        c = override["clip"]
        scene.clip(c["normal"], c.get("origin"), c.get("invert", False), c.get("crinkle", False))
    if "threshold" in override:
        t = override["threshold"]
        scene.threshold(t["array"], t["lower"], t["upper"], t.get("association", "POINTS"))
    if "contour" in override:
        c = override["contour"]
        scene.contour(c["array"], c["values"], tuple(c.get("color", (0.0, 0.0, 0.0))), c.get("width", 3.0))
    if "color_range" in override:
        scene.color_range(*override["color_range"])
    if "edges" in override and isinstance(override["edges"], bool) and not override["edges"]:
        scene.edges(False)
    if "camera" in override:
        scene.camera(**override["camera"])
    if "caption" in override:
        scene.caption(override["caption"])


def load_view_script(path):
    """The example's view module (``def view(scene): ...``)."""
    import importlib.util

    spec = importlib.util.spec_from_file_location("gallery_view_" + os.path.basename(path).replace(".", "_"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "view"):
        raise ValueError(f"{path}: a view script must define view(scene)")
    return module


def series_number(path):
    match = gs.SERIES.match(os.path.basename(path.replace("\\", "/")))
    return int(match["index"]) if match else None


def add_interface(view, reader, array, value, bounds):
    """The level-set contour ``array = value`` (the material interface) as a black line, lifted slightly
    toward the camera so it draws over the colored surface; it follows the view time in animations."""
    contour = pv.Contour(Input=pv.MergeBlocks(Input=reader), ContourBy=["POINTS", array], Isosurfaces=[value])
    lift = 1e-3 * max(bounds[1] - bounds[0], bounds[3] - bounds[2])
    lifted = pv.Transform(Input=contour)
    lifted.Transform.Translate = [0.0, 0.0, lift]
    rep = pv.Show(lifted, view)
    pv.ColorBy(rep, None)
    rep.DiffuseColor = [0.0, 0.0, 0.0]
    rep.AmbientColor = [0.0, 0.0, 0.0]
    rep.LineWidth = 3.0
    return lifted


MIN_FRAMES = 3                   # animate only series with at least this many steps
sample_indices = gs.sample_indices


def export_frames(view, source, choice, times, base, out, all_paths, readers):
    """One overlay-free PNG per sampled step, fixed camera, one color range spanning all frames.

    An ``.e-s.NNNN`` series steps through its files (one optimization iteration each), pointing every
    reader at the next file; a single file steps through its time steps. Returns ``{"frames",
    "labels", "unit", "steps", "legend"}`` or None when there are too few steps to animate;
    build_gallery.py stamps the labels and assembles the animation."""
    series = gs.series_paths(choice.path, all_paths)
    if len(series) >= MIN_FRAMES:
        indices = sample_indices(len(series))
        last = series_number(series[-1])
        step_labels = [f"iteration {series_number(p)} of {last}" for p in series]
        iterations, steps = len(series), len(series)

        def go(i):
            source.UpdatePipeline(set_file(readers, series[i], view))
    elif len(times) >= MIN_FRAMES:
        indices = sample_indices(len(times))
        step_labels, iterations = gs.step_labels(choice.path, all_paths, len(times))
        steps = len(times)

        def go(i):
            view.ViewTime = times[i]
            source.UpdatePipeline(times[i])
    else:
        return None
    name = (choice.label if choice.vector else choice.array) if choice.kind in ("result", "gen") else None
    legend = None
    if name:
        lo, hi = math.inf, -math.inf
        for i in indices:
            go(i)
            r = source.PointData[name].GetRange(-1)
            lo, hi = min(lo, r[0]), max(hi, r[1])
        lut = pv.GetColorTransferFunction(name)
        if lo < 0 < hi:
            bound = max(-lo, hi)
            lo, hi = -bound, bound
        lut.RescaleTransferFunction(lo, hi)
        lut.AutomaticRescaleRangeMode = "Never"
        legend = se._legend(name, lut)
    se.hide_overlays(view)
    frames_dir = os.path.join(out, "frames")
    os.makedirs(frames_dir, exist_ok=True)
    files, labels = [], []
    for k, i in enumerate(indices, start=1):
        go(i)
        name_k = f"{base}-frame-{k:03d}.png"
        pv.SaveScreenshot(os.path.join(frames_dir, name_k), view, ImageResolution=SIZE)
        files.append(name_k)
        labels.append(step_labels[i])
    return {"frames": files, "labels": labels, "unit": "iteration" if iterations else "time step",
            "steps": steps, "iterations": iterations, "legend": legend}


def export_example(example, out):
    override = example.get("override", {})
    if override.get("skip"):
        return None, "skipped by override"
    infos = [file_info(path) for path in example["files"]]
    choice = gs.choose(infos, override)
    if choice is None:
        return None, "no file with geometry"
    view, three_d, cells, step, by_block, source, times, interface, readers, shown, extra = render(example, choice, out)
    base = example["id"]
    file_name = os.path.basename(choice.path.replace("\\", "/"))
    caption = KIND_CAPTION[choice.kind].format(label=choice.label)
    if choice.kind == "mesh" and by_block:
        caption = "Mesh, colored by element block"
    if shown and shown[1]:
        bulk, hidden = shown
        caption += f", on the bulk phase{'s' if len(bulk) > 1 else ''} {', '.join(map(str, bulk))} only; " \
                   f"phase{'s' if len(hidden) > 1 else ''} {', '.join(map(str, hidden))} hidden " \
                   f"(no bulk physics for {choice.label} there)"
    if interface == "outline":
        caption += "; black lines: the boundaries of the bulk phases (interfaces and outer edges)"
    elif interface:
        caption += f"; black line: the level-set interface, {interface[0]} = {interface[1]:g}"
    if extra:
        caption = (caption + ". " + extra).rstrip(".")    # the page adds the final period itself
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
            "%Y-%m-%d", time.localtime(os.path.getmtime(choice.path))),
                   "run": "showcase" if example.get("showcase") else "regression"},
        "reproduce": example.get("reproduce", []),
    }
    se.hide_overlays(view)
    pv.SaveScreenshot(os.path.join(out, meta["thumb_raw"]), view, ImageResolution=SIZE)
    # the animation runs after the scene and images are written: its fixed color range changes the LUT
    if not override.get("no_animation"):
        animation = export_frames(view, source, choice, times, base, out, example["files"], readers)
        if animation:
            meta["animation"] = dict(animation, file=f"{base}-anim.webp")
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
            traceback.print_exc(file=sys.stdout)    # into export.log
            frames = traceback.extract_tb(exc.__traceback__)
            where = next((f for f in reversed(frames) if f.filename.endswith(("export_scenes.py", "gallery_select.py"))),
                         frames[-1])
            meta, reason = None, f"error: {exc} ({os.path.basename(where.filename)}:{where.lineno})"
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
    # run from the imported module, not from __main__: ParaView executes ProgrammableFilter scripts in
    # __main__'s namespace and star-imports vtk's numpy algorithms there, which replaces max/min/round
    import export_scenes
    sys.exit(export_scenes.main())
