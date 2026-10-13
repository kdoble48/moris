# Mesh generation from an OBJ surface (sphere)

A 3D mesh-generation-only run in which the geometry is a triangulated surface read from a
Wavefront OBJ file. MORIS is invoked with `--meshgen` and an XML parameter file instead of a
compiled deck; it builds a B-spline background mesh, evaluates a signed-distance field from
the surface, cuts the background mesh, and writes the foreground mesh. The example
demonstrates the `object_file` geometry type. Its siblings use an
{doc}`image <image_bear>` and {doc}`pre-defined lines <pre_defined_diamond>`.

```{exa-summary} projects/EXA/mesh_generation/obj_sphere
```

## Problem

| | |
|---|---|
| Input files | `OBJ_Sphere.xml` and `Sphere.obj` in this directory |
| Geometry | `Sphere.obj`: 1202 vertices, 2400 triangles (written by the Visualization Toolkit) |
| Surface extent | x ∈ [−60.53, 17.90], y ∈ [−49.62, 28.96], z ∈ [−39.37, 39.37] (about 78.7 in each direction) |
| Background grid | 17 × 17 × 17 elements on 80 × 80 × 80, origin (−60, −50, −40) |
| Physics | none |

The surface nearly fills the background box; its minimum x (−60.53) lies 0.53 outside the
grid origin.

## Setup in MORIS

The mesh-generation XML has three sections. The manual in `share/doc/mesh_generation`
(`Overview.tex`, `Tutorial.tex`) describes them as inputs to HMR, GEN and XTK respectively.

**Mesh (HMR) — `BackgroundMeshes`.** `BaseGrid` of size 17 × 17 × 17; one `MeshGrid`
(index 0) with no initial and no interface refinement; one B-spline mesh (index 0) of
polynomial order 1.

**Geometry and phases (GEN) — `Geometries`.** A single geometry with
`field_type="object_file"`, `FileName` `Sphere.obj`, `ObjectOrigin` 0.0, 0.0, 0.0 and
`SdfShift` 0.00. No `PhaseMap` is given.

**Cut and enrichment (XTK) — `ForegroundMesh`.** `DecompositionGrid` 0,
`FgPolynomialOrder` 1, `TriangulateAllFgElems` false, and `OutputExtractionOperators` false.
Per the manual, the foreground mesh is written as `foreground_mesh.exo`; with
`OutputExtractionOperators` false the extraction-operator HDF5 output is suppressed.

## What the test checks

The test copies `OBJ_Sphere.xml` and `Sphere.obj` from
`<MORIS root>/projects/EXA/mesh_generation/obj_sphere/` into the working directory and, only
when running on one processor, calls the workflow with `--meshgen OBJ_Sphere.xml`. No result
is checked: the return code is ignored and the output check is a `TODO` comment.

## Run it

```bash
ctest -R OBJ_Sphere-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The CMake registration uses `NO_INPUTS`: no deck `.so` is built for this example.

## Notes

- A comment at the top of `OBJ_Sphere.xml` reads "this example has a large mesh. Memory
  consumption ~10GB".
