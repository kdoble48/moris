# Mesh generation from an image (bear)

A mesh-generation-only run in which the geometry is a 2D image: a signed-distance field
sampled on an 880 × 528 pixel grid and stored in HDF5. MORIS is invoked with `--meshgen`
and an XML parameter file instead of a compiled deck; it builds a B-spline background mesh,
cuts it with the image level set, and writes the foreground mesh. The example demonstrates
the `image_file` geometry type of the mesh-generation workflow. Its siblings use an
{doc}`OBJ surface <obj_sphere>` and {doc}`pre-defined lines <pre_defined_diamond>`.

```{exa-summary} projects/EXA/mesh_generation/image_bear
```

## Problem

| | |
|---|---|
| Input files | `share/doc/mesh_generation/examples/Bear_Example.xml` and `Bear.hdf5` |
| Geometry | `field_type="image_file"`, `Bear.hdf5` (dataset `SDF` with 464 640 values, `Dimensions` = 880, 528) |
| Image extent | 8.80 × 5.28, origin (0, 0) — one pixel per 0.01 |
| Background grid | 220 × 132 elements on 8.80 × 5.28, origin (0, 0) — element edge 0.04 (4 pixels) |
| Physics | none |

The directory `projects/EXA/mesh_generation/image_bear` contains only `CMakeLists.txt` and
`example_test_case.cpp`; the XML and HDF5 inputs live in the mesh-generation manual's
example folder.

## Setup in MORIS

The mesh-generation XML has three sections. The manual in `share/doc/mesh_generation`
(`Overview.tex`, `Tutorial.tex`) describes them as inputs to HMR, GEN and XTK respectively.

**Mesh (HMR) — `BackgroundMeshes`.** `BaseGrid` of size 220 × 132; one `MeshGrid` (index 0)
with no initial and no interface refinement; one B-spline mesh (index 0) of polynomial order
1 on that grid.

**Geometry and phases (GEN) — `Geometries`.** A single image geometry with
`ImageDimensions` 8.80, 5.28 and `ImageOrigin` 0.0, 0.0. No `PhaseMap` is given.

**Cut and enrichment (XTK) — `ForegroundMesh`.** `DecompositionGrid` 0,
`FgPolynomialOrder` 1, `TriangulateAllFgElems` false, `ExtractionOperatorFormat` Elemental.
Per the manual, the foreground mesh is written as `foreground_mesh.exo` and the extraction
operators as one `Elemental_Extraction_Operators_B<i>.hdf5` file per B-spline mesh.

## What the test checks

The test copies `Bear_Example.xml` and `Bear.hdf5` from
`<MORIS root>/share/doc/mesh_generation/examples/` into the working directory and, only when
running on one processor, calls the workflow with `--meshgen Bear_Example.xml`. No result is
checked: the return code is ignored and the output check is a `TODO` comment.

## Run it

```bash
ctest -R Image_Bear-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The CMake registration uses `NO_INPUTS`: no deck `.so` is built for this example. Outside
the test, the same run is `moris --meshgen Bear_Example.xml` in a directory containing both
input files.
