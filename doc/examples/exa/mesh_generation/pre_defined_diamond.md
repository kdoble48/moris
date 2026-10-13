# Mesh generation with pre-defined lines (rotated square)

A 2D mesh-generation-only run in which a diamond (a square rotated by 45°) is described by
four pre-defined line level sets. MORIS is invoked with `--meshgen` and an XML parameter file
instead of a compiled deck. The example demonstrates the `pre_defined` geometry type, a
`PhaseMap` that merges the sign combinations of several geometries into two materials, two
B-spline background meshes of different order, and interface refinement. Its siblings use an
{doc}`image <image_bear>` and an {doc}`OBJ surface <obj_sphere>`.

```{exa-summary} projects/EXA/mesh_generation/pre_defined_diamond
```

## Problem

| | |
|---|---|
| Input file | `share/doc/mesh_generation/examples/Rotated_Square_Example.xml` |
| Geometry | four lines: through (−0.70710678118, 0) with normals (−1, 1) and (−1, −1); through (0.70710678118, 0) with normals (1, 1) and (1, −1) |
| Phase map | `1,0;2,0;3,0;4,0;5,0;8,0;10,0;12,0;0,1` |
| Background grid | 8 × 8 elements on 2.0 × 2.0, origin (−1, −1) |
| Physics | none |

The four lines (y = ±(x + 0.7071) and y = ±(x − 0.7071)) bound a square with corners at
(±0.7071, 0) and (0, ±0.7071). The mesh-generation manual (`Tutorial.tex`) defines the `PhaseMap` syntax as
semicolon-separated pairs "phase index, material index": here phase index 0 is mapped to
material 1 and the listed phase indices 1, 2, 3, 4, 5, 8, 10, 12 to material 0.

## Setup in MORIS

The mesh-generation XML has three sections. The manual in `share/doc/mesh_generation`
(`Overview.tex`, `Tutorial.tex`) describes them as inputs to HMR, GEN and XTK respectively.

**Mesh (HMR) — `BackgroundMeshes`.** `BaseGrid` of size 8 × 8. Two mesh grids:

| Mesh grid | Initial refinements | Interface refinements | B-spline mesh on it |
|---|---|---|---|
| 0 | 1 | 1 | B-spline mesh 0, order 2 |
| 1 | 0 | 0 | B-spline mesh 1, order 1 |

**Geometry and phases (GEN) — `Geometries`.** Four `field_type="pre_defined"`,
`geom="line"` geometries, each given by a `Point` and a `Normal`, followed by the `PhaseMap`
above.

**Cut and enrichment (XTK) — `ForegroundMesh`.** `DecompositionGrid` 0 (the refined grid),
`FgPolynomialOrder` 2, `TriangulateAllFgElems` false, `ExtractionOperatorFormat` Elemental.
Per the manual, the foreground mesh is written as `foreground_mesh.exo` and the extraction
operators as `Elemental_Extraction_Operators_B0.hdf5` and `…_B1.hdf5`, one per B-spline mesh.
The manual uses this example to illustrate the output block names (`HMR_dummy_c_p0`,
`HMR_dummy_n_p0`, …).

## What the test checks

The test copies `Rotated_Square_Example.xml` from
`<MORIS root>/share/doc/mesh_generation/examples/` into the working directory and, only when
running on one processor, calls the workflow with `--meshgen Rotated_Square_Example.xml`. No
result is checked: the return code is ignored and the output check is a `TODO` comment.

## Run it

```bash
ctest -R Pre_Defined_Diamond-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The CMake registration uses `NO_INPUTS`: no deck `.so` is built for this example.
