# Grain structure

A static heat-conduction problem on a 2D square partitioned into several grains by four
user-defined indicator fields. Each field returns ±1 rather than a distance, and the regions are
made mutually exclusive in the field definitions, so XTK cuts the background mesh into
multiple material phases. Temperature continuity across the grain boundaries is enforced by a
Nitsche interface condition. The test checks only the size of the resulting mesh.

```{exa-summary} projects/EXA/thermal/diffusion/Grain_Structure
```

## Problem

| | |
|---|---|
| Domain | 4.0 × 4.0 square, offset (−2, −2) |
| Grain 1 (`Colors_1`) | x > 0.1117 and y ≤ −0.11 |
| Grain 2 (`Colors_2`) | inside the circle of radius 1.2 about (2, 2), excluding grain 1 |
| Grain 3 (`Colors_3`) | x − y + 1.5 ≤ 0, excluding grains 1 and 2 |
| Grain 4 (`Colors_4`) | inside the circle of radius 2.6 about (−2, −2), excluding grains 1–3 |
| Material | density 1.0, heat capacity 1.0, conductivity 0.00005 in every grain |
| Right edge | temperature 2.0, Nitsche-enforced, on phases 0, 4, 8 |
| Left edge | heat flux 2.0 on phases 1, 2 |
| Grain boundaries | temperature continuity, symmetric Nitsche |
| Analysis | static, one solve |

## Setup in MORIS

**Mesh (HMR).** 11 × 11 elements, linear Lagrange and B-spline meshes, refinement and staircase
buffers 1.

**Geometry and phases (GEN).** Four `USER_DEFINED` geometries (`Colors_1` … `Colors_4`), each
returning +1 inside its region and −1 elsewhere, refined twice on mesh 0, not discretized
(`discretization_mesh_index = −1`), with `isocontour_tolerance = 1e-12`. No phase function is
given; the deck uses the default sign-pattern phases 0, 1, 2, 4, 8. GEN output: `gen_output.exo`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
off (`isGhost = false`), integration-mesh Exodus output, `write_cluster_measures_to_exo = false`.

**Physics (FEM).** Two `DIFF_LIN_ISO` models with identical properties (`CMDiffusion`,
`CMDiffusion_2`). `SPATIALDIFF_BULK` on the bulk sets of phases 0, 1, 2, 4, 8;
`SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 100.0) on `SideSet_2` of phases
0, 4, 8; `SPATIALDIFF_NEUMANN` on `SideSet_4` of phases 1, 2; and
`SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, 100.0) on six double-sided
interface sets: `dbl_iside_p0_0_p1_{1,2,4,8}` and `dbl_iside_p0_1_p1_{2,8}`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
A ghost IWG on `ghost_p160` is defined inside `if (isGhost)` and is therefore inactive. One
`DOF` IQI for TEMP.

**Solver (SOL).** `BELOS_IMPL` with preconditioner 0, an `IFPACK` ILU preconditioner; Newton with
relative drop 10⁻⁴, relaxation 1.0, and one iteration; one static step. MSI sets
`order_adofs_by_host = false`.

**Output (VIS).** `Grain_Structure.exo` on the bulk sets with nodal TEMP and `IQIBulkTEMP`.

## What the test checks

The single Catch2 case runs the deck (`gInterpolationOrder = 1`, `gDim = 2`) and checks that
`Grain_Structure.exo` has 2 dimensions, 2451 nodes, and 1875 elements. It compares no field
values and does not check the workflow return code.

## Run it

```bash
ctest -R Grain_Structure-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The directory has no `example_documentation.dox`.
- The Catch2 tags are `[structure],[linear]` although the example lives under
  `thermal/diffusion`.
