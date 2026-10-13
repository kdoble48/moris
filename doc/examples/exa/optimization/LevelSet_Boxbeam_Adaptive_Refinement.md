# Level-set box beam with adaptive refinement

A two-phase level-set topology optimization of a 2D cantilever box beam in which the
background mesh is refined once around both level-set geometries and a remeshing parameter
list (`ab_initio` mode) names those fields. A fixed frame surrounds a design region seeded
with a 6 × 3 hole pattern; the objective is stiffness plus a perimeter penalty under a mass
constraint. The deck demonstrates geometry-driven HMR refinement (`number_of_refinements` in
GEN); the .dox lists the feature as "Adaptive Refinement". Two decks run the same problem with
linear and quadratic interpolation.

```{exa-summary} projects/EXA/optimization/LevelSet_Boxbeam_Adaptive_Refinement
```

## Problem

| | |
|---|---|
| Domain | 3.0 × 1.0 rectangle |
| Phases | 0 void, 1 design interior, 2 frame |
| Material | linear elastic, isotropic: E = 1.0, ν = 0.0, density 1.0, plane problem in UX, UY |
| Bedding | 10⁻⁶ on frame and interior |
| Supports | left edge of the frame (side set 4), zero displacement, symmetric Nitsche |
| Load | traction (0, 1.0) on the right edge of the frame (side set 2), applied only below y = 0.2 |
| Objective | (SE_frame + SE_interior) / 9.123573 + 0.2 × perimeter / 41.5976 |
| Constraint | volume of the interior phase / 1.0 − 1 ≤ 0 |
| Optimizer | GCMMA, step size 0.05, penalty 5.0, 14 iterations |

The strain-energy normalization is written as 4.306294 + 4.817279. The frame is the
user-defined function `Box_2D3D`: a superellipse-type box centred at (1.5, 0.5) with
semi-axes 1.45 × 0.45 (wall thickness 0.05), exponent 24, scaled by √(1.45 · 0.45). The holes
come from `Hole_Pattern_2D3D`: superellipses of exponent 6 and radius 0.2471 / 3 ≈ 0.0824 on a
grid of spacing 0.5 in x and 1/3 in y (6 × 3 holes in the domain), combined by a maximum over
all holes; the loop extends two holes beyond the domain on each side.

## Setup in MORIS

**Mesh (HMR).** 60 × 20 elements, no initial refinement, refinement and staircase buffers of 1.
One Lagrange mesh and one B-spline mesh on pattern 0, of order 1 (linear deck) or 2 (quadratic
deck). The B-spline mesh carries both the level-set ADVs and the displacement DOFs (`UX`, `UY`
set to 0 in MSI). The Lagrange mesh is written to `HMRLagrangeMesh.exo`
(`HMRLagrangeMesh_Quadratic.exo` for the quadratic deck).

**Geometry and phases (GEN).** Two `USER_DEFINED` level sets, `Box` (`Box_2D3D`) and
`Level_Set_Field` (`Hole_Pattern_2D3D`), each with `number_of_refinements = 1` and
`refinement_mesh_index = 0`. `Level_Set_Field` is discretized on B-spline mesh 0 with bounds
±0.0824 (`tBsplineLimit = tHoleRadius`); its coefficients are the ADVs. The phase table
`[0, 1, 2, 2]` maps the geometry signs to void, interior, and frame. `IQI_types` is
`IQIBulkStrainEnergy_Frame`, `IQIBulkStrainEnergy_Interior`, `IQIBulkVolume_Interior`,
`IQIPerimeter_InterfaceVoid`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, Exodus output of the integration mesh.

**Physics (FEM).** The same IWG and IQI set as {doc}`LevelSet_Boxbeam <LevelSet_Boxbeam>`: two
`STRUC_LIN_ISO` models with `STRUC_LINEAR_BULK` IWGs and bedding;
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` support (`DIRICHLET_NITSCHE` 100.0);
`STRUC_LINEAR_NEUMANN` traction through `Func_Neumann_U`;
`STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` on `dbl_iside_p0_1_p1_2` (`NITSCHE_INTERFACE`
100.0); `GHOST_NORMAL_FIELD` on `ghost_p1` and `ghost_p2` with `ghost_order` equal to the
displacement order. Criteria: `STRAIN_ENERGY` and `VOLUME` per phase, perimeter as a `VOLUME`
IQI on `iside_b0_2_b1_0` and `iside_b0_1_b1_0`.

**Solver (SOL).** Direct `AMESOS_IMPL` solve, no preconditioner, one Newton iteration
(`NLA_rel_res_norm_drop = 1.0`).

**Optimization (OPT).** GCMMA with `step_size = 0.05`, `penalty = 5.0`, `max_its = 14`,
`restart_index = 0`.

**Remeshing (MORIS general).** A remeshing parameter list with `mode = "ab_initio"`,
`remeshing_field_names = "Box,Level_Set_Field"`, `remeshing_levels_of_refinement = "1"`, and
`remeshing_refinement_pattern = "0"`.

**Output (VIS).** `Levelset_Boxbeam_Adaptive_Refinement[_Quadratic].exo` on the frame and
interior sets with nodal UX, UY and the four criteria as global fields, written every
iteration with a time offset of 10.0.

## The two decks

| | `Levelset_Boxbeam_Adaptive_Refinement` | `..._Quadratic` |
|---|---|---|
| Lagrange and B-spline order | 1 | 2 |
| `GHOST_DISPL` parameter (both phases) | 0.005 | 1.0 |
| GEN `isocontour_tolerance`, `intersection_tolerance` | default | 10⁻⁸ on both geometries |
| HMR mesh output file | `HMRLagrangeMesh.exo` | `HMRLagrangeMesh_Quadratic.exo` |

## What the test checks

Two Catch2 test cases each run one deck, require a return code of 0, and then open the HMR
Lagrange mesh file and require exact matches of dimension (2), node count, and element count:

| Deck | 1 proc (nodes / elements) | 4 procs, ranks 0–1 | 4 procs, ranks 2–3 |
|---|---|---|---|
| linear | 4495 / 4182 | 1276 / 1157 | 1267 / 1148 |
| quadratic | 17545 / 4236 | 5553 / 1320 | 5553 / 1320 |

Any other processor count raises an error. Objective values are not compared. The first test
case removes `*exo*` and `*hdf5*` files before running.

## Run it

```bash
ctest -R Levelset_Boxbeam_Adaptive_Refinement- --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest names are `Levelset_Boxbeam_Adaptive_Refinement-1-procs` and
`Levelset_Boxbeam_Adaptive_Refinement-4-procs`; each runs both decks.

## Notes

- `tIs3D` is hard-coded to `false`; the 3D branches (45 × 15 × 15 mesh, `UZ`) are not
  exercised.
- `tLoadLimitY` is computed but `Func_Neumann_U` uses the literal limit 0.2.
- Both bulk IWGs are named `IWGBulkU_Frame`; the second one acts on the interior sets.
- {doc}`LevelSet_Boxbeam_Adaptive_Refinement_Restart <LevelSet_Boxbeam_Adaptive_Refinement_Restart>`
  is a related deck.
