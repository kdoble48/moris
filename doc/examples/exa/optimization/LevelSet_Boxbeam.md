# Level-set box beam

A two-phase level-set topology optimization of a 2D cantilever box beam: a fixed outer frame
surrounds a design region whose material layout is optimized for stiffness under a mass
constraint. This is the legacy-format deck — per-module parameter-list functions with
hand-written OPT callbacks — of the problem that
{doc}`LevelSet_Boxbeam_NewIO <LevelSet_Boxbeam_NewIO>` writes with Deck API v2.

```{exa-summary} projects/EXA/optimization/LevelSet_Boxbeam
```

## Problem

| | |
|---|---|
| Domain | 3.0 × 1.0 rectangle |
| Phases | 0 void, 1 design interior, 2 frame |
| Material | linear elastic, isotropic: E = 1.0, ν = 0.0, density 1.0, plane problem in UX, UY |
| Bedding | 10⁻⁶ on frame and interior |
| Supports | left edge of the frame (side set 4), zero displacement, symmetric Nitsche |
| Load | traction (0, 1.0) on the right edge of the frame (side set 2), applied only below y = `tLoadLimitY` (≈ 0.2) |
| Objective | (SE_frame + SE_interior) / 3.26439 + 0.2 × perimeter / 29.0669 |
| Constraint | volume of the interior phase / 1.0 − 1 ≤ 0 |
| Optimizer | GCMMA, step size 0.01, penalty 5.0, 2 iterations |

`tLoadLimitY` is `floor(0.2 / h) · h` with element edge h = 1/15 / 2² = 1/60. The frame is a
superellipse centred at (1.5, 0.5) with semi-diameters 1.45 × 0.45 (wall thickness 0.05) and
exponent 24; the design region is seeded with a 15 × 5 field array of superellipse holes
(semi-diameter 0.065, exponent 8) spanning x ∈ [0.12, 2.88], y ∈ [0.12, 0.88].

## Setup in MORIS

**Mesh (HMR).** 45 × 15 elements with two levels of initial refinement (element edge 1/60),
refinement and staircase buffers of 3. One linear Lagrange mesh and two linear B-spline
meshes: mesh 0 carries the level-set ADVs, mesh 1 the displacement DOFs (`UX`, `UY` set to 1 in
MSI) ([coupling C6](../../../user-guide/interlinks.md#c6-mesh-indices-hmr-xtk-gen)).

**Geometry and phases (GEN).** A `SUPERELLIPSE` level set for the frame and a `SUPERELLIPSE`
field array for the holes. The array is discretized on B-spline mesh 0
(`discretization_mesh_index = 0`) with bounds ±0.06; its coefficients are the ADVs. The phase
table `[0, 1, 2, 2]` maps the two geometry signs to void, interior, and frame. `IQI_types` is
`IQIBulkStrainEnergy_Frame`, `IQIBulkStrainEnergy_Interior`, `IQIBulkVolume_Interior`,
`IQIPerimeter_InterfaceVoid`; the OPT callbacks index `aCriteria` 0–3 in that order
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on, Exodus output of the integration mesh.

**Physics (FEM).** Two `STRUC_LIN_ISO` models (frame and interior, same properties), each with
a `STRUC_LINEAR_BULK` IWG carrying the bedding property. The support uses
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` with a `DIRICHLET_NITSCHE` penalty of 100.0. The
traction uses `STRUC_LINEAR_NEUMANN` with the user function `Func_Neumann_U`. The frame–interior
interface `dbl_iside_p0_1_p1_2` uses `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` with a
`NITSCHE_INTERFACE` penalty of 100.0. Ghost terms use `GHOST_NORMAL_FIELD` with `GHOST_DISPL`
0.01 on `ghost_p1` and `ghost_p2`. Criteria: `STRAIN_ENERGY` per phase, `VOLUME` per phase
weighted by density, and the perimeter as a `VOLUME` IQI on the material–void interfaces
`iside_b0_2_b1_0` and `iside_b0_1_b1_0`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** Direct `AMESOS_IMPL` solve with the default Amesos solver, no preconditioner;
one Newton iteration (`NLA_rel_res_norm_drop = 1.0`); one time step over a time frame of 1.0.

**Optimization (OPT).** GCMMA with `step_size = 0.01`, `penalty = 5.0`, `max_its = 2`, on the
user-defined objective and constraint callbacks above. The normalization constants
(3.26439, 29.0669) are hard-coded.

**Output (VIS).** `Levelset_Boxbeam.exo` with nodal UX and UY on the frame and interior sets,
written every iteration with a time offset of 10.0; GEN writes `GEN_Levelset_Boxbeam.exo`.

## Second deck: `Levelset_Boxbeam_Lin_Quad_LevelSet.cpp`

The directory also contains `Levelset_Boxbeam_Lin_Quad_LevelSet.cpp`, which CMake does not
build and whose test case is commented out. It uses the older
`Vector< Vector< ParameterList > >` interface, a quadratic B-spline mesh for the level set and
a linear one for displacements (initial refinement 2 and 3), user-defined frame
(`Box_2D3D`) and hole-pattern (`Hole_Pattern_2D3D`, 9 × 3 holes) geometries on a 15 × 5 base
mesh, and GCMMA with step size 0.05 and 3 iterations.

## What the test checks

The test runs `Levelset_Boxbeam.so` through the workflow interface and requires a return code
of 0. It checks that the deck loads, cuts, solves, and completes its GCMMA iterations; it does
not compare objective values or mesh sizes. (`check_results` exists but is not called, and it
returns before its `REQUIRE` statements.)

## Run it

```bash
ctest -R '^Levelset_Boxbeam-' --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest names are `Levelset_Boxbeam-1-procs` and `Levelset_Boxbeam-3-procs` (lowercase
"s"). The anchored pattern excludes the `Levelset_Boxbeam_Adaptive_Refinement` tests.

## Notes

- Both bulk IWGs are named `IWGBulkU_Frame`; the second one acts on the interior sets.
- `PropFlux` and `tFrameFreeSSets` are defined but unused; the latter concatenates two set
  lists without a separating comma.
- The CMake comment notes that the legacy test was additionally gated on `MORIS_HAVE_GCMMA`,
  and that `moris_add_example` registers it under `MORIS_HAVE_PARALLEL_TESTS` only.
- Related decks: {doc}`LevelSet_Boxbeam_Adaptive_Refinement <LevelSet_Boxbeam_Adaptive_Refinement>`
  (user-defined geometries with refinement) and
  {doc}`LevelSet_Boxbeam_NewIO <LevelSet_Boxbeam_NewIO>`.
