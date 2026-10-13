# Level-set beam with SIMP hole seeding

A 2D cantilever beam optimized with a combined level-set / SIMP formulation. The design field
starts uniform at 0.8 with the level-set threshold at 0.5, so the initial design contains no
holes; inside the material the same field is mapped to a SIMP density that scales Young's
modulus and the mass. A density-shift schedule and a heat-method regularization act on this
field. The example demonstrates nucleating holes in a level-set design through the SIMP
interpolation instead of a geometric hole seed.

```{exa-summary} projects/EXA/optimization/Level_Set_Beam_SIMP_Hole_Seeding
```

## Problem

| | |
|---|---|
| Domain | 60.0 × 20.0 rectangle |
| Phases | 0 void (field < 0.5), 1 material (field > 0.5); no `phase_table` is set |
| Material | linear elastic, ν = 0.0, E = 1.0 · ρ³ with the SIMP density ρ below; bedding 10⁻⁸ |
| Supports | displacement (0, 0) on `SideSet_2_n_p1`, symmetric Nitsche (penalty 100) |
| Load | traction (0, 10.0) on `SideSet_4_n_p1`, applied only where y < 0.2 (`Func_Neumann_U`) |
| Auxiliary fields | THETA and PHID, the two diffusion problems of the heat method |
| Objective | SE / 1437.01 + 0.02 · perimeter / 60 + 0.02 · heat-method penalty / 60 |
| Constraint | density-weighted volume / 400 − 1 ≤ 0 |
| Optimizer | GCMMA, step size 0.021, penalty 5.0, 20 iterations |

The SIMP density is computed from the level-set value φ (PDV `LS1`) as
ρ = s + (1 − s) · clamp((φ − 0.5)/(1 − 0.5)), where values below 0 are set to 10⁻⁴ and values
above 1 to 1 (`tYoungsFunc`, `tDensityFunc`). The density shift s is updated in
`Output_Criterion`, registered as the time solver's `TSA_Output_Criteria`, which increments an
internal counter n on each call: s = 0 for n < 20, s = 10⁻⁴ + (1 − 10⁻⁴) · ((n − 20)/100)² for 20 ≤ n < 120, and s = 1
afterwards.

## Setup in MORIS

**Mesh (HMR).** 60 × 20 elements on the 60 × 20 domain. Pattern 0 has one level of initial
refinement, pattern 1 none. Lagrange meshes: quadratic on pattern 0, linear on pattern 1.
B-spline meshes: linear on pattern 0, quadratic on pattern 1. MSI places UX, UY on B-spline
mesh 0.

**Geometry and phases (GEN).** One user-defined level-set geometry, `Level_Set_Field`
(`Const_Geometry`, value 0.8 everywhere) with `isocontour_threshold = 0.5`. Its coefficients on
B-spline mesh 1 are the ADVs (`discretization_mesh_index = 1`), bounded to [0, 1]. Two
`SCALED_FIELD` properties (factor 1.0) expose the field to FEM as PDVs `LS1` and `DENSITY` on
the bulk sets of both phases and on `ghost_p1`
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types`: `IQIBulkStrainEnergy`, `IQIBulkVolume`, `IQIPerimeter_InterfaceVoid`,
`IQIHeatMethodPenalty` — `aCriteria` 0–3 in the OPT callbacks
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on, `high_to_low_dbl_side_sets = true`, Exodus output of the integration mesh.

**Physics (FEM).** Structure on the material phase: `STRUC_LIN_ISO` with the SIMP-scaled
`PropYoungs` (`function_parameters = "1.0,3.0"`), `STRUC_LINEAR_BULK` with bedding,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 100), `STRUC_LINEAR_NEUMANN`,
and `GHOST_NORMAL_FIELD` on `ghost_p1` (`GHOST_DISPL`, 0.005). Heat method on both phases:
`DIFF_LIN_ISO` models for THETA (conductivity 100.0, capacity 0.001) and PHID (capacity 0.0,
eigenstrain 1.0), `SPATIALDIFF_BULK`, `TIME_CONTINUITY_DOF` for THETA, and
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` prescribing THETA = 1.0 and PHID = 0.0 on the
material–void interfaces `iside_b0_1_b1_0` and `iside_b0_0_b1_1` (penalty 100). Criteria: the
`STRAIN_ENERGY` IQI on the material, the `VOLUME` IQI weighted by `PropDensity` (ρ), the
perimeter as `VOLUME` on `iside_b0_1_b1_0`, and `HEAT_METHOD_PENALTY` on the void sets with
reference properties built from φ − 0.5 and parameters
`1.0 / 0.88222 / 4.69867 / 0.1 / 1.0 / 1.0 / 1.0` (B-spline limit 1.0; the second and third
values are computed from a bandwidth of 3 element lengths as ln(199)/6 and
2 ln 10 / (2/(e^(−2·3·0.88222) + 1) − 1)²). `H1_ERROR`, `DOF` and `PROPERTY` IQIs are defined
for output.

**Solver (SOL).** Amesos direct solve, `Amesos_Mumps` if MORIS is built with
`MORIS_USE_MUMPS`, otherwise `Amesos_Superludist`. Newton solvers (up to 5 iterations,
relative residual drop 10⁻⁷) for THETA, PHID and UX,UY are staged in that order as sub-solvers
of an `NLBGS_SOLVER`. The solver warehouse sets `SOL_save_operator_to_matlab = "Mat.dat"`.

**Optimization (OPT).** `user_defined` GCMMA problem with hand-written objective, constraint
and `dcriteria` callbacks (the `dadv` derivatives are zero);
`reinitialize_interface_iter = 10000`.

**Output (VIS).** `Level_Set_Beam_SIMP_Hole_Seeding.exo` on the material and void bulk sets
with nodal UX, UY, VOL, PHID, THETA, LVLSET, HEATMETHOD and LEVELSETHEAT, written every
iteration with `Time_Offset = 10.0`
([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)).

## What the test checks

The test runs `Level_Set_Beam_SIMP_Hole_Seeding.so`, requires a return code of 0, and opens
`Level_Set_Beam_SIMP_Hole_Seeding.exo.e-s.0018`. The check functions print the number of
dimensions, nodes and elements and then return unconditionally (`if ( true ) { …; return; }`),
so the reference counts in the file (serial: 10477 nodes, 8178 elements; per-rank values for 4
processors) are not compared. Any processor count other than 1 or 4 raises an error.

## Run it

```bash
ctest -R Level_Set_Beam_SIMP_Hole_Seeding-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Level_Set_Beam_SIMP_Hole_Seeding-4-procs --output-on-failure
```

The ctest entries are registered only when `MORIS_HAVE_GCMMA` is on; the `.so` and the test
executable are built either way.

## Notes

- The variable names for the boundary sets are swapped relative to their use: the Dirichlet
  IWG uses `tLoadSets` (`SideSet_2_n_p1`) and the traction IWG uses `tDirichletSets`
  (`SideSet_4_n_p1`). Each string lists the same set twice and only the `_n_` variant.
- `tDerYoungsFunc` and `tDerDensityFunc`, registered as `dv_derivative_functions`, abort with
  `MORIS_ERROR( false, "Do not need this one" )` if called.
- The remeshing parameter list names the fields `Box,Level_Set_Field` with zero levels of
  refinement; the deck defines no `Box` geometry.
