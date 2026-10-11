# Level-set box beam with restarted adaptive refinement

A two-phase level-set topology optimization of a 2D cantilever box beam with a heat-method
regularization, run in two stages. The first deck optimizes with adaptive remeshing and writes
the ADV vector and the HMR refinement pattern to HDF5 every iteration; the second deck restarts
the optimization from the iteration-11 files. The example tests OPT restart from an ADV file
together with HMR restart from a stored refinement pattern. See also
{doc}`LevelSet_Boxbeam_Adaptive_Refinement <LevelSet_Boxbeam_Adaptive_Refinement>`.

```{exa-summary} projects/EXA/optimization/LevelSet_Boxbeam_Adaptive_Refinement_Restart
```

## Problem

| | |
|---|---|
| Domain | 3.0 × 1.0 rectangle (`tIs3D = false`; a 3D branch exists but is off) |
| Phases | 0 void, 1 design interior, 2 frame |
| Material | linear elastic, isotropic, E = 1.0, ν = 0.0, in UX, UY; bedding 10⁻⁶ in the bulk IWGs |
| Supports | displacement fixed to (0, 0) on side set 4 of the frame, symmetric Nitsche |
| Load | traction (0, 1.0) on side set 2 of the frame, applied only where y < 0.2 (`Func_Neumann_U`) |
| Auxiliary fields | THETA and PHID, two diffusion problems of the heat method on the interior and void phases |
| Objective | (SE_frame + SE_interior)/4.40520 + 0.05 · perimeter/14.0771 + 0.01 · heat-method penalty/0.243945 |
| Constraint | interior volume / 1.0 − 1 ≤ 0 |
| Optimizer | GCMMA, step size 0.05, penalty 5.0 |

The frame is a superellipse (exponent 24) centred at (1.5, 0.5) with a wall thickness of 0.05.
The design region is seeded by `Hole_Pattern_2D3D`, a periodic array of superellipse holes
(exponent 6). The normalization 4.40520 is the sum of the two initial strain energies
2.04222 + 2.36298 hard-coded in the deck.

## The two decks

| | `Levelset_Boxbeam_Create_File` | `Levelset_Boxbeam_Restart` |
|---|---|---|
| Hole array (x × y) | 3 × 2 | 7 × 3 |
| Hole radius | 0.4771 / 2 | 0.2771 / 3 |
| GCMMA `max_its` | 12 | 5 |
| OPT `restart_file` | `""` | `ADV_Alg_0_Iter_11.hdf5` |
| GCMMA `restart_index` | 0 | 10 |
| `reinitialize_interface_iter` | 10 | 50 |
| HMR | `write_refinement_pattern_file = true` | `restart_refinement_pattern_file = HMR_Background_Refinement_Iter_11.hdf5`, Lagrange mesh written to `HMRLagrangeMesh.exo` |

All other parameter lists are the same in both files. The two ctest cases run in declaration
order in one executable; the restart case depends on the HDF5 files left by the first.

## Setup in MORIS

**Mesh (HMR).** 30 × 10 elements, no initial refinement (`pattern_initial_refinement = 0,0`).
Two linear Lagrange meshes on patterns 0 and 1 and two linear B-spline meshes on the same
patterns; `lagrange_to_bspline = "0,1;-1"`. Refinement and staircase buffers are 1. In MSI the
displacements UX, UY use B-spline mesh 0 and the heat-method fields THETA, PHID use mesh 1
([coupling C6](../../../user-guide/interlinks.md#c6-mesh-indices-hmr-xtk-gen)).

**Adaptive refinement (MORIS general).** A remeshing parameter list in `ab_initio` mode remeshes
on the GEN fields `Box` and `Level_Set_Field` with `remeshing_levels_of_refinement = "2,1;2,1"`
on `remeshing_refinement_pattern = "0,1;0,1"`. In addition, each GEN geometry requests
`number_of_refinements = 1,0` with `refinement_mesh_index = 0,1`.

**Geometry and phases (GEN).** Two user-defined level-set geometries, `Box` (`Box_2D3D`) and
`Level_Set_Field` (`Hole_Pattern_2D3D`). The hole field is discretized on B-spline mesh 1
(`discretization_mesh_index = 1`), whose coefficients are the ADVs, bounded to ± the hole
radius. The phase table `[0, 1, 2, 2]` maps the geometry signs to void, interior and frame. A
`SCALED_FIELD` property `LevelsetField` (factor 1.0) exposes the hole field to FEM as the PDV
`LS1` on all bulk sets
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types` lists, in order, `IQIBulkStrainEnergy_Frame`, `IQIBulkStrainEnergy_Interior`,
`IQIBulkVolume_Interior`, `IQIPerimeter_InterfaceVoid`, `IQIHeatMethodPenalty`; the OPT
callbacks index `aCriteria` 0–4 in this order
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on, `high_to_low_dbl_side_sets = true`, Exodus output of the integration mesh.

**Physics (FEM).** Structure: `STRUC_LIN_ISO` constitutive models for frame and interior,
`STRUC_LINEAR_BULK` on both phases, `STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` support
(`DIRICHLET_NITSCHE`, penalty 100), `STRUC_LINEAR_NEUMANN` traction,
`STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` on the frame–interior set `dbl_iside_p0_1_p1_2`
(`NITSCHE_INTERFACE`, penalty 100), and `GHOST_NORMAL_FIELD` on `ghost_p1` and `ghost_p2`
(`GHOST_DISPL`, 0.005). Heat method: `DIFF_LIN_ISO` models for THETA (conductivity 1.0, capacity
0.1) and PHID (capacity 0.0, eigenstrain 1.0), `SPATIALDIFF_BULK` and `TIME_CONTINUITY_DOF` IWGs
on the interior and void sets, and `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` prescribing
THETA = 1.0 and PHID = 0.0 on the material–void interfaces (penalty 100). The regularization
criterion is the `HEAT_METHOD_PENALTY` IQI on PHID, with reference properties built from `LS1`
(`tLevelSetFunc`, `tLevelSetGradxFunc` and their `dv_derivative_functions`) and
`function_parameters = "1.6 / 1.0 / 36.0 / 0.2 / 0.2 / 0.333 / 0.333"`. The perimeter IQI is
`VOLUME` on the interface sets `iside_b0_2_b1_0,iside_b0_1_b1_0`; `H1_ERROR` and `PROPERTY`
(`IQILevelSet`) IQIs are defined for output.

**Solver (SOL).** `Amesos_Umfpack` direct solve. Three Newton solvers (one iteration each) for
THETA, PHID and UX,UY are staged as sub-solvers of an `NLBGS_SOLVER`, in that order; one time
solver over `THETA;PHID;UX,UY`.

**Optimization (OPT).** `user_defined` problem with GCMMA; objective, constraint and their
`dcriteria` derivatives are hand-written callbacks (the `dadv` derivatives are zero). The
restart keys are listed in the table above.

**Output (VIS).** `Levelset_Boxbeam_Restart.exo` (and the equivalent for the first deck),
standard mesh on all three phases, nodal UX, UY, PHID, THETA and LEVELSET, written every
iteration with `Time_Offset = 10.0`
([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)). GEN writes
`GEN_<name>.exo`.

## What the test checks

- **Create File case:** deletes `*exo*` and `*hdf5*` files in the working directory, runs
  `Levelset_Boxbeam_Create_File.so`, and requires a return code of 0.
- **Restart case:** runs `Levelset_Boxbeam_Restart.so`, requires a return code of 0, then opens
  `Levelset_Boxbeam_Restart.exo.e-s.0016` and compares nodal field index 2 at nodes 0, 100, …,
  800 against nine reference values (−0.000135399, 7.40848, 14.3591, 15.4603, 0.861916,
  −0.499474, −7.59165, −13.1653, −15.2731) with relative tolerance 10⁻⁴.

Both cases run on 1 processor.

## Run it

```bash
ctest -R Levelset_Boxbeam_Adaptive_Refinement_Restart-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name keeps the lowercase `Levelset` spelling of the legacy test.

## Notes

- `CMakeLists.txt` notes that the legacy test was also gated on `MORIS_HAVE_GCMMA`;
  `moris_add_example` registers it under `MORIS_HAVE_PARALLEL_TESTS` only.
- `tLoadLimitY` is computed in both decks but not used; `Func_Neumann_U` uses the literal 0.2.
