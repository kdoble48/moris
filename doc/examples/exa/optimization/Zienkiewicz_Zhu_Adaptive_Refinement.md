# Zienkiewicz–Zhu adaptive refinement in level-set topology optimization

A 2D level-set topology optimization of a stadium-shaped plate with two bolt holes, in which the
background mesh is adapted during the run. HMR remeshes on two criteria: a band around each
level-set interface, and an elemental Zienkiewicz–Zhu (ZZ) von Mises stress error field computed
by FEM. The design layout enters the stiffness through a density-like interpolation of the
level-set field, and a heat-method penalty regularizes the level set. The directory has no
`example_documentation.dox`.

```{exa-summary} projects/EXA/optimization/Zienkiewicz_Zhu_Adaptive_Refinement
```

## Problem

| | |
|---|---|
| Domain | 8.0 × 4.0 rectangle centred at the origin |
| Outer shape | stadium: two half-discs of radius 1.237 centred at (±2.0, 0) joined by straight edges at y = ±1.237 |
| Holes | at (±2.0, 0): inner radius 0.517, outer radius 0.7173 (two separate level sets) |
| Phases | phase table of 16 entries; FEM uses phase 9 (design interior), 8 (void) and 11 (frame) through raw set names |
| Material | linear elastic, ν = 0.0; frame E = 1.0; interior E = 1.0 · ρ³ with ρ derived from the level set; bedding 10⁻⁸ in the interior |
| Supports | Nitsche Dirichlet u = 0 on `iside_b0_11_b1_15`, restricted to x < 0 by `Func_Select` |
| Load | traction (0, 1.0) on the same side set, applied where x > 0 and y < 0 (`Func_Neumann_U`) |
| Objective | 0.1 · (SE_frame + SE_interior)/90 + mass + 0.1 · perimeter/8 + 0.2 · heat-method penalty/250 |
| Constraint | criterion 5, the maximum-stress IQI (one constraint, type 1) |
| Optimizer | GCMMA, step size 0.02, penalty 100, 19 iterations |

The interior density is ρ = ((φ/2.4) − 0.5) / (1 − 0.5), shifted by the density-shift value and
clipped to [10⁻⁴, 1]; 0.5 is `phi_sh_scale = 1.2 / (3 · 0.8)`.

## Setup in MORIS

**Mesh (HMR).** 60 × 30 elements on [−4, 4] × [−2, 2], no initial refinement. A quadratic
Lagrange mesh on pattern 0, a linear B-spline mesh on pattern 1 and a quadratic B-spline mesh on
pattern 2; `lagrange_to_bspline = "0,1"`. MSI places UX, UY, STRESS_DOF, THETA and PHID on
B-spline mesh 0.

**Adaptive remeshing (MORIS general parameters).** A remeshing parameter list in `advanced`
mode remeshes on five fields: `Level_Set_Field`, `Outer_shape`, `Inner_Holes`, `Outer_Holes`
and `FieldZienkiewiczZhu`. The first four use `refinement_function`, which refines elements
whose nodal values straddle the band 1.2 ± 0.2 up to level 2 (hold at level 2) and drops the
others. The ZZ field uses `refinement_function_stress`, which compares the elemental value to a
fixed reference 0.1: ratio ≥ 0.7 refines (up to level 2), ratio ≥ 0.2 holds, smaller values
coarsen. Old patterns 0, 1, 2 are copied to 3, 4, 5; `remeshing_refinement_pattern` is
`0;0;0;0;0,1`; `minimum_refinement_level` is
`"0, 0, 30, 10, 100, 0, 500; 1, 0, 30, 0, 500;2, 0, 500"`; remeshed meshes are written
(`output_meshes = true`). OPT sets `reinitialize_interface_iter = 18`.

**Geometry and phases (GEN).** Four `USER_DEFINED` level sets, each offset so that the
interface sits at `isocontour_threshold = 1.2`: `general_shape_outer` (`Outer_shape`),
`Hole_inner`, `Hole_outer`, and `Const_Geometry` (`Level_Set_Field`, initial value 2.0
everywhere). `Level_Set_Field` is discretized on B-spline mesh 1 with bounds [0.0, 2.4]; its
coefficients are the ADVs. A `SCALED_FIELD` property (scaling 1.0) maps it to the PDV `LS1` on
the phase-8 and phase-9 sets
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
The phase table maps the 16 sign combinations one to one except that combination 10 maps to
phase 11. `IQI_types`: `IQIBulkStrainEnergy_Frame`, `IQIBulkStrainEnergy`, `IQIBulkVolume`,
`IQIPerimeter_InterfaceVoid`, `IQIHeatMethodPenalty`, `IQIMaxStress`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization, `high_to_low_dbl_side_sets`, Exodus output of the integration mesh, verbose
output on.

**Physics (FEM).** The deck calls `hack_for_legacy_fem()` and addresses sets by name
(`HMR_dummy_*_p9`, `iside_b0_9_b1_8`, `dbl_iside_p0_9_p1_11`, `ghost_p9`, …) rather than by
FEM phases ([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

- Elasticity: `STRUC_LIN_ISO` CMs for interior (E from `tYoungsFunc`, parameters 1.0, 3.0) and
  frame; `STRUC_LINEAR_BULK` in both; `STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE`
  (`DIRICHLET_NITSCHE`, 100) and `STRUC_LINEAR_NEUMANN` on the frame side set;
  `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, 100) between interior and
  frame.
- Stress: `STRUC_VON_MISES_STRESS` solves a von Mises stress field `STRESS_DOF` in the interior.
- Heat method: THETA and PHID diffusion (`DIFF_LIN_ISO`, `SPATIALDIFF_BULK`) on interior and
  void, Dirichlet values (`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE`, 100) on both sides of the
  interior–void interface, and `TIME_CONTINUITY_DOF` on THETA (weights 10.0).
- Ghost (`GHOST_NORMAL_FIELD`, order 1): interior displacement 0.005, frame displacement 0.001,
  `STRESS_DOF` 0.001, THETA and PHID 0.001.

Criteria: `STRAIN_ENERGY` on frame and interior, `VOLUME` of the interior weighted by the
density property, `VOLUME` on the interior–void side set (perimeter), `HEAT_METHOD_PENALTY` on
PHID (parameters 1.2 / φ-gradient / γ / 0.0 / 0.2 / 1.0 / 1.0), and `MAX_DOF` of `STRESS_DOF`
(reference 5.0, exponent 2.0, shift 1.0). The IQI `IQIBulkZienkiewiczZhu`
(`ZIENKIEWICZ_ZHU_VON_MISES_STRESS`) feeds the FEM field `FieldZienkiewiczZhu`
(`ELEMENTAL`, `FIELD_1`), which is the stress remeshing field.

**Solver (SOL).** One `AMESOS_IMPL` linear algorithm (no `Solver_Type` set). The time solver
uses an `NLBGS_SOLVER` (one pass, drop 10⁻⁷) over four Newton sub-solvers of one iteration each:
THETA, PHID, (UX, UY) with `STRESS_DOF` as secondary DOF, and `STRESS_DOF`. No separate
sensitivity solver is set.

**Optimization (OPT).** GCMMA (`step_size = 0.02`, `penalty = 100`, `max_its = 19`) on the
objective and constraint above. The density shift (`compute_density_shift`, called from
`Output_Criterion`) starts at output count 20 and ramps quadratically over 150 counts from 10⁻⁴
to 1.

**Output (VIS).** `Zienkiewicz_Zhu_Adaptive_Refinement.exo` on interior, void and frame sets:
nodal UX, UY, VOL, PHID, THETA, the level set, heat-method fields, truncated PHID, STRESS, and
the ZZ field as `ELEMENTAL_AVG`, time offset 10.

## What the test checks

The test runs the workflow and requires a return code of 0. On 1 processor it then opens
`Remeshing_Lagrange_Mesh_Pattern_1_Iter_19.exo` and requires 2 dimensions, 2475 nodes and 2295
elements. A 4-processor branch exists but its check is commented out, and the test calls
`MORIS_ERROR` for any other processor count. No field values or objective values are compared.

## Run it

```bash
ctest -R Zienkiewicz_Zhu_Adaptive_Refinement --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest `Zienkiewicz_Zhu_Adaptive_Refinement-1-procs` is registered only when
`MORIS_HAVE_GCMMA` is set; the input `.so` and test executable are built either way.

## Notes

- `tDerYoungsFunc` and `tDerDensityFunc` call `MORIS_ERROR( false, "Do not need this one" )`.
- With `max_its = 19` and a density-shift start of 20 output calls, the code sets the shift to
  0.0 until the counter reaches 20.
- `pdv_mesh_set_names` adds `ghost_p9` only in the branch where `tUseGhost` is false; with
  ghost enabled (the deck default) the PDV is defined on the phase-8 and phase-9 sets only.
