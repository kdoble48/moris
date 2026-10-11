# Shape sensitivity of a thermoelastic two-material problem

A finite-difference verification of shape sensitivities for a static thermoelastic problem with
two materials separated by a circular level set. Both phases carry material; temperature and
displacement are coupled across the material interface with Nitsche terms, and both phases are
ghost-stabilized. The OPT `SWEEP` algorithm compares the analytic gradients of the strain energy
and volume of one phase with respect to the circle's center and radius against forward,
backward, and central finite differences, and the test compares the final adjoint vector with a
stored reference. One deck solves monolithically, the other with a staggered (NLBGS) solver.
The single-material transient counterpart is
{doc}`Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient <Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient>`.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Two_Material_Sweep_Thermoelastic
```

## Problem

| | |
|---|---|
| Domain | 6.0 × 4.0 rectangle |
| Geometry | `CIRCLE` level set, center (3.0, 2.21), radius 1.4 |
| Phases | material 1 on phase 1, material 2 on phase 0 |
| Conductivity | 1.0 (material 1), 2.0 (material 2) |
| Structure | `STRUC_LIN_ISO`, `PLANE_STRESS`, E = 1.0, ν = 0.0 in both; CTE 1.0 (material 1), 2.0 (material 2); reference temperature 0.0 |
| Thermal BCs | TEMP = 0 on `SideSet_4_n_p1` (x-min side, symmetric Nitsche); surface flux 5.0 on `SideSet_2_n_p1` (x-max side) |
| Structural BCs | zero displacement on `SideSet_4_n_p1` and `SideSet_2_n_p1` (unsymmetric Nitsche) |
| Interface | `dbl_iside_p0_1_p1_0`, symmetric Nitsche for TEMP and for UX, UY |
| Time | one time step over a time frame of 1.0 |
| Objective | `IQIBulkStrainEnergyDISP`: strain energy of material 2 on phase 0 |
| Constraint | `IQIBulkVolume`: volume of phase 0 |

## Setup in MORIS

**Mesh (HMR).** 6 × 4 elements, linear Lagrange and linear B-spline meshes on pattern 0,
`pattern_initial_refinement` 1, buffers 1.

**Geometry and phases (GEN).** One `CIRCLE` field with three ADVs (`center_x`, `center_y`,
`radius`), bounds equal to the initial values. `IQI_types` is
`IQIBulkStrainEnergyDISP, IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on (`tUseGhost = true`), `high_to_low_dbl_side_sets` on.

**Physics (FEM).** Per phase: `STRUC_LIN_ISO` and `DIFF_LIN_ISO` (conductivity only) CMs,
`STRUC_LINEAR_BULK` and `SPATIALDIFF_BULK` IWGs. Boundary IWGs on phase-1 side sets: two
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE`, `SPATIALDIFF_NEUMANN`, and
`SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE`, each Nitsche parameter `DIRICHLET_NITSCHE` with penalty
100.0. Interface IWGs `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` and
`STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE` with `NITSCHE_INTERFACE` parameters (penalty 100.0).
Ghost: four `GHOST_NORMAL_FIELD` IWGs (TEMP and UX, UY on `ghost_p1` and `ghost_p0`) with
`GHOST_DISPL` parameters of 0.01. IQIs: `DOF` (TEMP), `STRAIN_ENERGY` with `CMStrucLinIso2` on
phase 0, `VOLUME` on phase 0, and `MAX_DOF` on phase 1. The computation block sets
`POINT_3_CENTRAL` with perturbation 1.0e-4. The set names follow the XTK grammar
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** `AMESOS_IMPL`, no preconditioner, relaxation 1.0, relative residual drop
1.0e-12, `NLA_combined_res_jac_assembly = false`; the final adjoint vector is saved with
`SOL_save_final_adjoint_vec_to_file`. The decks differ in the nonlinear solver (see below).

**Optimization (OPT).** `SWEEP`, one evaluation per ADV, `finite_difference_type = "all"`,
`finite_difference_epsilons = 1e-6`, output `shape_opt_test.hdf5`.

**Output (VIS).** TEMP (nodal) and `MAX_DOF` (global) on the phase-1 sets; `Output_Criterion`
returns `false`.

## Deck variants

| | `..._Thermoelastic` | `..._Thermoelastic_Staggered` |
|---|---|---|
| `NLA_max_iter` | 3 | 4 |
| Nonlinear algorithms | one Newton | Newton (0) and `NLBGS_SOLVER` (1) |
| Nonlinear solvers | one on `UX,UY,TEMP` | Newton on `UX,UY` (0); Newton on `TEMP` with secondary DOFs `UX,UY` (1); NLBGS over `UX,UY;TEMP` with `NLA_Sub_Nonlinear_Solver = "1,0"` (2) |
| Time solver | default, `TSA_DofTypes = "UX,UY,TEMP"` | `TSA_Nonlinear_Solver = 2`, `TSA_DofTypes = "UX,UY;TEMP"` |
| Adjoint file | `Shape_Sensitivity_Two_Material.hdf5` | `Shape_Sensitivity_Two_Material_Staggered.hdf5` |

## What the test checks

The ctest runs two Catch2 test cases, one per deck (both set `gPrintReferenceValues = false`).
Each requires a return code of 0, then:

- compares `LHS/Values` from the deck's adjoint file with `Shape_Sensitivity_Two_Material_Ref.hdf5`
  (the same reference for both decks): equal element counts and `computed − reference < 3e-10`
  per entry (signed difference);
- compares the analytic objective and constraint gradients in `shape_opt_test.hdf5` with the
  `fd_forward`, `fd_backward`, and `fd_central` gradients per ADV: relative difference below
  1e-4 when the FD value exceeds 1e-6 in magnitude, otherwise an absolute margin of 1e-6.

## Run it

```bash
ctest -R Shape_Sensitivity_Two_Material_Sweep_Thermoelastic-1-procs --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- Both interface IWGs carry the name `IWGInterface`, and the phase-0 structural ghost IWG is
  named `SPGPStruct2`, the same as its stabilization parameter.
- Density, heat-capacity, heat-load, and time-continuity properties are defined but not
  referenced by any CM or IWG.
- The staggered deck's VIS file name is `ShapeSensitivitiesTransientCircle.exo`.
