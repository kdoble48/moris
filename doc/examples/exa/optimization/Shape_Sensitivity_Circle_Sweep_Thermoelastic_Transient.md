# Shape sensitivity of a transient thermoelastic circle problem

A finite-difference verification of shape sensitivities for a transient, coupled
thermoelastic problem in a 2D domain cut by a circular level set. The OPT `SWEEP` algorithm
compares the analytic gradients of a structural strain-energy objective and a volume constraint
with respect to the circle's center and radius against forward, backward, and central finite
differences, and the test compares the final adjoint vector with a stored reference. One deck
solves the coupled system monolithically, the other with a staggered (NLBGS) nonlinear solver.
The thermal-only version is
{doc}`Shape_Sensitivity_Circle_Sweep_Transient <Shape_Sensitivity_Circle_Sweep_Transient>`.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient
```

## Problem

| | |
|---|---|
| Domain | 6.0 × 4.0 rectangle |
| Geometry | `CIRCLE` level set, center (3.0, 2.21), radius 1.4 |
| Phases | bulk phase 1 is the material |
| Thermal material | `DIFF_LIN_ISO`: conductivity 1.0, density 1.0, heat capacity 1.0 |
| Structural material | `STRUC_LIN_ISO`, `model_type FULL`: E = 1.0, ν = 0.0, CTE 1.0, reference temperature 0.0 |
| Thermal BCs | TEMP = 0 on `SideSet_4_n_p1` (x-min side); surface flux 5.0 on `SideSet_2_n_p1` (x-max side); body heat load 1.0 |
| Structural BCs | zero displacement on both `SideSet_4_n_p1` and `SideSet_2_n_p1` |
| Initial condition | UX = UY = TEMP = 0 |
| Time | transient, 2 time steps over a time frame of 1.0; two time levels for UX, UY, TEMP |
| Objective | `IQIBulkStrainEnergyDISP` (`STRAIN_ENERGY` of the structural CM) |
| Constraint | `IQIBulkVolume` |

## Setup in MORIS

**Mesh (HMR).** 6 × 4 elements, linear Lagrange and linear B-spline meshes on pattern 0,
`pattern_initial_refinement` 1, buffers 1.

**Geometry and phases (GEN).** One `CIRCLE` field with three ADVs (`center_x`, `center_y`,
`radius`), each with bounds equal to its initial value. `IQI_types` is
`IQIBulkStrainEnergyDISP, IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization off, Exodus output of the integration mesh.

**Physics (FEM).** The structural CM depends on `UX,UY;TEMP` (Displacement, Temperature), so the
temperature drives thermal strain. IWGs on phase 1: `STRUC_LINEAR_BULK`, `SPATIALDIFF_BULK` (with
the heat load), two `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` (x-min and x-max side sets,
`DIRICHLET_NITSCHE` penalty 100.0), `SPATIALDIFF_NEUMANN` (x-max),
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (x-min, penalty 10.0), and `TIME_CONTINUITY_DOF` on TEMP
(weights 100.0, initial condition 0.0). IQIs: `DOF` (TEMP), `STRAIN_ENERGY`, `VOLUME`, and
`MAX_DOF` (function parameters `1.0/2.0`). DOF types agree across FEM and SOL
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Solver (SOL).** `AMESOS_IMPL`, no preconditioner, Newton with relative residual drop 1.0e-8,
relaxation 1.0, at most 2 iterations, `NLA_combined_res_jac_assembly = false`. The final adjoint
vector is saved with `SOL_save_final_adjoint_vec_to_file`. The decks differ in the nonlinear
solver (see below).

**Optimization (OPT).** `SWEEP`, one evaluation per ADV, `finite_difference_type = "all"`;
`finite_difference_epsilons` is not set (default `1E-8`). Results go to `shape_opt_test.hdf5`.

**Output (VIS).** TEMP (nodal) and `MAX_DOF` (global) on the phase-1 sets; `Output_Criterion`
returns `false`.

## Deck variants

| | `..._Thermoelastic_Transient` | `..._Thermoelastic_Transient_Staggered` |
|---|---|---|
| Nonlinear algorithms | one Newton | Newton (0) and `NLBGS_SOLVER` (1) |
| Nonlinear solvers | one Newton on `UX,UY,TEMP` | Newton on `UX,UY` (0); Newton on `TEMP` with secondary DOFs `UX,UY` (1); NLBGS over `UX,UY;TEMP` with `NLA_Sub_Nonlinear_Solver = "1,0"` (2) |
| Time solver | default nonlinear solver, `TSA_DofTypes = "UX,UY,TEMP"` | `TSA_Nonlinear_Solver = 2`, `TSA_DofTypes = "UX,UY;TEMP"` |
| Strain-energy IQI DOF dependencies | `UX,UY` | `UX,UY;TEMP` |
| FEM `finite_difference_scheme` | `POINT_3_CENTRAL` | not set |
| Adjoint file | `Shape_Sensitivity_Thermoelastic_Transient.hdf5` | `Shape_Sensitivity_Thermoelastic_Transient_Staggered.hdf5` |
| VIS file | `ShapeSensitivitiesThermoelasticTransient.exo` | `ShapeSensitivitiesThermoelasticTransientStaggered.exo` |

## What the test checks

The ctest runs two Catch2 test cases, one per deck. Each requires a return code of 0, then:

- compares `LHS/Values` from the deck's adjoint file with
  `Shape_Sensitivity_Thermoelastic_Transient_Ref.hdf5` (the same reference for both decks):
  equal element counts and `computed − reference < 1e-10` per entry (signed difference);
- compares the analytic objective and constraint gradients in `shape_opt_test.hdf5` with the
  `fd_forward`, `fd_backward`, and `fd_central` gradients per ADV: relative difference below
  1e-4 when the FD value exceeds 1e-6 in magnitude, otherwise an absolute margin of 1e-6.

## Run it

```bash
ctest -R Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient-1-procs --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The constraint log line in `example_test_case.cpp` prints the analytic value in the
  finite-difference slot and the objective's relative error; the `CHECK` statements use the
  correct constraint values.
