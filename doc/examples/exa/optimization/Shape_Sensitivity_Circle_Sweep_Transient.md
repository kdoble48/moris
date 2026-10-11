# Shape sensitivity of a transient heat-conduction circle problem

A finite-difference verification of shape sensitivities for transient linear heat conduction in
a 2D domain cut by a circular level set. The OPT `SWEEP` algorithm computes the analytic
gradients of a diffusion strain-energy objective and a volume constraint with respect to the
circle's center and radius and compares them with forward, backward, and central finite
differences. The test also compares the final adjoint vector against a stored reference. Two
decks solve the same problem with different solver configurations. The thermoelastic
counterpart is
{doc}`Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient <Shape_Sensitivity_Circle_Sweep_Thermoelastic_Transient>`.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Circle_Sweep_Transient
```

## Problem

| | |
|---|---|
| Domain | 6.0 × 4.0 rectangle (`sL = 6.0`) |
| Geometry | `CIRCLE` level set, center (3.0, 2.21), radius 1.4 |
| Phases | bulk phase 1 is the material (`HMR_dummy_n_p1,HMR_dummy_c_p1`) |
| Material | `DIFF_LIN_ISO`: conductivity 1.0, density 1.0, heat capacity 1.0 |
| Body load | volumetric heat load 1.0 |
| Boundary conditions | TEMP = 0 on `SideSet_4_n_p1` (x-min side, unsymmetric Nitsche); surface flux 5.0 on `SideSet_2_n_p1` (x-max side) |
| Initial condition | TEMP = 0 |
| Time | transient, 2 time steps over a time frame of 1.0 |
| Objective | `IQIBulkStrainEnergy` (`STRAIN_ENERGY` evaluated with the diffusion CM) |
| Constraint | `IQIBulkVolume` |

## Setup in MORIS

**Mesh (HMR).** 6 × 4 elements, linear Lagrange and linear B-spline meshes on pattern 0,
`pattern_initial_refinement` 1, `refinement_buffer` and `staircase_buffer` 1.

**Geometry and phases (GEN).** One `CIRCLE` field with the default phase table. `center_x`,
`center_y`, and `radius` are each given as (lower bound, initial value, upper bound) with all
three values equal, so the sweep has three ADVs. `IQI_types` is
`IQIBulkStrainEnergy, IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization off (`tUseGhost = false`), verbose output, Exodus output of the integration mesh.

**Physics (FEM).** IWGs on the phase-1 sets: `SPATIALDIFF_BULK` (with the heat load),
`SPATIALDIFF_NEUMANN` on the x-max side set, `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` with a
`DIRICHLET_NITSCHE` parameter (penalty 10.0) on the x-min side set, and `TIME_CONTINUITY_DOF`
(current and previous weights 100.0, initial condition 0.0, `time_continuity = true`). IQIs:
`DOF` (TEMP), `STRAIN_ENERGY`, `VOLUME`, and `MAX_DOF` (function parameters `1.0/2.0`). Side-set
and bulk-set names follow the XTK grammar
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** Nonlinear relative residual drop 1.0e-8, relaxation 1.0, at most 2 iterations;
two time levels for TEMP; the final adjoint vector is saved to an HDF5 file
(`SOL_save_final_adjoint_vec_to_file`). The two decks differ here (see below).

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1` and
`finite_difference_type = "all"`; `finite_difference_epsilons` is not set, so the default
`1E-8` from `fn_PRM_OPT_Parameters.hpp` applies. Results go to `shape_opt_test.hdf5`. The
callbacks return `aCriteria(0)` as the objective and `aCriteria(1)` as the constraint.

**Output (VIS).** TEMP (nodal) and `MAX_DOF` (global) on the phase-1 sets. `Output_Criterion`
returns `false` in both decks.

## Deck variants

| | `Shape_Sensitivity_Circle_Sweep_Transient` | `Shape_Sensitivity_Circle_Sweep_Transient_Vol_2` |
|---|---|---|
| Linear algorithms | `AMESOS_IMPL` | `AMESOS_IMPL` (0) and `BELOS_IMPL` (1) with preconditioner 0 |
| Preconditioner | `NONE` | `IFPACK`, `ILU` |
| Linear solvers | one | two; the second uses `DLA_Linear_solver_algorithms = "1"` |
| Nonlinear algorithms | one, `NLA_combined_res_jac_assembly = false` | two, both with combined assembly `true`; the second sets `NLA_Linear_solver = 1` |
| Nonlinear solvers | one on TEMP | two on TEMP; the second sets `NLA_Sub_Nonlinear_Solver = "1"` |
| Time solver | default | `TSA_Nonlinear_Sensitivity_Solver = 1` |
| Adjoint file | `Shape_Sensitivity_Transient.hdf5` | `Shape_Sensitivity_Transient_V2.hdf5` |
| VIS file | `ShapeSensitivitiesTransientCircle.exo` | `ShapeSensitivitiesTransientCircle_V2.exo` |

Geometry, physics, and sweep settings are identical.

## What the test checks

The ctest runs two Catch2 test cases in sequence, one per deck. Each requires a return code of
0 and then:

- reads `LHS/Values` from the deck's adjoint file and from
  `Shape_Sensitivity_Transient_Ref.hdf5` (the same reference for both decks), checks equal
  element counts, and checks `computed − reference < 1e-10` for every entry (a signed
  difference, not an absolute value);
- reads the analytic objective and constraint gradients and the `fd_forward`, `fd_backward`,
  and `fd_central` gradients from `shape_opt_test.hdf5`, and for every ADV requires a relative
  difference below 1e-4 when the FD value exceeds 1e-6 in magnitude, otherwise agreement within
  an absolute margin of 1e-6.

## Run it

```bash
ctest -R Shape_Sensitivity_Circle_Sweep_Transient-1-procs --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- `sL1`, `sL2`, `tInterfaceRefinement`, `tPhase1Ghost`, and the ghost stabilization branches are
  defined in both decks but not used, since `tUseGhost = false`.
