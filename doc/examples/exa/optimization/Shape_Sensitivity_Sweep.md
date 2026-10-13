# Shape sensitivity sweep with two line cuts

A finite-difference verification of shape sensitivities for a 2D linear-elastic problem cut by
two straight level-set lines. The OPT `SWEEP` algorithm evaluates one design point, computes the
analytic (adjoint) gradients of a strain-energy objective and a volume constraint with respect
to the line parameters, and compares them with forward, backward, and central finite
differences. {doc}`Shape_Sensitivity_Sweep_NewIO <Shape_Sensitivity_Sweep_NewIO>` is the
{doc}`Deck API v2 <../../../user-guide/deck-api-v2>` twin of this deck and passes the same checks.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Sweep
```

## Problem

| | |
|---|---|
| Domain | 2.0 × 2.0 square, offset to [−1, 1] × [−1, 1] |
| Geometry | two `LINE` level sets: a vertical line through (0.8, 0.3) with normal (1, 0), and an oblique line through (−0.6, −0.3) with normal (1, 1) |
| Material | bulk phase 1, linear elastic, isotropic, E = 1.0, ν = 0.0 |
| Support | zero displacement on the interface between phase 1 and phase 3 (`iside_b0_1_b1_3`), symmetric Nitsche |
| Load | traction (0, 0.1) on the interface between phase 1 and phase 0 (`iside_b0_1_b1_0`) |
| Objective | strain energy of phase 1 (`IQIBulkStrainEnergy`) |
| Constraint | volume of phase 1 (`IQIBulkVolume`, density property 1.0) |
| Analysis | static, one time step |

## Setup in MORIS

**Mesh (HMR).** 2 × 2 elements on the 2.0 × 2.0 domain, linear Lagrange and linear B-spline
meshes on pattern 0, `refinement_buffer` and `staircase_buffer` 3, `adaptive_refinement_level` 1.

**Geometry and phases (GEN).** The two lines use the default phase table (four phases for two
geometries). Parameters given as three values are ADVs in the order (lower bound, initial value,
upper bound); all bounds here equal the initial value. The ADVs are: `center_x` and `normal_y` of
the vertical line, and `center_x`, `center_y`, `normal_x`, `normal_y` of the oblique line (six in
total). `IQI_types` is `IQIBulkStrainEnergy, IQIBulkVolume`, so `aCriteria(0)` is the strain energy and
`aCriteria(1)` the volume
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, `high_to_low_dbl_side_sets` on, Exodus output of the integration mesh.

**Physics (FEM).** One `STRUC_LIN_ISO` constitutive model on the phase-1 bulk sets
(`HMR_dummy_n_p1,HMR_dummy_c_p1`). IWGs: `STRUC_LINEAR_BULK`,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` with a `DIRICHLET_NITSCHE` parameter (penalty 100.0),
and `STRUC_LINEAR_NEUMANN`. IQIs: `DOF` for UX and UY (output), `STRAIN_ENERGY`, and `VOLUME`.
The computation block sets `finite_difference_scheme = POINT_3_CENTRAL` with perturbation
`1.0e-7`. Set names follow the XTK grammar
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** `BELOS_IMPL` with an `IFPACK` ILU preconditioner, default Newton on UX, UY,
one time step over a time frame of 1.0, and `Sensitivity_Analysis_Type = ADJOINT`.

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1`, `include_bounds = false`,
objective and constraint gradients evaluated, `finite_difference_type = "all"` (forward,
backward, central), and `finite_difference_epsilons = 1.0e-7`. Results go to
`shape_opt_test.hdf5`. The callbacks return the objective as `aCriteria(0)` and the constraint as
`aCriteria(1)`; `get_constraint_types` returns 1; the explicit ADV derivatives are zero. See
[workflow step 4](../../../user-guide/workflow.md#4-sensitivity-sweep) for the sweep procedure.

**Output (VIS).** `shape_sensitivities.exo` on the phase-1 bulk sets with nodal UX and UY,
written every step.

## What the test checks

The test runs `Shape_Sensitivity_Sweep.so` through the workflow interface and requires a return
code of 0. On one processor it then reads `shape_opt_test.hdf5`:

- the analytic gradients `objective_gradients eval_1-1 analytical` and
  `constraint_gradients eval_1-1 analytical`, and requires both to have the same length;
- for each of `fd_forward`, `fd_backward`, `fd_central`, the gradients at
  `eval_1-1 epsilon_1-1`, and requires matching lengths;
- for every ADV, `CHECK(analytical == Approx(fd).margin(1e-16))` for the objective and the
  constraint. With the 1e-16 margin, the comparison is governed by Catch2's default relative
  tolerance for `Approx` (100 × single-precision machine epsilon, ≈ 1.2e-5).

The percent error per ADV is written to the log for each FD scheme.

## Run it

```bash
ctest -R '^Shape_Sensitivity_Sweep-1-procs$' --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```

The anchored pattern keeps `-R` from also selecting `Shape_Sensitivity_Sweep_NewIO-1-procs`.
