# Two-bar truss sensitivity sweep

Two tapered bars, described by one user-defined level-set function, form a truss that is clamped
on the left edge and loaded on the right. The deck runs the OPT `SWEEP` algorithm, which evaluates
the objective and constraint and their gradients analytically and by forward, backward and
central finite differences, and writes them to HDF5. It demonstrates a parametric user geometry
with a user-supplied sensitivity function, OPT callbacks indexing GEN `IQI_types`, and an
iterative Belos linear solve with an ILU preconditioner.

```{exa-summary} projects/EXA/structure/linear/Two_Bar_Truss
```

## Problem

| | |
|---|---|
| Domain | unit square 1.0 × 1.0 |
| Bars | bar 1 from (0, 0) to (1, 0.5), bar 2 from (0, 1) to (1, 0.5), end radii 0.1 |
| Phases | `p0` (inside the bars) is the only analyzed phase |
| Material | E = 1.0, ν = 0.3, density 1.0; bedding 10⁻⁶ |
| Support | left edge (`SideSet_4_n_p0`, `SideSet_4_c_p0`), displacement (0, 0) by Nitsche |
| Load | traction (0, 1.0) on the right edge (`SideSet_2_*_p0`) for 0.475 < y < 0.52, zero elsewhere |
| Objective | strain energy / 2.7 |
| Constraint | volume / 1.9 − 1.0 (constraint type 1) |

## Setup in MORIS

**Mesh (HMR).** 20 × 20 elements, linear Lagrange and B-spline meshes (pattern 0), refinement and
staircase buffer 1.

**Geometry and phases (GEN).** One `USER_DEFINED` geometry with `field_function_name = "Bars"`
and `sensitivity_function_name = "BarsGrad"`. `Bars` takes 6 parameters per bar (end-point A,
radius at A, end-point B, radius at B) and combines the bars' signed distances with a
Kreisselmeier–Steinhauser aggregate (β = −100); beyond an elongation factor of 1.1 of the bar
length the distance to the end point is used. The 12 parameters are inserted as
`variable_01` … `variable_12`, each a `Design_Variable` whose three arguments are equal.
`IQI_types` is `IQIBulkStrainEnergy, IQIBulkVolume`, so `aCriteria(0)` is the strain energy and
`aCriteria(1)` the volume
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, Exodus output of the integration mesh.

**Physics (FEM).** `STRUC_LIN_ISO` on `UX,UY`. IWGs: `STRUC_LINEAR_BULK` on the bars,
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on the left edge (`DIRICHLET_NITSCHE` 100.0),
`STRUC_LINEAR_NEUMANN` on the right edge with traction from `Func_Neumann`, and
`GHOST_NORMAL_FIELD` on `ghost_p0` (`GHOST_DISPL` 0.01). IQIs: `DOF` UX and UY, `STRAIN_ENERGY`,
`MAX_DOF` on UY (`function_parameters = "1.0/2.0"`), and `VOLUME`. The computation block sets a
3-point central finite-difference scheme with perturbation 10⁻⁴.

**Solver (SOL).** `BELOS_IMPL` with preconditioner 0, an `IFPACK` preconditioner of type `ILU`,
default Newton, DOFs `UX,UY`.

**Optimization (OPT).** `problem = "user_defined"`, `library = "./Two_Bar_Truss.so"`, algorithm
`SWEEP` with `hdf5_path = "TwoBarTruss.hdf5"`, objective and constraint gradients evaluated,
`num_evaluations_per_adv = "1"`, `include_bounds = false`, `finite_difference_type = "all"`,
`finite_difference_epsilons = "1e-6"`
([coupling C2](../../../user-guide/interlinks.md#c2-problem-declaration-opt)). The callbacks
return d(objective)/d(criteria) = (1/2.7, 0) and d(constraint)/d(criteria) = (0, 1/1.9) and zero
explicit ADV derivatives.

**Output (VIS).** `TwoBarTruss.exo` on the bars: nodal UX, UY and global `STRAIN_ENERGY`, `MAXUY`,
`VOLUME`.

## What the test checks

After the workflow call (its return code is not checked), the test opens `TwoBarTruss.hdf5` and
reads the analytical objective and constraint gradients (`objective_gradients eval_1-1
analytical`, `constraint_gradients eval_1-1 analytical`) and, for each of `fd_forward`,
`fd_backward` and `fd_central` at epsilon 1-1, the finite-difference gradients. It requires only
that all these vectors have equal length. Per-ADV percent differences are written to the log;
the value comparisons (`CHECK( ... == Approx( ... ) )`) are commented out, so gradient agreement
is not asserted.

## Run it

```bash
ctest -R Two_Bar_Truss-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

To run a sweep of this kind on your own deck, see {doc}`../../../user-guide/workflow`.

## Notes

- In `BarsGrad`, the perturbed value is applied to the local `tNominalValue` but not written back
  into `tGeometryParameters` before `Bars` is called, so the positive and negative evaluations
  use the same parameters.
- The logged constraint comparison prints `tConstraintsAnalytical` in place of the
  finite-difference value.
- The `.dox` page lists "Belos Solver" as the only feature.
