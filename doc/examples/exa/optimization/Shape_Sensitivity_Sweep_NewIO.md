# Shape sensitivity sweep with two line cuts (Deck API v2)

The {doc}`Deck API v2 <../../../user-guide/deck-api-v2>` port of
{doc}`Shape_Sensitivity_Sweep <Shape_Sensitivity_Sweep>`: the same 2D linear-elastic problem
cut by two level-set lines, verified by the same `SWEEP` comparison of analytic and
finite-difference shape sensitivities. The deck has a single `MORIS_DECK` entry point, no
`extern "C"` module functions, no `Func_Const` or `Output_Criterion` callbacks (built-ins are
used), and no OPT callback functions: the objective and constraint are criteria expressions,
the mesh-set strings come from the typed vocabulary, and the FEM block is emitted by the
linear-elastic preset. The source header records 416 lines for the legacy twin; this deck has 187.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Sweep_NewIO
```

## Problem

Identical to {doc}`Shape_Sensitivity_Sweep <Shape_Sensitivity_Sweep>`:

| | |
|---|---|
| Domain | 2.0 × 2.0 square, offset to [−1, 1] × [−1, 1] |
| Geometry | vertical `LINE` through (0.8, 0.3), normal (1, 0); oblique `LINE` through (−0.6, −0.3), normal (1, 1) |
| Material | bulk phase 1, linear elastic, preset defaults E = 1.0, ν = 0.0, density 1.0 |
| Support | zero displacement on the phase 1 / phase 3 interface, symmetric Nitsche (preset penalty 100.0) |
| Load | traction (0, 0.1) on the phase 1 / phase 0 interface |
| Objective | `IQIBulkStrainEnergy` |
| Constraint | `IQIBulkVolume` ≤ 0 |
| Analysis | static, one time step |

## What the v2 version replaces

| Legacy deck | This deck |
|---|---|
| eight `*ParameterList` functions under `extern "C"` | one `MORIS_DECK( aDeck )` body with `aDeck.hmr()`, `aDeck.gen()`, … |
| `compute_objectives`, `compute_constraints`, four derivative callbacks, `get_constraint_types` | `aDeck.objective( deck::criterion( "IQIBulkStrainEnergy" ) )` and `aDeck.constraint( deck::criterion( "IQIBulkVolume" ) <= 0.0 )` |
| hand-written `IQI_types` in GEN | derived from the expressions ([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)) |
| literal set strings | `deck::Phase( 1 ).bulk()`, `deck::interface( tMaterial, deck::Phase( 3 ) )`, `deck::interface( tMaterial, deck::Phase( 0 ) )` |
| hand-written properties, CM, SP, three IWGs, strain-energy and volume IQIs | `fem::presets::linear_elastic()` with `mBulkSets`, `mDirichletSets`, `mNeumannSets`, `mTraction = "0.0;0.1"` |
| `library` key in OPT | omitted; the loaded deck is handed to OPT |

## Setup in MORIS

**Mesh (HMR).** Same as the legacy deck: 2 × 2 elements, linear Lagrange and B-spline meshes,
buffers 3, `adaptive_refinement_level` 1.

**Geometry and phases (GEN).** The two `LINE` geometries with the same six ADVs (bounds equal to
the initial values): `center_x`, `normal_y` of the vertical line and `center_x`, `center_y`,
`normal_x`, `normal_y` of the oblique line.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, `high_to_low_dbl_side_sets` on.

**Physics (FEM).** The preset emits `STRUC_LIN_ISO`, `STRUC_LINEAR_BULK`,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE`, `STRUC_LINEAR_NEUMANN`, and the
`STRAIN_ENERGY` and `VOLUME` IQIs; no ghost or bedding sets are passed. Raw parameter lists add
the `DOF` IQIs for UX and UY and the computation block (`POINT_3_CENTRAL`, perturbation
`1.0e-7`). Set names come from `deck::Phase` handles
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** `BELOS_IMPL` with `IFPACK` ILU, one static step, adjoint sensitivities, and the
built-in always-true `Output_Criterion`.

**Optimization (OPT).** `SWEEP` with the legacy settings: one evaluation per ADV, bounds
excluded, `finite_difference_type = "all"`, epsilon `1.0e-7`, output `shape_opt_test.hdf5`.
Criteria gradients come from reverse-mode differentiation of the expressions.

**Output (VIS).** `shape_sensitivities_newio.exo` with nodal UX and UY on the phase-1 bulk sets.

## What the test checks

The test case is the legacy test with renamed identifiers: it runs
`Shape_Sensitivity_Sweep_NewIO.so`, requires return code 0, and on one processor compares the
analytic objective and constraint gradients in `shape_opt_test.hdf5` against the `fd_forward`,
`fd_backward`, and `fd_central` gradients per ADV with `Approx(...).margin(1e-16)` (Catch2's
default relative tolerance, ≈ 1.2e-5, governs). See
{doc}`Shape_Sensitivity_Sweep <Shape_Sensitivity_Sweep>` for the full list of checks.

## Run it

```bash
ctest -R Shape_Sensitivity_Sweep_NewIO --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```
