# Shape sensitivity of a thermoelastic plate with a circular hole

A coupled thermoelastic plate with a circular immersed hole, whose center and radius are the
ADVs. Two decks run a SWEEP that compares analytical and finite-difference gradients of a
strain-energy objective and a volume constraint, one with a monolithic Newton solve and one
with a staggered (block Gauss–Seidel) solve of the temperature and displacement fields. A third
deck runs two GCMMA iterations on the same problem with non-degenerate ADV bounds.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Circle_Sweep_Thermoelastic
```

## Problem

| | |
|---|---|
| Domain | 6.0 × 4.0 rectangle |
| Geometry | circle, center (3.0, 2.21), radius 1.4 |
| Phases | phase 1 outside the circle (material), phase 0 inside (no physics) |
| Structure | linear elastic, E = 1.0, ν = 0.0, CTE 1.0, reference temperature 0.0 |
| Heat conduction | conductivity 1.0, density 1.0, heat capacity 1.0, body heat load 1.0 |
| Thermal BCs | temperature 0.0 on the left edge (side set 4), heat flux 5.0 on the right edge (side set 2) |
| Structural BCs | zero displacement on the left and right edges, unsymmetric Nitsche |
| Analysis | steady, one time step, time frame 1.0; temperature drives the displacement through thermal expansion |
| Objective | strain energy of the structure |
| Constraint | volume of phase 1 |

## Setup in MORIS

**Mesh (HMR).** 6 × 4 elements (`sDim = 2`), linear Lagrange and B-spline meshes (pattern 0),
one level of initial refinement (`pattern_initial_refinement = 1`), refinement and staircase
buffers 1.

**Geometry and phases (GEN).** One `CIRCLE` level-set geometry. `center_x`, `center_y` and
`radius` are each given as (lower bound, initial value, upper bound), which makes them three
ADVs. The default phase table puts the positive side (outside the circle) in phase 1.
`IQI_types` is `IQIBulkStrainEnergyDISP`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0,
`high_to_low_dbl_side_sets`, Exodus output of the integration mesh. Ghost stabilization is off
(`tUseGhost = false`).

**Physics (FEM).** All IWGs and IQIs act on `HMR_dummy_{n,c}_p1` and its side sets
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
The `STRUC_LIN_ISO` constitutive model depends on `UX,UY` and `TEMP` and takes Young's modulus,
Poisson ratio, CTE and reference temperature; the `DIFF_LIN_ISO` model depends on `TEMP`. IWGs:
`STRUC_LINEAR_BULK` and `SPATIALDIFF_BULK` (heat load 1.0);
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on `SideSet_4_n_p1` and on `SideSet_2_n_p1`
(`DIRICHLET_NITSCHE`, 100.0); `SPATIALDIFF_NEUMANN` (flux 5.0) on `SideSet_2_n_p1`;
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (temperature 0.0, `DIRICHLET_NITSCHE` 10.0) on
`SideSet_4_n_p1`. IQIs: `STRAIN_ENERGY` of the structural model, `VOLUME`, nodal `DOF` for TEMP,
and `MAX_DOF` of TEMP (parameters 1.0/2.0).

**Solver (SOL).** `AMESOS_IMPL` direct solve, no preconditioner; Newton with relative residual
drop 10⁻¹⁰ (10⁻⁸ in the GCMMA deck), relaxation 1.0, at most 2 iterations; one time step. The
staggered variant is described below.

**Optimization (OPT).** The sweep decks use `SWEEP` with `num_evaluations_per_adv = 1`,
`finite_difference_type = "all"`, `finite_difference_epsilons = 1e-6` and output to
`shape_opt_test.hdf5`; `include_bounds` is left at its default. Objective = criterion 0,
constraint = criterion 1.

**Output (VIS).** `ShapeSensitivitiesTransientCircle.exo` on phase 1 with nodal TEMP and global
`MAX_DOF`. The deck's `Output_Criterion` returns `false`.

## The three decks

| | Base | `_Staggered` | `_GCMMA` |
|---|---|---|---|
| Algorithm | `SWEEP` | `SWEEP` | `GCMMA`, `max_its = 2` |
| ADV bounds (center_x, center_y, radius) | (3.0, 3.0, 3.0), (2.21, 2.21, 2.21), (1.4, 1.4, 1.4) | same as base | (2.0, 3.0, 4.0), (2.0, 2.21, 3.0), (1.0, 1.4, 2.0) |
| Nonlinear solve | monolithic Newton on `UX,UY,TEMP` | `NLBGS_SOLVER` over two Newton sub-solvers | monolithic Newton |
| Structural `model_type` | `PLANE_STRESS` | `PLANE_STRESS` | `FULL` |
| Newton residual drop | 10⁻¹⁰ | 10⁻¹⁰ | 10⁻⁸ |

**Staggered solve.** The staggered deck defines two nonlinear algorithms (Newton with the
settings above, and `NLBGS_SOLVER`) and three nonlinear solvers: solver 0, Newton on `UX,UY`;
solver 1, Newton on `TEMP` with `UX,UY` as secondary DOFs; solver 2, `NLBGS_SOLVER` on
`UX,UY;TEMP` with sub-solvers `"1,0"` (temperature, then displacement) and algorithm 1. The time
solver uses solver 2 (`TSA_Nonlinear_Solver = 2`) with `TSA_DofTypes = "UX,UY;TEMP"`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

In the GCMMA deck the strain-energy IQI depends only on `UX,UY` and the volume IQI only on
`TEMP`; in the other two decks both depend on `UX,UY;TEMP`.

## What the test checks

- **`Shape_Sensitivity_Circle_Sweep_Thermoelastic`** (runs only on 1 processor): return code 0;
  for every ADV and each FD type (forward, backward, central), if the |FD gradient| exceeds 10⁻⁶
  the relative difference must be below 10⁻⁴, otherwise the analytical value must match within
  an absolute margin of 10⁻⁶ (`CHECK`, for objective and constraint).
- **`..._Staggered`** (runs only on 1 processor): return code 0; analytical and FD gradients
  must agree under Catch2 `Approx(...).scale(1000)` for objective and constraint.
- **`..._GCMMA`**: return code 0 only, and only when MORIS is built with `MORIS_HAVE_GCMMA`.

## Run it

```bash
ctest -R Shape_Sensitivity_Circle_Sweep_Thermoelastic-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Shape_Sensitivity_Circle_Sweep_Thermoelastic-2-procs --output-on-failure
```

On 2 processors the two sweep test cases skip their body, so only the GCMMA case executes.

## Notes

- The deck defines time-continuity properties (`PropWeightCurrent`, `PropWeightPrevious`,
  `PropInitialCondition`) that no IWG references; there is no time-continuity IWG.
- The length constants `sL1` (4.1) and `sL2` are defined but unused.
