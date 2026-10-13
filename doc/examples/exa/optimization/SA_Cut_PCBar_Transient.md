# Cut bar with phase change, transient heat conduction: sensitivity sweep

A finite-difference check of shape sensitivities for a transient heat-conduction problem with a
phase-change material. The geometry is the bar of
{doc}`SA_Cut_Bar_Static <SA_Cut_Bar_Static>`: a horizontal immersed interface, whose height is
the single ADV, separates a conductive fin on top from a phase-change material below. The bar
is heated through its top for 20 time steps; the SWEEP algorithm compares the analytical
(adjoint) gradients of a maximum-temperature objective and a volume constraint with forward,
backward and central finite differences.

```{exa-summary} projects/EXA/optimization/SA_Cut_PCBar_Transient
```

## Problem

| | |
|---|---|
| Background domain | 0.16 × 0.6, offset (−0.0342356, −0.031345) |
| Interface | horizontal line y = s, with s the ADV: initial value and lower bound 0.425, upper bound 0.475 |
| Phases | phase 3 above the interface (fin), phase 2 below (phase-change material) |
| Fin material | conductivity 50.0, density 1.0, heat capacity 5.0 |
| Phase-change material | conductivity 5.0, density 1.0, heat capacity 50.0, latent heat 500.0, phase-change temperature 20.0, melting range 6.0, phase-state function parameter 2.0 |
| Heat flux | 1000.0 on side set 3 (top) of the fin |
| Initial temperature | 10.0, imposed weakly through the time-continuity IWG |
| Analysis | transient, 20 time steps over a time frame of 0.1 |
| Objective | `MAX_DOF` of TEMP, parameters reference 1.0 / exponent 2.0 |
| Constraint | volume of the phase-change material |

Unlike the static deck, there is no prescribed-temperature boundary condition.

## Setup in MORIS

**Mesh (HMR).** Same as the static deck: 1 × 40 linear elements, refinement and staircase
buffers 0.

**Geometry and phases (GEN).** Same as the static deck: a `LINE` outside the domain (center
x = −1.0, normal (1, 0)) and a horizontal `LINE` whose `center_y` is the ADV, bounded by
(0.425, 0.475). `IQI_types` is `IQIMaxTemp`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization, `high_to_low_dbl_side_sets`, Exodus output of the enriched and integration meshes.

**Physics (FEM).** The fin uses `DIFF_LIN_ISO`; the storage material uses `DIFF_LIN_ISO_PC`
(phase change) with the latent heat, phase-change temperature, phase-state function and
melting-range properties (`tIsPhaseChange = true`; the `false` branch falls back to
`DIFF_LIN_ISO`). Both bulk IWGs are `SPATIALDIFF_BULK` with `GGLS_DIFFUSION` stabilization; the
phase-change GGLS parameter receives the real latent-heat properties. The interface uses
`SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, 100.0) on
`dbl_iside_p0_3_p1_2`; the flux uses `SPATIALDIFF_NEUMANN` on `SideSet_3_{n,c}_p3`. A
`TIME_CONTINUITY_DOF` IWG on both materials (`time_continuity = true`) uses weights 100.0 for the
current and previous step and the initial-condition property `Func_Initial_Temperature`. Ghost
penalty is `GHOST_NORMAL_FIELD` with `GHOST_DISPL` (0.01) on `ghost_p3` and `ghost_p1`. IQIs:
`MAX_DOF`, nodal `DOF`, and `VOLUME` on the storage material. The `COMPUTATION` list sets
`POINT_3_CENTRAL` differencing with perturbation 10⁻⁸ for the element-level sensitivities.

**Solver (SOL).** `BELOS_IMPL` with an `IFPACK` ILU preconditioner. Two nonlinear algorithms:
the forward Newton solve (relative residual drop 5 × 10⁻⁷, relaxation 0.5 with the
`InvResNormAdaptive` strategy and damping 0.5, at most 25 iterations) and the adjoint solve (at
most 1 iteration, relaxation 1.0). The time solver selects them through
`TSA_Nonlinear_Solver = 0` and `TSA_Nonlinear_Sensitivity_Solver = 1`, and keeps two time levels
for TEMP (`TSA_time_level_per_type = "TEMP,2"`).

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1`, `include_bounds = false`,
`finite_difference_type = "all"` and `finite_difference_epsilons = 1E-6`; results go to
`SA_Cut_PCBar_Transient.hdf5`.

**Output (VIS).** `SA_Cut_PCBar_Transient.exo` on both materials: nodal TEMP and global
`MAX_DOF`.

## Differences from the static deck

| | {doc}`SA_Cut_Bar_Static <SA_Cut_Bar_Static>` | SA_Cut_PCBar_Transient |
|---|---|---|
| Time integration | default (steady) | 20 steps, frame 0.1, time-continuity IWG |
| Storage-material model | `DIFF_LIN_ISO`, k = 1.0 | `DIFF_LIN_ISO_PC`, k = 5.0, latent heat 500.0 |
| Fin heat capacity | 50.0 | 5.0 |
| Bottom temperature BC | 10.0 (Nitsche) | none |
| Newton | 10 iterations, relaxation 1.0 | 25 iterations, adaptive relaxation 0.5; separate adjoint solver |
| Sweep FD epsilon | default 10⁻⁸ | 10⁻⁶ |
| `MAX_DOF` reference | 10.0 | 1.0 |

## What the test checks

The test requires a return code of 0, then reads the analytical and FD gradients (forward,
backward, central) from `SA_Cut_PCBar_Transient.hdf5`. For every ADV it checks (`CHECK`) that the
relative difference is below 0.05 for the objective and below 0.001 for the constraint. It runs
on 1 processor.

## Run it

```bash
ctest -R SA_Cut_PCBar_Transient-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The PCM ghost IWG uses the set `ghost_p1`, while the PCM bulk sets are on phase 2 (the static
  deck uses `ghost_p2`).
- `Func_Neumann_BC`, which switches the flux on only for y ≥ 0.45 (0.9 × `tYlength`), is defined but not
  used; `PropNeumannBC` uses `Func_Const`.
- `tHaveDirichlet`, `tHaveGhost` and `tIsPhaseChange` at the top of the deck toggle the interface
  Nitsche IWG, the ghost IWGs and the phase-change model.
