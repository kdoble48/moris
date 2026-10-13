# Solver configurations for a transient thermoelastic phase-change problem

One deck, run under several solver configurations: a 2D two-material, transient thermoelastic
problem in which the second material has a phase-change (latent heat) diffusion model. The deck
performs one GCMMA iteration, so each configuration runs the forward and the sensitivity
analysis. Global flags set by the test select a monolithic or staggered (NLBGS) forward solve, a
monolithic or staggered adjoint solve, and single or double time levels for the displacement
DOFs. The test checks only that each configuration completes.

```{exa-summary} projects/EXA/optimization/Solver_Examples
```

## Problem

| | |
|---|---|
| Domain | 0.8 × 0.2 rectangle |
| Geometry | `LINE` level set with normal (1, 0) at `center_x` = 0.4, splitting the domain into phase 0 (`Phase1`, material 1) and phase 1 (`Phase2`, material 2) |
| Material 1 | `DIFF_LIN_ISO`: density 1.0, heat capacity 5.0, conductivity 50.0; `STRUC_LIN_ISO` (`PLANE_STRESS`): E = 100.0, ν = 0.31, CTE 1.0e-2 |
| Material 2 (PCM) | `DIFF_LIN_ISO_PC`: density 1.0, heat capacity 50.0, conductivity 5.0, latent heat 500.0, phase-change temperature 20.0, phase-change constant 6.0, phase-state function parameter 2.0; `STRUC_LIN_ISO` (`PLANE_STRESS`): E = 10.0, ν = 0.31, CTE 1.0e-3 |
| Reference temperature | 1.0 |
| Thermal load | heat flux 1000.0 on side ordinal 4 (x-min) of material 1 |
| Structural load | pressure 1.0 on side ordinal 3 (y-max) of both materials |
| Support | zero displacement on side ordinal 4 of material 1, unsymmetric Nitsche |
| Initial condition | `TSA_Initialize_Sol_Vec` sets UX = UY = 0, TEMP = 1.0; the time-continuity initial-condition property returns 10.0 (`Func_Initial_Condition`) |
| Time | 2 time steps over a time frame of 0.1 |
| Objective | `IQIMaxTemp_1 + IQIMaxTemp_2 + IQIStrainEnergy_1 + IQIStrainEnergy_2` |
| Constraint | `IQIVolume_1` (`get_constraint_types` returns 1) |

## Setup in MORIS

**Mesh (HMR).** 32 × 8 elements, linear Lagrange and linear B-spline meshes on pattern 0.

**Geometry and phases (GEN).** One `LINE` geometry; `center_x` is the only ADV, with lower bound
0.1, initial value 0.4, and upper bound 0.7. `isocontour_tolerance` is 1e-12 and GEN writes
`GEN_Solver_Examples.exo`. `IQI_types` lists five criteria in the order `IQIMaxTemp_1`,
`IQIMaxTemp_2`, `IQIStrainEnergy_1`, `IQIStrainEnergy_2`, `IQIVolume_1`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, `keep_all_opt_iters` on, enriched-mesh printing and Exodus output of the
integration mesh.

**Physics (FEM).** The deck declares FEM phases (`FEM::PHASES`: `Phase1` = index 0, `Phase2` =
index 1) and attaches every CM, SP, IWG, and IQI by `leader_phase_name` / `follower_phase_name`,
side sets by `side_ordinals`, instead of writing mesh-set strings
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
IWGs:

- bulk: `SPATIALDIFF_BULK` per phase with `GGLS_DIFFUSION` stabilization, `STRUC_LINEAR_BULK` per phase;
- boundary: `SPATIALDIFF_NEUMANN` (heat flux), `STRUC_LINEAR_NEUMANN` on both phases (pressure),
  `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE` (penalty 100.0);
- interface: `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` and `STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE`
  with `NITSCHE_INTERFACE` (penalty 100.0);
- ghost: `GHOST_NORMAL_FIELD` for TEMP and UX, UY in each phase, `GHOST_DISPL` 0.01;
- `TIME_CONTINUITY_DOF` on TEMP (weights 100.0), declared twice on `Phase2`.

IQIs per phase: `DOF` for TEMP, UX, UY; `MAX_DOF` (function parameters `1.0/2.0`);
`STRAIN_ENERGY`; `VOLUME`.

**Solver (SOL).** The configuration is assembled in `SOLParameterList` from four `extern` globals that
the test sets (declared in `EXA_Globals.hpp`):

| Global | Effect |
|---|---|
| `gUseBelosWithILUT` | `true`: `BELOS_IMPL` with an `IFPACK` `ILUT` preconditioner (drop tolerance 1e-10, level of fill 25.0); `false`: `AMESOS_IMPL`, no preconditioner |
| `gHaveStaggeredFA` | `false`: one Newton algorithm (relative residual drop 5.0e-6, relaxation 1.0, at most 15 iterations) on `UX,UY;TEMP`. `true`: Newton on `UX,UY` using a one-iteration Newton algorithm, Newton on `TEMP` (no algorithm index set), and an NLBGS solver (one iteration) with `NLA_Sub_Nonlinear_Solver = "1,0"`; the time solver uses `TSA_Nonlinear_Solver = 2` |
| `gHaveStaggeredSA` | with staggered FA, `false` adds `TSA_Nonlinear_Sensitivity_Solver = 3`, a separate monolithic Newton solver on all DOFs for the adjoint; `true` leaves the adjoint on the staggered solver |
| `gUseMixedTimeElements` | `true` sets one time level for UX and UY (`sTLSD = "1"`); TEMP keeps two |

All configurations write the operator with `SOL_save_operator_to_matlab = "Mat.dat"`.

**Optimization (OPT).** `GCMMA` with `max_its = 1` and `step_size = 0.2`. The objective
derivative with respect to the criteria is 1 for the first four criteria and 0 for the volume;
the constraint derivative is 1 for the volume only.

**Output (VIS).** Two output meshes, one per material (`Solver_Examples.exo_material_1` and
`_material_2`, `Time_Offset` 100.0), each with nodal TEMP, UX, UY and global `MAX_DOF`,
`STRAIN_ENERGY`, `VOLUME`.

## Configurations exercised

| Test case | `gHaveStaggeredFA` | `gHaveStaggeredSA` | `gUseMixedTimeElements` | `gUseBelosWithILUT` |
|---|---|---|---|---|
| `Standard_Monolithic` | false | false | false | false |
| `Staggered_FA_Monolithic_SA` | true | false | false | false |
| `Staggered_FA_and_SA` | true | true | false | false |
| `Staggered_FA_and_SA_Mixed_Time_Elements` | true | true | true | false |

A fifth case, `Monolithic_Mixed_Time_Elements`, is commented out with the note that this
configuration returns an error in matrix assembly due to the mixed time elements. No test case
enables `gUseBelosWithILUT`.

## What the test checks

Each test case sets the four globals, runs `Solver_Examples.so` through the workflow interface,
and requires a return code of 0. No solution values, criteria, or gradients are compared.

## Run it

```bash
ctest -R '^Solver_Examples$' --output-on-failure    # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name has no `-procs` suffix (`NO_PROCS_SUFFIX`) and runs on one processor; all four
test cases run in one ctest.

## Notes

- The legacy `add_test` was also gated on `MORIS_HAVE_GCMMA`; the current registration depends
  only on `MORIS_HAVE_PARALLEL_TESTS` (comment in `CMakeLists.txt`).
- The deck's header comments show a standalone run:
  `create_shared_object.sh . build_dbg Solver_Examples`, then `moris ./Solver_Examples.so`.
