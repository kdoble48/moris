# U-bend turbulent flow: pseudo-time continuation sensitivity test

A finite-difference sensitivity check for Spalart–Allmaras turbulent flow through a 2D U-bend
channel with a field of solid inclusions whose level-set B-spline coefficients are the design
variables. The forward flow problem is solved by a staggered NLBGS iteration with pseudo-time
continuation and a viscosity ramp; the wall distance for the turbulence model comes from a
heat-method solve on the immersed geometry. The `SWEEP` compares adjoint and finite-difference
gradients of three criteria with respect to one ADV. The directory has no
`example_documentation.dox`.

```{exa-summary} projects/EXA/optimization/UBend_Pseudo_Time_Continuation_Sensitivity_Test
```

## Problem

Lengths are scaled by the reference length 0.1 m (so the channel unit sH = 0.1 m becomes 1.0),
velocities by the inlet velocity, and densities by the fluid density.

| | |
|---|---|
| Design box | x ∈ [2, 12] sH, y ∈ [0, 10] sH, bounded by four `LINE` level sets |
| Walls | two `Func_Thin_Wall` geometries centred at (6.5, 5.0) sH: inner (radius 0.5 sH, between y = 4.5 and 5.5 sH) and outer (radius 2.5 sH, between y = 2.5 and 7.5 sH) |
| Inclusions | 3 × 3 superellipses (exponent 12, radius sH) at x ∈ {3, 7, 11} sH, y ∈ {1, 5, 9} sH |
| Phases | 0 void, 1 fluid, 2 inclusion (`get_phase_index_sharp`) |
| Fluid | ρ = 1.0 kg/m³, ν = 4.0·10⁻⁵ m²/s |
| Inlet | parabolic velocity u = U (1 − ((y − 6.5 sH)/sH)²) for y ∈ [5.5, 7.5] sH on the left background side (side ordinal 4), U = 0.01 m/s; turbulent viscosity 3ν |
| Walls and inclusions | zero velocity and zero turbulent viscosity (Nitsche) |
| Analysis | steady, solved by pseudo-time continuation |

The deck computes the Reynolds number as Re = L U / ν = 0.1 · 0.01 / 4.0·10⁻⁵ = 25. A comment
records the original values (Re = 5000 with U = 2.0).

## Setup in MORIS

**Mesh (HMR).** Background domain 15.0214 × 12.0214 (scaled units) offset by
(−2.0214, −1.0214), with `ceil(dim / 5.0) + 1` = 5 × 4 coarse elements. B-spline mesh 0 is
linear on pattern 0, refined three levels (element edge 0.625); B-spline mesh 1 is quadratic
on pattern 1 with no refinement. The Lagrange mesh is linear on pattern 0.
`use_advanced_T_matrix_scheme = 1`. MSI places all DOFs (THETA, PHID, VX, VY, P, VISCOSITY, L2)
on mesh 0.

**Geometry and phases (GEN).** Seven geometries: four `LINE` fields (y = 0, y = 10 sH,
x = 2 sH, x = 12 sH), two `Func_Thin_Wall` fields and the `Func_Inclusion` field named
`Level_Set_Field`. The last is discretized on B-spline mesh 1 (`discretization_mesh_index = 1`)
with bounds ±10.0 (twice the coarse element size), and its coefficients are the ADVs
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
A `SCALED_FIELD` property (scaling 1.0) maps it to the PDV `LS1` on the fluid and inclusion
sets; FEM reads `LS1` in the heat-method penalty IQIs. `IQI_types` lists
`IQIVolumePowDisp`, `IQIPerimeterItf`, `IQIHeatMethodPenaltyFluid`,
`IQIHeatMethodPenaltyInclusion`, `IQIFluidVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0 only, ghost
stabilization, `high_to_low_dbl_side_sets`, and Exodus output of the integration mesh.

**Physics (FEM).** Flow in the fluid phase uses `FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS`
with `INCOMPRESSIBLE_NS_VELOCITY_BULK`/`INCOMPRESSIBLE_NS_PRESSURE_BULK` (pressure spring 10⁻¹²),
`INCOMPRESSIBLE_FLOW` (parameters `36.0/0.0`), and `SPALART_ALLMARAS_TURBULENCE_BULK` with
the `SPALART_ALLMARAS_TURBULENCE` CM (ft2 = 0.0, alpha = 10.0) and
`SUPG_SPALART_ALLMARAS_TURBULENCE` (2.0/1.0). The inlet uses
`INCOMPRESSIBLE_NS_VELOCITY/PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE` and
`SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE` with select and upwind properties;
walls and inclusions use the same IWG types with zero values. Nitsche penalties are 10.0;
ghost penalties (`GHOST_NORMAL_FIELD`) are 5·10⁻³ for velocity, convection, viscosity, L2,
THETA and PHID, and 5·10⁻⁴ for pressure.

Pseudo-time continuation enters FEM through `TIME_CONTINUITY_DOF` IWGs on VX, VY and on
VISCOSITY. Their weight, `Func_Time_Weights`, is |v_prev| / (CFL · h_e) using the NLBGS
`PseudoTimeStep` from the logger, set to zero when that step is below 1 and during the
sensitivity analysis. `Func_Ramp_Property` multiplies the dynamic and kinematic viscosities by
up to 100 and returns linearly to the nominal value over the first 3 NLBGS iterations.

The wall distance is built in three stages: THETA (diffusion with time continuity, value 1.0
prescribed on fluid–inclusion and fluid–void boundaries), PHID (diffusion with a
THETA-dependent eigenstrain, value 0.0 prescribed on the same boundaries), and an `L2`
projection of log₁₀(1/φ²) computed from PHID. The SA model's `PropWallDistance` is
√(10^(−L2)), floored at 10⁻⁹.

Criteria: `IQIVolumePowDisp` (`POWER_DISSIPATION_BULK`, restricted by `Func_Select_Int_Domain`
to x ∈ [2, 12] sH, y ∈ [0, 10] sH), `IQIPerimeterItf` (`VOLUME` on the fluid–inclusion side
set), two `HEAT_METHOD_PENALTY` IQIs (fluid and inclusion), and `IQIFluidVolume` (`VOLUME`).
FEM computes its sensitivities by finite differences (`is_analytical_sensitivity = false`,
`POINT_3_CENTRAL`, perturbation 10⁻⁵).

**Solver (SOL).** One Amesos linear solver (`Amesos_Mumps` or `Amesos_Superludist`). The deck
documents six tested configurations in a comment table; the active one (case 4, "tested via
ctest") is:

| Setting | Value |
|---|---|
| `tSolverType` / `tFluidViscSolver` | `NLBGS` / `Staggered`: flow (VX, VY, P) and VISCOSITY are separate NLBGS sub-solvers |
| Sub-Newton | `InvResNormAdaptive` relaxation (1.0, damping 0.5), relative drop 0.25, at most 15 iterations |
| Flow NLBGS | at most 35 iterations, relative residual 10⁻⁹, reference iteration 4 |
| Pseudo-time control | `SolverPseudoTimeControlType::Expur`: 3 initial steps, initial and steady step 10⁶, at most 35 steps, maximum step 10⁶, residual drop 10⁻⁹, update threshold 10, residual factor 1.8, exponent 0.33 |
| Load control | off |

The forward problem is an outer `NLBGS_SOLVER` (one pass) over THETA, PHID, L2 and the flow
NLBGS. The sensitivity problem uses an outer NLBGS over THETA, PHID, L2 and a monolithic
single-step Newton solve of VX, VY, P, VISCOSITY. Initial values: VX = 0.001, VISCOSITY = 3ν
(scaled), L2 = 1.0, others 0.

**Optimization (OPT).** With `tCheckSensitivities = true`, a `SWEEP` with
`finite_difference_type = all`, epsilon 10⁻⁵, and `finite_difference_adv_indices = 30`, writing
`uBend_2D_Optimization.hdf5`. The three constraints are criteria 0, 1 and 4 (power dissipation,
perimeter, fluid volume). The objective is power dissipation normalized by its value at the
first iteration (perimeter and heat-method weights are 0.0). With `tCheckSensitivities = false`
the deck instead uses GCMMA (step 0.01, penalty 50, `max_its = 0`) and a fluid-volume
constraint against 0.3 of the available volume.

**Output (VIS).** `uBend_2D_Optimization.exo` with nodal THETA, PHID, L2, VX, VY, P, VISCOSITY,
turbulent and effective viscosity, the level set, the wall distance, and projected level-set
and alpha fields from the heat-method IQIs, on the fluid and inclusion sets, time offset 10.

## What the test checks

The ctest runs on 1 processor (case 0) and 4 processors (case 1; rank 0's output). Both cases
write `uBend_2D_Optimization.exo.e-s.0001` and `uBend_2D_Optimization.hdf5`.

| Case | Procs | Nodes | Elements | Reference node |
|---|---|---|---|---|
| 0 | 1 | 2916 | 2607 | 1991 |
| 1 | 4 | 848 | 772 | 735 |

The test requires the counts above, the reference-node coordinates (6.99144, 3.486625) to a
relative error of 10⁻⁸, the time value 11.0, and VX/VY, P and VISCOSITY at that node to a
relative error of 10⁻⁴. It requires the three constraints in `constraints eval_1-1` to match
references (≈ 11.586, 21.208, 102.33) to 10⁻⁴, and checks (`CHECK`) the analytic constraint
gradients against the forward, backward and central finite-difference gradients to a relative
error of 10⁻³.

## Run it

```bash
ctest -R uBend_2D_Optimization --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest names are `uBend_2D_Optimization-1-procs` and `uBend_2D_Optimization-4-procs`.

## Notes

- The deck defines its own `gInterpolationOrder = 1` and does not read `gTestCaseIndex`; the
  test sets both globals only for the harness.
- `CMakeLists.txt` keeps `EXTRA_SO_INCLUDES ${MORIS_PACKAGE_DIR}/GEN/GEN_CORE/src`; its comment
  notes that this directory no longer exists.
- The comment next to the third `LINE` geometry says x = 0, but its centre is x = 2 sH.
- Related decks with the same wall-distance construction:
  {doc}`Thermal_Turbulent_Flow_About_Sphere_2D`.
