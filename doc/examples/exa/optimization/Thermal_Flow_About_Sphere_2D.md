# Thermal flow about a sphere (2D)

Conjugate heat transfer in a 2D channel: incompressible laminar flow passes a heated circular
solid inclusion, and the fluid carries the heat downstream by advection. The inclusion radius is
the single ADV, and the deck runs a `SWEEP` that compares adjoint sensitivities of nine flow and
thermal criteria against finite differences. It tests the immersed, stabilized
Navier–Stokes/advection–diffusion formulation, a staged NLBGS forward solve with a monolithic
sensitivity solve, and phase assignment through a user phase function instead of a phase table.

```{exa-summary} projects/EXA/optimization/Thermal_Flow_About_Sphere_2D
```

## Problem

All inputs are given in SI units and scaled to a nondimensional system with reference length
1 mm, reference pressure equal to the inlet pressure, and reference density 0.3289 kg/m³.

| | |
|---|---|
| Channel | 5.0 mm × 1.0 mm (fluid region, shifted by 0.02 mm from the background origin) |
| Inclusion | circle at (0.75, 0.5) mm, radius 0.225 mm |
| Phases | 0 fluid, 1 solid, 2 inlet void, 3 outlet void, 4 lateral void |
| Fluid | air at 800 °C: μ = 43.32·10⁻⁶ N s/m², ρ = 0.3289 kg/m³, c_p = 1210 J/(kg K), k = 76.26·10⁻³ W/(m K) |
| Solid | ρ = 7000 kg/m³, c_p = 200 J/(kg K), k = 150 W/(m K) |
| Inlet | imposed pressure 0.001 × 101325 Pa = 101.325 Pa; temperature 0.0 |
| Outlet | no condition imposed |
| Walls | zero velocity on the lateral channel walls and on the inclusion surface (Nitsche) |
| Heat source | volumetric load 250·10⁴ (deck units W/m², scaled as a volumetric load) in the solid |
| Analysis | steady, one time step |

The reference velocity is √(101.325 / 0.3289) ≈ 17.55 m/s, which gives
Re = 10⁻³ · 17.55 · 0.3289 / 43.32·10⁻⁶ ≈ 133. The deck prints this value and its
nondimensional equivalent, 1/μ*.

## Setup in MORIS

**Mesh (HMR).** A 5.25 mm × 1.25 mm background domain offset by (−0.125, −0.125) mm, with
`ceil(5.25/0.2) = 27` × `ceil(1.25/0.2) = 7` elements and no initial refinement. Two B-spline
meshes of order `gInterpolationOrder` (pattern 0 for the level sets, pattern 1 for the state
fields) and one Lagrange mesh of the same order; refinement and staircase buffers are 1.
MSI places VX, VY and TEMP on B-spline mesh 1.

**Geometry and phases (GEN).** Five `USER_DEFINED` level-set fields: `Func_Sphere` for the
inclusion and four instances of `Func_Plane` for the inlet, outlet, lower and upper planes, each
shifted by 0.1 element lengths (0.02 mm) off the background grid. The sphere radius is a
`Design_Variable` with bounds 0.9r and 1.1r around r = 0.225; `Func_Sphere_Deriv` supplies
the field sensitivity (∂φ/∂r = 1). Instead of a phase table, `phase_function_name =
get_phase_index` maps the five geometry signs to five phases, and FEM refers to them by name
(`PhaseFluid`, `PhaseSolid`, `PhaseVoidFront`, `PhaseVoidBack`, `PhaseVoidLateral`, `PhaseAll`).
`IQI_types` lists nine criteria, in the order the OPT callbacks index them
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on, `high_to_low_dbl_side_sets`, and Exodus output of the integration mesh.

**Physics (FEM).** The fluid uses `FLUID_INCOMPRESSIBLE` with
`INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK` (pressure spring 10⁻⁶),
stabilized by `INCOMPRESSIBLE_FLOW` (SUPG/PSPG, C_inv = 36 for linear and 60 for quadratic
interpolation). Fluid temperature uses `DIFF_LIN_ISO` with `SPATIALDIFF_BULK` and
`ADVECTION_BULK` (`SUPG_ADVECTION`); the solid uses `SPATIALDIFF_BULK` with the volumetric heat
load. The fluid–solid interface is coupled with `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE`
(`NITSCHE_INTERFACE`, penalty 100). No-slip on the inclusion and the lateral voids uses the
`INCOMPRESSIBLE_NS_VELOCITY/PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE` pair
(`VELOCITY_DIRICHLET_NITSCHE`, penalty 100); the inlet uses `INCOMPRESSIBLE_NS_IMPOSED_PRESSURE`
and `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` for temperature. Ghost penalty
(`GHOST_NORMAL_FIELD`, order = interpolation order) is applied to velocity (`VISCOUS_GHOST`
and `CONVECTIVE_GHOST`, 0.05), pressure (`PRESSURE_GHOST`, 0.005) and fluid and solid
temperature (`GHOST_DISPL`, 0.05). FEM computes its own derivatives by `POINT_3_CENTRAL` finite
differences with perturbation 10⁻⁴.

The nine criteria are, in order:

| Index | IQI | Type |
|---|---|---|
| 0, 1 | `IQIInletThermalEnergy`, `IQIOutletThermalEnergy` | `THERMAL_ENERGY_CONVECTIVE_FLUX` on inlet/outlet |
| 2, 3 | `IQIInletTotalPressure`, `IQIOutletTotalPressure` | `TOTAL_PRESSURE` on inlet/outlet |
| 4 | `IQIPerimeterItf` | `VOLUME` on the solid–fluid side set |
| 5, 6 | `IQIInletMassFlow`, `IQIOutletMassFlow` | `MASS_FLOW` on inlet/outlet |
| 7 | `IQIMaxTemp` | `MAX_DOF` of TEMP in the solid, reference 250, exponent 8 |
| 8 | `IQISolidVolume` | `VOLUME` of the solid |

**Solver (SOL).** Two Amesos linear solvers (`Amesos_Mumps` when `MORIS_USE_MUMPS` is
defined, otherwise `Amesos_Superludist`). The forward problem is a one-way staged `NLBGS_SOLVER`
(one iteration): a Newton solve of VX, VY, P (relative residual drop 10⁻⁸, at most 20
iterations), then one Newton step for TEMP. The sensitivity problem uses a monolithic Newton
solver over VX, VY, P, TEMP. The initial guess is VX = 0.01, all other fields 0.

**Optimization (OPT).** `SWEEP` with `finite_difference_type = all` (forward, backward and
central), epsilon 10⁻⁴, one evaluation per ADV, results in `Thermal_Flow_About_Sphere_2D_<case>_SEN.hdf5`.
The objective is criterion 0; the nine constraints are the nine criteria, one to one. See the
{doc}`sensitivity-sweep step <../../../user-guide/workflow>` of the workflow guide.

**Output (VIS).** `Thermal_Flow_About_Sphere_2D_<case>.exo` with nodal VX, VY, P and TEMP on the
fluid and solid sets, time offset 10; GEN writes `Thermal_Flow_About_Sphere_2D_<case>_GEN.exo`.

## What the test checks

The deck reads `gInterpolationOrder` (B-spline and Lagrange order) and `gTestCaseIndex` (output
file suffix). The test runs two cases per processor count:

| Case | Procs | Order | Nodes | Elements |
|---|---|---|---|---|
| 0 | 1 | 1 | 654 | 443 |
| 1 | 1 | 2 | 1043 | 443 |
| 2 | 3 | 1 | 296 | 195 |
| 3 | 3 | 2 | 413 | 195 |

For each case it reads the Exodus file (on 3 processors, rank 0's file only) and requires: the
dimension, node and element counts above; the coordinates of a reference node (IDs 79, 613, 28,
85; all at (1.0417, 0.5893)) to a relative error of 10⁻⁸; the time value 11.0; and VX/VY,
P and TEMP at that node to a relative error of 10⁻⁴. It then requires the nine constraint
values from `constraints eval_1-1` in the HDF5 file to match stored references to 10⁻⁴, and
checks (with `CHECK`, not `REQUIRE`) that the analytic constraint gradients agree with each of
the forward, backward and central finite-difference gradients to a relative error of 10⁻³.

## Run it

```bash
ctest -R Thermal_Flow_About_Sphere_2D --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

This runs `Thermal_Flow_About_Sphere_2D-1-procs` and `Thermal_Flow_About_Sphere_2D-3-procs`.

## Notes

- The {doc}`3D variant <Thermal_Flow_About_Sphere_3D>` uses the same setup in a channel of
  square cross-section; {doc}`Thermal_Turbulent_Flow_About_Sphere_2D` adds a Spalart–Allmaras
  turbulence model.
- Several comments in the deck do not match the values (the inlet pressure is commented as
  0.6 atm and the inlet temperature as 500 K); the table above gives the values the deck uses.
