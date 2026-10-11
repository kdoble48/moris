# Thermal turbulent flow about a sphere (2D)

Conjugate heat transfer in a 2D channel with a Spalart–Allmaras turbulence model: incompressible
flow passes a heated circular solid inclusion, and the turbulent eddy viscosity enters both the
momentum and the heat equations. The wall distance that the turbulence model needs is computed
on the immersed geometry from two auxiliary diffusion problems (THETA and PHID). The inclusion
radius is the single ADV, and a `SWEEP` compares adjoint and finite-difference sensitivities of
eleven criteria. The deck extends {doc}`Thermal_Flow_About_Sphere_2D` and, like it, assigns
phases through a user phase function.

```{exa-summary} projects/EXA/optimization/Thermal_Turbulent_Flow_About_Sphere_2D
```

## Problem

Inputs are in SI units and scaled with reference length 1 mm, reference pressure equal to the
inlet pressure, and reference density 0.3289 kg/m³.

| | |
|---|---|
| Channel | 5.0 mm × 1.0 mm (fluid region, shifted by 0.01 mm from the background origin) |
| Inclusion | circle at (0.75, 0.5) mm, radius 0.225 mm |
| Phases | 0 fluid, 1 solid, 2 inlet void, 3 outlet void, 4 lateral void |
| Fluid | μ = 43.32·10⁻⁶ N s/m², ρ = 0.3289 kg/m³, ν = μ/ρ ≈ 1.317·10⁻⁴ m²/s, c_p = 1210 J/(kg K), Pr = 0.7, k = μ c_p / Pr ≈ 0.0749 W/(m K) |
| Solid | ρ = 7000 kg/m³, c_p = 200 J/(kg K), k = 150 W/(m K) |
| Inlet | imposed pressure 101.325 Pa; temperature 0.0; turbulent viscosity 3ν with a profile across the height (`Func_Inlet_V`) |
| Walls | zero velocity and zero turbulent viscosity on the lateral walls and the inclusion (Nitsche) |
| Heat source | volumetric load 250·10⁴ (deck units W/m², scaled as a volumetric load) in the solid |
| Analysis | steady, one time step |

The deck computes Re = L v_ref / ν = 10⁻³ · 17.55 / 1.317·10⁻⁴ ≈ 133. `Func_Inlet_V` scales the
inlet viscosity by sin(η) / √(0.05² + sin²(η)) with η = π (y − shift) / H, so it vanishes at
both walls.

## Setup in MORIS

**Mesh (HMR).** A 5.25 mm × 1.25 mm background domain offset by (−0.125, −0.125) mm, element
size 0.1 mm, giving 53 × 13 elements, no initial refinement. Two B-spline meshes of order
`gInterpolationOrder` (pattern 0 for the level sets, pattern 1 for the state fields) and one
Lagrange mesh. MSI places THETA, PHID, VX, VY, P, VISCOSITY and TEMP on B-spline mesh 1.

**Geometry and phases (GEN).** Five `USER_DEFINED` level sets, as in the laminar deck:
`Func_Sphere` (radius ADV, bounds 0.9r–1.1r around r = 0.225, sensitivity
`Func_Sphere_Deriv`) and four `Func_Plane` instances for inlet, outlet, lower and upper planes,
each 0.1 element lengths (0.01 mm) off the grid. `get_phase_index` maps the geometry signs to the
five phases, which FEM uses by name. `IQI_types` holds eleven criteria
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization, `high_to_low_dbl_side_sets`, and Exodus output of the integration mesh.

**Physics (FEM).** Flow uses `FLUID_INCOMPRESSIBLE_TURBULENCE_SPALART_ALLMARAS` (DOFs VX, VY, P,
VISCOSITY) with `INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK`
(pressure spring 0.0), stabilized by `INCOMPRESSIBLE_FLOW` (parameters `C_inv/0.0`, C_inv = 36
or 60 for linear or quadratic interpolation). The turbulent viscosity solves
`SPALART_ALLMARAS_TURBULENCE_BULK` with the `SPALART_ALLMARAS_TURBULENCE` CM (ft2 = 1.0,
alpha = 10.0) and `SUPG_SPALART_ALLMARAS_TURBULENCE` (exponent 2.0, source weight 1.0). Fluid
temperature uses `DIFF_LIN_ISO_TURBULENCE_SPALART_ALLMARAS` with `SPATIALDIFF_BULK` and
`ADVECTION_BULK` (`SUPG_ADVECTION`); the solid uses `DIFF_LIN_ISO` with the volumetric heat load.
Boundary and interface terms:

- fluid–solid temperature: `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, 100);
- no-slip on inclusion and lateral walls: `INCOMPRESSIBLE_NS_VELOCITY/PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE`
  (`VELOCITY_DIRICHLET_NITSCHE`, 100);
- zero turbulent viscosity on the same walls: `SPALART_ALLMARAS_TURBULENCE_DIRICHLET_UNSYMMETRIC_NITSCHE`
  (`TURBULENCE_DIRICHLET_NITSCHE`, 100);
- inlet: `INCOMPRESSIBLE_NS_IMPOSED_PRESSURE`, `SPALART_ALLMARAS_TURBULENCE_DIRICHLET_SYMMETRIC_NITSCHE`
  with `PropInletV`, and `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` for temperature;
- ghost penalty (`GHOST_NORMAL_FIELD`) on velocity (0.05), pressure (0.005), fluid and solid
  temperature (0.05), and VISCOSITY (0.05).

**Wall distance.** The property `PropWallDistance` returns the PHID field
(`Func_Wall_Distance`, with `Func_Wall_Distance_Der` as its DOF derivative) and feeds the
turbulence CM. THETA solves `SPATIALDIFF_BULK` with `TIME_CONTINUITY_DOF` terms (weights 10.0,
initial condition 0.0; capacity 0.1) and a prescribed value 1.0 on the inclusion surface (both
sides) and the lateral walls. PHID solves `SPATIALDIFF_BULK` with a THETA-dependent CM
(`EigenStrain` 1.0) and a prescribed value 0.0 on the same surfaces. Both use Nitsche penalty
100 and ghost penalty 0.01.

The eleven criteria, in order: `IQIOutletThermalEnergy` (`THERMAL_ENERGY_CONVECTIVE_FLUX`),
`IQIInletTotalPressure` and `IQIOutletTotalPressure` (`TOTAL_PRESSURE`), `IQIPerimeterItf`
(`VOLUME` on the solid–fluid side set), `IQIInletMassFlow` and `IQIOutletMassFlow`
(`MASS_FLOW`), `IQIMaxTemp` (`MAX_DOF`, reference 250, exponent 8), `IQISolidVolume`,
`IQIInletPowDisp` and `IQIOutletPowDisp` (`POWER_DISSIPATION` on inlet/outlet), and
`IQIVolumePowDisp` (`POWER_DISSIPATION_BULK`). FEM derivatives use `POINT_3_CENTRAL` finite
differences with perturbation 10⁻⁴.

**Solver (SOL).** One Amesos linear solver (`Amesos_Mumps` if `MORIS_USE_MUMPS`, else
`Amesos_Superludist`). The forward problem is an `NLBGS_SOLVER` (one pass) over four blocks in
order THETA, PHID, (VX, VY, P, VISCOSITY), TEMP. The flow block uses Newton with
`InvResNormAdaptive` relaxation (parameter 0.5, damping 0.5), relative residual drop 10⁻⁸ and at
most 100 iterations; the other blocks take one Newton step. The sensitivity solver is a second
`NLBGS_SOLVER` over the same blocks in which the flow block also takes one Newton step. The
initial guess is VX = 0.01 and VISCOSITY = 3ν (scaled); all other fields 0.

**Optimization (OPT).** `SWEEP` with `finite_difference_type = all`, epsilon 10⁻⁴, one
evaluation per ADV, results in `ChannelTurbulent2D_<case>_SEN.hdf5`. The objective is
criterion 0; the eleven constraints are the criteria, one to one.

**Output (VIS).** `ChannelTurbulent2D_<case>.exo` with nodal THETA, PHID, VX, VY, P, VISCOSITY,
TEMP, and the turbulent dynamic viscosity, effective viscosity and effective conductivity
(`TURBULENCE_FLUID_COEFFICIENT`, `TURBULENCE_HEAT_COEFFICIENT` IQIs), on the fluid and solid
sets, time offset 10.

## What the test checks

The deck file is `ChannelTurbulent2D.cpp`; it reads `gInterpolationOrder` and `gTestCaseIndex`.

| Case | Procs | Order | Nodes | Elements |
|---|---|---|---|---|
| 0 | 1 | 1 | 1440 | 1307 |
| 1 | 1 | 2 | 2785 | 1307 |
| 2 | 3 | 1 | 637 | 544 |
| 3 | 3 | 2 | 1058 | 544 |

For each case, from the Exodus file (rank 0's file on 3 processors), the test requires the
dimension, node and element counts above; the coordinates of a reference node (IDs 350, 1266,
119, 408; all at (1.6580, 0.7404)) to a relative error of 10⁻⁸; the time value 11.0; and THETA
and PHID, VX and VY, P, VISCOSITY and TEMP at that node to a relative error of 10⁻⁴. It then
requires the eleven constraint values in `constraints eval_1-1` to match stored references to
10⁻⁴, and checks (`CHECK`) the analytic constraint gradients against the central
finite-difference gradients: relative error below 10⁻³ when the finite-difference value exceeds
10⁻⁶, otherwise an absolute margin of 10⁻⁸.

## Run it

```bash
ctest -R ChannelTurbulent2D --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest names are `ChannelTurbulent2D-1-procs` and `ChannelTurbulent2D-3-procs` (`TEST_BASE`
differs from the directory name).

## Notes

- `PropFluidTurbPrandtl` is set from `prandtl` (0.7); the variable `turbPr = 0.9` is defined but
  not used.
- The comment above the sensitivity solver says "monolithic", but the solver it selects
  (index 6) is an `NLBGS_SOLVER`.
- In the gradient check, the branch condition and the small-value comparison index the matrices
  with the ADV index alone (`tConstraintsFD( tADVIndex )`), not with the constraint and ADV
  indices.
