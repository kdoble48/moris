# Heated compressible channel (transient)

A transient compressible Navier–Stokes simulation of a long 2D channel of air at rest,
heated by a volumetric heat source concentrated at mid-length. Pressure is held at both
channel ends and the side walls allow only tangential flow. The deck exercises the
compressible bulk IWG with a body heat load, Nitsche boundary terms for the compressible
system, and time-continuity IWGs for space–time time stepping. The radially heated
counterpart is {doc}`Heated_Bubble_2D <Heated_Bubble_2D>`; the steady boundary-term case is
{doc}`Channel_2D_Compressible <Channel_2D_Compressible>`.

```{exa-summary} projects/EXA/fluid/compressible/Heated_Channel_2D
```

## Problem

| | |
|---|---|
| Domain | 100 × 1 rectangle, x ∈ [0, 100], y ∈ [0, 1] |
| Phases | one fluid phase (phase 0) covering the whole domain |
| Material | perfect gas: c_v = 718, R = 287.058, μ = 1.716·10⁻⁵, κ = 24.35·10⁻³ |
| Heat load | Q(x) = 10⁵ · exp(−200 (2x/L − 1)²), L = 100 |
| Side walls | zero normal velocity on bottom and top (side ordinals 1, 3) |
| Channel ends | pressure fixed to the initial pressure on sides 2 and 4, with pressure upwinding (weight 1.0) |
| Initial state | ρ₀ = 1.0, T₀ = 273, P₀ = R ρ₀ T₀ ≈ 78 367, velocity 0 |
| Time | 10 time steps of Δt = 10/c_s ≈ 0.0302, total ≈ 0.302 |

The speed of sound is computed in the deck as c_s = √(γ R T₀) ≈ 331.2 with
γ = 1 + R/c_v ≈ 1.400; the step size is `20 L / 100 / c_s / 2`. The HMR callback prints the
step size and the time-continuity weight at start-up.

## Setup in MORIS

**Mesh (HMR).** 100 × 1 elements with one quadratic Lagrange mesh (`tIpOrder = 2`). No
B-spline mesh is added because `tUseLagrange = true`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Func_Dummy_Plane`
(φ = y − 10000), which places the whole domain in phase 0; the fluid sets are
`HMR_dummy_c_p0, HMR_dummy_n_p0`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, no ghost
stabilization; the enriched integration mesh is printed and written to Exodus.

**Physics (FEM).** Material model `PERFECT_GAS` and constitutive model
`FLUID_COMPRESSIBLE_NEWTONIAN`, as in the sibling decks. IWGs:

- `COMPRESSIBLE_NS_BULK` on `P;VX,VY;TEMP` with `PropHeatLoad` as `BodyHeatLoad`
  (`Func_Heat_Load_Distribution`); GLS is off (`tHaveGLS = false`).
- `COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE` on sides "1,3" with `PropZeroU` as the
  prescribed velocity and `PropSelectY` (matrix [0 0; 0 1]) as `SelectVelocity`, so only the
  y component is constrained.
- `COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE` on sides "2,4" with `PropInitialPressure` as
  `PrescribedDof1` and `PropUpwind` (1.0) as `PressureUpwind`.
- Three `TIME_CONTINUITY_DOF` IWGs for P, VX/VY and TEMP, each with weights
  `tTimePenalty / Δt` = 0.1/Δt ≈ 3.31 and the initial-condition properties.

The Nitsche parameter `COMPRESSIBLE_DIRICHLET_NITSCHE` uses penalties 0, 100, 100, 0 for P, VX,
VY, TEMP, scaled by viscosity and conductivity. IQIs: `DOF` for P, VX, VY, TEMP and `PROPERTY`
IQIs for the Mach number, the Reynolds number (length scale = channel height) and the heat
load.

**Solver (SOL).** `AMESOS_IMPL` direct solve, Newton with `NLA_rel_res_norm_drop` 10⁻⁹,
relaxation 1.0, at most 10 iterations. The time solver runs `TSA_Num_Time_Steps` 10 over
`TSA_Time_Frame` ≈ 0.302 with `TSA_time_level_per_type` "P,2;VX,2;VY,2;TEMP,2" and initial
values P = P₀, VX = VY = 0, TEMP = 273
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). No
preconditioner.

**Output (VIS).** `Heated_Channel_2D.exo` on the fluid sets with nodal fields
`P, VX, VY, TEMP, Ma, Re, Q`, written every step (`Save_Frequency` 1, `Time_Offset` 10.0).

## What the test checks

The test case calls the workflow with `./Heated_Channel_2D.so` and requires a return code of
0 (`REQUIRE( tRet == 0 )`). It sets `gInterpolationOrder = 0` and
`gPrintReferenceValues = false`; no field values are compared.

## Run it

```bash
ctest -R Heated_Channel_2D-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- A `COMPRESSIBLE_NS_BOUNDARY` outlet IWG is present only as a commented-out block.
- Setting `tHaveFixedEnds = true` replaces the end-pressure condition by zero velocity on sides
  2 and 4. Debug switches `tWriteJacToMatlab`, `tWriteSolVecToHDF5`, `tWriteLhsToHDF5` and
  `tOrderAdofsByHost` are all `false`.
