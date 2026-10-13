# Compressible channel, wall boundary terms only

A single-element 2D channel for the compressible Navier–Stokes formulation (perfect gas,
unknowns pressure, velocity and temperature). In its checked-in configuration the deck
assembles only the symmetric Nitsche no-slip condition on the channel walls; the bulk
residual and the inlet and outlet conditions are present but switched off by flags. The deck
writes the Jacobian, the LHS and the solution vectors to file, which marks it as a debugging
and regression case for the compressible boundary terms rather than a flow simulation. Its
transient, heated siblings are {doc}`Heated_Channel_2D <Heated_Channel_2D>` and
{doc}`Heated_Bubble_2D <Heated_Bubble_2D>`.

```{exa-summary} projects/EXA/fluid/compressible/Channel_2D_Compressible
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 0.1 rectangle, x ∈ [0, 1], y ∈ [−0.05, 0.05] |
| Phases | one fluid phase (phase 0) covering the whole domain |
| Material | perfect gas: c_v = 718, R = 287.06, μ = 1.716·10⁻³, κ = 24.35·10⁻³ |
| Walls | no-slip velocity (0, 0) on bottom and top (side ordinals 1, 3), symmetric Nitsche |
| Inlet / outlet | defined (inlet pressure (1 + 10⁻⁴)·10⁵, inlet temperature 300, outlet pressure 10⁵) but disabled |
| Bulk residual | disabled (`tHaveBulk = false`) |
| Initial state | P = 10⁵, VX = √(10⁻⁴·P/ρ₀) ≈ 2.93 with ρ₀ = P/(R·300) ≈ 1.16, VY = 0, TEMP = 300 |
| Time | static (`tIsPseudoTransient = false`; no time steps set) |

The deck's switches are file-scope booleans: `tHaveBulk`, `tHaveTopBottomVelNitsche`,
`tHavePressureInletBC`, `tHavePressureOutletBC`, `tHaveTemperatureInletBC`, `tAllowSlip`,
`tHaveGLS`, `tIsPseudoTransient`. Only `tHaveTopBottomVelNitsche` and `tHaveGLS` are `true`,
and `tHaveGLS` only acts on the disabled bulk IWG.

## Setup in MORIS

**Mesh (HMR).** 1 × 1 element, quadratic Lagrange and B-spline orders
(`tIpOrder = 2`), `domain_offset` (0, −0.05). `lagrange_to_bspline = "-1"`
(`tUseLagrange = true`), so the Lagrange basis carries the DOFs directly.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Func_Dummy_Plane`
(φ = y − 10000), which is negative everywhere in the domain. Every element therefore lies in
phase 0 and the fluid sets are `HMR_dummy_c_p0, HMR_dummy_n_p0`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, no ghost
stabilization; the enriched integration mesh is printed and written to Exodus.

**Physics (FEM).** Material model `PERFECT_GAS` (`MMFluid`, DOF dependencies P, TEMP) feeds the
constitutive model `FLUID_COMPRESSIBLE_NEWTONIAN` (`CMFluid`). Stabilization parameters:
`COMPRESSIBLE_DIRICHLET_NITSCHE` (`NitscheSP`) with penalties `0; γ; γ; γ` for P, VX, VY, TEMP,
where γ = 100·√(h_x² + h_y²) ≈ 100.5, scaled by the viscosity and conductivity properties
(`tAutoScaleNitschePenalty = true`); and a `DIRICHLET_NITSCHE` dummy (`DummySP`, value 1.0)
reserved for GLS. The active IWG is `COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE` on the
sideset with `side_ordinals` "1,3" and prescribed velocity `PropZeroU`. The disabled IWGs are
`COMPRESSIBLE_NS_BULK`, `COMPRESSIBLE_NS_DIRICHLET_UNSYMMETRIC_NITSCHE` for inlet pressure
and inlet temperature (side 4), and `COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE` with
pressure upwinding for the outlet (side 2). IQIs: `DOF` for P, VX, VY and TEMP, and `PROPERTY`
IQIs for the local Mach number and Reynolds number, computed in `Func_Mach_Number` and
`Func_Reynolds_Number` from the field values.

**Solver (SOL).** `AMESOS_IMPL` direct solve, Newton with `NLA_rel_res_norm_drop` 2·10⁻⁶,
relaxation 1.0 and at most 10 iterations on DOFs `P;VX,VY;TEMP`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). Debug output:
`DLA_LHS_output_filename` "LHS", `SOL_save_operator_to_matlab` "Channel_Compressible",
`TSA_Save_Sol_Vecs_to_file` "SolVec". No preconditioner.

**Output (VIS).** `Channel_2D_Compressible.exo` on the fluid sets with nodal fields
`P, VX, VY, TEMP, InletVelX, InletVelY, Ma, Re`, `Save_Frequency` 1, `Time_Offset` 10.0.

## What the test checks

The test case calls the workflow with `./Channel_2D_Compressible.so` and requires a return
code of 0 (`REQUIRE( tRet == 0 )`). It sets `gInterpolationOrder = 0` and
`gPrintReferenceValues = false`; no field values are compared.

## Run it

The example is not built: `add_subdirectory(Channel_2D_Compressible)` is commented out in
`projects/EXA/fluid/compressible/CMakeLists.txt`, so the ctest below is not registered. After
re-enabling that line, in a build with `MORIS_USE_EXAMPLES=ON`:

```bash
ctest -R Channel_2D_Compressible-1-procs --output-on-failure
```

## Notes

- VIS `IQI_Names` lists `IQIinletVelX` and `IQIinletVelY`, which the FEM block does not define;
  the other six IQI names match FEM IQIs
  ([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)).
- `Func_Mach_Number` and `Func_Reynolds_Number` read both velocity components from the `VX`
  field interpolator (the vector field VX, VY).
