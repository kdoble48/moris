# Two channels with a separation wall, transient

A conjugate heat-transfer problem in 2D: two parallel fluid channels with counter-flowing
parabolic inlets at temperatures 1.0 and 0.0 are separated by a solid wall. Incompressible flow is solved in both channels, diffusion in the wall, and
the fluid–solid temperature coupling is a Nitsche interface condition on the immersed interfaces.
The deck is written in the pre-modern parameter-list API and is currently disconnected from the
build (see Notes below).

```{exa-summary} projects/EXA/thermal/advection/Two_Channels_with_Separation_Wall_Transient
```

## Problem

| | |
|---|---|
| Background domain | 3.0 × 1.0 rectangle, offset (−0.835, −0.186) |
| Physical domain | 0 ≤ x ≤ 1.845, 0 ≤ y ≤ 0.41 |
| Wall | 0.18 ≤ y ≤ 0.23 (mid-height 0.205, half-thickness 0.025) |
| Channels | bottom 0 ≤ y ≤ 0.18 and top 0.23 ≤ y ≤ 0.41, each of half-width 0.09 |
| Fluid | viscosity 1.0, density 1.0, conductivity 500.0, heat capacity 0.001 |
| Solid | density 2.0, conductivity 1.0, heat capacity 0.001 |
| Inlet velocity | u_x = −d (y − (y_c + r))(y − (y_c − r)) / (2 r²), peak magnitude 0.5; d = −1 (top), +1 (bottom) |
| Inlet temperature | 1.0 (top), 0.0 (bottom) |
| Channel walls | no-slip |
| Time | 20 time steps over a time frame of 0.4; initial VX = VY = 10⁻⁴, P = TEMP = 0 |

## Setup in MORIS

**Mesh (HMR).** 62 × 31 elements, linear Lagrange and B-spline meshes, refinement and staircase
buffers of 3, `adaptive_refinement_level = 0`.

**Geometry and phases (GEN).** Six `USER_DEFINED` planes through one function `Func_Plane`,
parameterized by `constant_parameters` (n_x, x_c, n_y, y_c): y = 0, 0.41, 0.18, 0.23 and
x = 0, 1.845. Phases are the default sign-pattern indices; the deck's fluid sets are `p34` and
`p46` and its solid set is `p42`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`, enriched integration mesh printed.

**Physics (FEM).** Constitutive models `FLUID_INCOMPRESSIBLE`, `DIFF_LIN_ISO` for the fluid, and
`DIFF_LIN_ISO` for the solid. In the fluid: `INCOMPRESSIBLE_NS_VELOCITY_BULK` and
`INCOMPRESSIBLE_NS_PRESSURE_BULK` with `INCOMPRESSIBLE_FLOW` (36.0), `SPATIALDIFF_BULK`, and
`ADVECTION_BULK` with `SUPG_ADVECTION`; in the solid: `SPATIALDIFF_BULK`. Inlets and channel
walls use `INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE` /
`INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE` (`VELOCITY_DIRICHLET_NITSCHE`
`100.0/1.0`) and `SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE` 100.0) for the
inlet temperatures. The fluid–solid coupling is `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` with
`NITSCHE_INTERFACE` (100.0) on `dbl_iside_p0_34_p1_42,dbl_iside_p0_46_p1_42`. Ghost IWGs
(`GHOST_NORMAL_FIELD`) use `VISCOUS_GHOST` 0.05, `CONVECTIVE_GHOST` 0.05, `PRESSURE_GHOST`
`0.05/1.0`, and `GHOST_DISPL` 0.05 for fluid and solid temperature. Four `DOF` IQIs report VX,
VY, P, TEMP. All sets are named explicitly
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
The deck contains no `TIME_CONTINUITY_DOF` IWG.

**Solver (SOL).** `AMESOS_IMPL`, one monolithic nonlinear solver over `VX,VY;P;TEMP` with
relative residual drop 10⁻⁷, 20 time steps over 0.4, `TSA_time_level_per_type` 2 for every
DOF type.

**Output (VIS).** `2_Channels_Temp.exo` on the fluid and solid sets with nodal VX, VY, P, TEMP.

## What the test checks

The Catch2 case loads `./bin/Two_Channels_with_Separation_Wall_Transient.so` and calls the
workflow. The return-code check is disabled (`(void) tRet;`) and all Exodus comparisons are
commented out, so the case asserts nothing about the result.

## Run it

The leaf registers `Two_Channels_with_Separation_Wall_Transient` (one process,
`NO_PROCS_SUFFIX`), but `projects/EXA/thermal/advection/CMakeLists.txt` has the
`add_subdirectory` line commented out, so the test does not exist in a build. After the deck is
ported and the line is restored:

```bash
ctest -R Two_Channels_with_Separation_Wall_Transient --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The parent CMakeLists states: "these decks use the pre-modern ParameterList API and no longer
  compile". The module functions take `Vector< Vector< ParameterList > >&`.
- Lines 235 and 241 of the deck (the inlet-velocity `function_parameters`) lack a `+` between
  `"/"` and the following `std::to_string(...)`.
- The test loads the `.so` from `./bin/`, unlike the `./` path used by the other thermal tests.
