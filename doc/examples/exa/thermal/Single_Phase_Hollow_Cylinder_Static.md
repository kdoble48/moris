# Steady conduction in a hollow cylinder (2D annulus)

A steady 2D heat-conduction problem in an annulus defined by two circular level sets on a
square background mesh. The temperature is prescribed on the inner circle and a heat flux on the
outer circle, both on immersed boundaries. The deck codes the radial analytic solution and
evaluates the L2 and H1 semi-norm errors, the maximum temperature and the volume. It exercises
single-phase physics on a region bounded by two immersed curves, with Nitsche-enforced
Dirichlet conditions, a Neumann condition on an immersed boundary, ghost stabilization, and
adaptive refinement at the interfaces.

```{exa-summary} projects/EXA/thermal/diffusion/Single_Phase_Hollow_Cylinder_Static
```

## Problem

| | |
|---|---|
| Background domain | 5.0 × 5.0 square, offset (−2.37, −2.32) |
| Physical domain | annulus 1.0 ≤ r ≤ 2.0 around (0, 0), phase 1 |
| Conductivity | κ = 1.0 |
| Inner circle (r = 1) | T = 5.0, Nitsche-enforced |
| Outer circle (r = 2) | heat flux q = 20.0 (Neumann) |
| Source | none |
| Time | static, default time-solver settings |

`AnalyticTemperatureFunc` evaluates T(r) = T_in + Q ln(r / R_in) / (2πκ) with
Q = q · 2π R_out, T_in = 5.0, R_in = 1.0 and R_out = 2.0; `AnalyticdTemperaturedxFunc` returns
the corresponding gradient.

## Setup in MORIS

**Mesh (HMR).** 21 × 21 elements of order `gInterpolationOrder` (set to 1 by the test) for both
the Lagrange and the B-spline mesh, refinement and staircase buffer 3, and
`adaptive_refinement_level` = 2. With the offset of 0.13 and 0.18 from a centred domain, the
origin does not lie on a grid line of the 21 × 21 background mesh (element size 5/21).

**Geometry and phases (GEN).** Two `CIRCLE` level sets; the deck sets only their radii, 2.0 and
1.0, while the analytic solution assumes the centre (0, 0). FEM uses phase 1 for the annulus,
with the default phase table.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, Exodus output of the integration mesh.

**Physics (FEM).** Legacy set-name input
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
`DIFF_LIN_ISO` (conductivity only) with `SPATIALDIFF_BULK` on `HMR_dummy_c_p1,HMR_dummy_n_p1`;
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE` (penalty 100.0) on
`iside_b0_1_b1_0`; `SPATIALDIFF_NEUMANN` with the flux 20.0 on `iside_b0_1_b1_3`;
`GHOST_NORMAL_FIELD` with `GHOST_DISPL` (0.001) on `ghost_p1`. IQIs: `DOF` temperature,
`PROPERTY` analytic temperature, `L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC`, `MAX_DOF`
(function parameters `5.0/30.0`) and `VOLUME`.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type, Newton with relative residual drop
10⁻⁷, otherwise default nonlinear and time-solver settings.

**Output (VIS).** `Single_Phase_Hollow_Cylinder_Static.exo` on the phase-1 bulk sets: nodal
`TEMP` and `TEMP_ANALYTIC`, global `L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC`, `MAX_DOF` and
`VOLUME`.

## What the test checks

The Catch2 case `Single_Phase_Hollow_Cylinder_Static` sets `gInterpolationOrder = 1` and calls
the workflow. It neither checks the return code nor reads the output, so the test fails only if
the run aborts.

## Run it

```bash
ctest -R Single_Phase_Hollow_Cylinder_Static-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The `.dox` file contains only empty section headings.
- `MSIParameterList` is empty, so the DOF-to-mesh assignment uses defaults.
