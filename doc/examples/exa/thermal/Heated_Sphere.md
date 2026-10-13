# Heated sphere in a conducting cube (3D)

A steady 3D heat-conduction problem: a sphere with a uniform heat source sits inside a cube of
lower-conductivity material, and the exact temperature is imposed on the six outer faces. The
deck provides an analytic temperature field and its gradient, and evaluates the L2 and H1
semi-norm errors, the sphere surface area, the heat flux through the sphere surface from each
side, and the flux jump across the interface. It exercises a two-material immersed interface
with symmetric Nitsche coupling and ghost stabilization in 3D.

```{exa-summary} projects/EXA/thermal/diffusion/Heated_Sphere
```

## Problem

| | |
|---|---|
| Domain | cube [−1, 1]³ |
| Phases | 1 inner (sphere), 0 outer (cube minus sphere) |
| Sphere | centre (0, 0, 0), radius 0.5 |
| Conductivity | 1.0 inner, 0.125 outer |
| Density, heat capacity | 0.0 in both phases |
| Volumetric heat source | 1.0 inner, 0.0 outer |
| Outer faces | exact temperature imposed on `SideSet_1` … `SideSet_6`, Nitsche-enforced |
| Time | one time step, time frame 1.0 |

The exact solution coded in `Func_Exact_Temperature`, with r the distance from the centre, is

- T = 3/8 − r²/6 for r ≤ 0.5,
- T = (1/3)(1/r − 1) for r > 0.5,

and `Func_Exact_TemperatureGradient` returns its spatial gradient.

## Setup in MORIS

**Mesh (HMR).** 8 × 8 × 8 elements on the cube (`domain_offset` = −1 in each direction), one
level of initial refinement, refinement and staircase buffer 1. One Lagrange and one B-spline
mesh, both of order `gInterpolationOrder` (set to 1 by the test).

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Inclusion` (φ = 0.5 − r, values
with |φ| < 10⁻⁸ replaced by 10⁻⁸), with the default phase assignment.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets`.

**Physics (FEM).** Legacy set-name input
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)):
`DIFF_LIN_ISO` per phase with `SPATIALDIFF_BULK` IWGs on `HMR_dummy_{n,c}_p1` and
`HMR_dummy_{n,c}_p0`; `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE`
(penalty 100.0) on `SideSet_1_n_p0` … `SideSet_6_n_p0`; `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE`
with `NITSCHE_INTERFACE` (penalty 100.0) on `dbl_iside_p0_1_p1_0`; and `GHOST_NORMAL_FIELD` with
`GHOST_DISPL` (0.01) on `ghost_p1` and `ghost_p0`. The time-continuity IWG is present only as a
commented-out block. IQIs: `DOF` temperature, `PROPERTY` analytic temperature,
`L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC`, `VOLUME` of the domain and of the sphere surface
(`iside_b0_1_b1_0`), `TRACTION` heat flux on `iside_b0_1_b1_0` (inner CM) and
`iside_b0_0_b1_1` (outer CM), and `JUMP_TRACTION` on the double-sided interface.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type, Newton with at most 2 iterations
and relative residual drop 10⁻⁸, one time step over a time frame of 1.0.

**Output (VIS).** `Heated_Sphere.exo` on the bulk, outer-boundary and interface sets: nodal
`TEMP` and `TEMP_ANALYTIC`; global L2 and H1 errors, volume, sphere surface area and total
surface heat fluxes; faceted fields for the elemental surface area, the two surface heat fluxes
and the flux jump (`FACETED_INT`, `FACETED_AVG`).

## What the test checks

The Catch2 case `Heated_Sphere_Linear` sets `gInterpolationOrder = 1` and calls the workflow. It
does not check the return code, and the calls to `check_results` are commented out for both the
1- and 4-processor runs. The test therefore fails only if the run aborts. The disabled
`check_results` function and a commented-out `Heated_Sphere_Quadratic` case remain in the file;
the reference node and mesh counts in `check_results` are placeholders (0).

## Run it

```bash
ctest -R Heated_Sphere-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Heated_Sphere-4-procs --output-on-failure
```

## Notes

- The two ghost IWGs share the name `IWGGPInnerTemp`.
- The `.dox` documentation file is absent from this directory.
- A sibling problem with a 1D two-material analytic solution is
  {doc}`Two_Material_Bar <Two_Material_Bar>`.
