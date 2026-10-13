# Channel with four cylinders, static

A steady, coupled incompressible-flow and heat-transport problem in a 2D channel obstructed by
four circular cylinders. The channel and the cylinders are immersed in a rectangular background
mesh as level-set geometries; velocity, pressure, and temperature are solved monolithically on
the fluid phase, with all boundary conditions on the cut boundaries imposed by Nitsche's method
and ghost stabilization on the cut elements. The deck exercises the SUPG/PSPG-stabilized
incompressible Navier–Stokes IWGs together with SUPG-stabilized advection–diffusion. The
diffusion-only variant is {doc}`Channel_with_Four_Cylinders_Static_Temp_Only <Channel_with_Four_Cylinders_Static_Temp_Only>`;
the transient variant is {doc}`Channel_with_Four_Cylinders_Transient <Channel_with_Four_Cylinders_Transient>`.

```{exa-summary} projects/EXA/thermal/advection/Channel_with_Four_Cylinders_Static
```

## Problem

| | |
|---|---|
| Background domain | 4.0 × 2.0 rectangle, offset (−1.24, −0.86) |
| Physical domain | channel 0 ≤ x ≤ 2.2, 0 ≤ y ≤ 0.41 |
| Cylinders | radius 0.05, centers (0.4, 0.1), (0.8, 0.3), (1.2, 0.1), (1.6, 0.3) |
| Fluid | viscosity 0.001, density 1.0, heat capacity 1.0, conductivity 0.00005 |
| Inlet (x = 0) | parabolic velocity u_x = 4 · 1.5 · y (0.41 − y) / 0.41², temperature 0.0 |
| Walls and cylinders | no-slip (zero velocity) |
| Cylinder surfaces | heat flux 2.0 |
| Outlet (x = 2.2) | no condition written in the deck |
| Analysis | static, one solve; initial guess VX = 1.0, VY = P = TEMP = 0 |

## Setup in MORIS

**Mesh (HMR).** 176 × 88 elements with one level of initial refinement
(`pattern_initial_refinement = 1`), refinement and staircase buffers of 1. Lagrange and B-spline
orders are both set from the global `gInterpolationOrder`, which the test sets to 1.

**Geometry and phases (GEN).** Eight level-set geometries: four `LINE` fields for the bottom,
top, left, and right channel edges and four `CIRCLE` fields for the cylinders, each refined once
on mesh 0. A user phase function `get_phase_index` maps the eight geometry signs to five phases:
0 fluid, 1 void inlet, 2 cylinder, 3 void walls, 4 void outlet (commented "not used"). The GEN field output is written to
`Channel_with_Four_Cylinders_Static_GEN.exo`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on (`isGhost = true`), and Exodus output of the integration mesh.

**Physics (FEM).** The block is written against named FEM phases (`PhaseFluid`, `PhaseInlet`,
`PhaseCylinder`, `PhaseWall`) and `neighbor_phases`, so no mesh-set names appear in the deck;
MORIS derives them from the phase indices
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
Constitutive models are `FLUID_INCOMPRESSIBLE` (VX, VY; P) and `DIFF_LIN_ISO` (TEMP). Bulk IWGs:
`INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK` with the
`INCOMPRESSIBLE_FLOW` (SUPG/PSPG) parameter (`36.0/1.0`; the deck notes that it uses a time
step in τ although the problem is static), `SPATIALDIFF_BULK`, and `ADVECTION_BULK` with
`SUPG_ADVECTION` (`1.0`). Boundary IWGs on the fluid side of each interface:

| Boundary | Velocity | Pressure | Temperature |
|---|---|---|---|
| inlet | `INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE` | `INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE` | `SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` |
| walls, cylinders | same pair, zero velocity | | — |
| cylinders | | | `SPATIALDIFF_NEUMANN`, flux 2.0 |

Nitsche penalties are `VELOCITY_DIRICHLET_NITSCHE` `100.0/1.0` and `DIRICHLET_NITSCHE` `100.0`.
Four `GHOST_NORMAL_FIELD` IWGs on the fluid ghost facets use `VISCOUS_GHOST` (0.05),
`CONVECTIVE_GHOST` (0.05), `PRESSURE_GHOST` (`0.005/1.0`), and `GHOST_DISPL` (0.005) for
temperature. Four `DOF` IQIs report VX, VY, P, and TEMP.

**Solver (SOL).** Direct solve with `AMESOS_IMPL`, no preconditioner. One Newton solver over
`VX,VY;P;TEMP` with combined residual/Jacobian assembly, relative residual drop 10⁻⁴,
relaxation 0.8, at most 20 iterations, and the default single static time step
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Output (VIS).** `Channel_with_Four_Cylinders_Static.exo` on the fluid sets
`HMR_dummy_n_p0,HMR_dummy_c_p0`: nodal VX, VY, P, TEMP and the four IQIs as global variables.

## What the test checks

Two Catch2 test cases. `Channel_with_Four_Cylinders_Static_Linear` runs the deck with
`gInterpolationOrder = 1` and reads the Exodus output:

- Serial: 2 dimensions, 16887 nodes, 15533 elements; reference node 3307 at (0.76, 0.14).
- 4 processors: rank 1 only checks 3490 nodes, 3244 elements, and reference node 2628.
- In both cases: node coordinates to 10⁻⁸, time value 1.0 to 10⁻⁸, and, with relative
  tolerance 10⁻⁵ (marked FIXME as loose because of serial/parallel differences), the nodal values
  VX ≈ 1.5935, VY ≈ −0.38320, P ≈ 1.3582, TEMP ≈ 1.3719 and the global IQIs
  VX ≈ 0.90158, VY ≈ −0.015818, P ≈ 0.76451, TEMP ≈ 13.729.

`Channel_with_Four_Cylinders_Static_Quadratic` sets `gInterpolationOrder = 2`, but its body is
wrapped in `if (false)` ("quadratic case currently not working") and runs nothing.

## Run it

```bash
ctest -R Channel_with_Four_Cylinders_Static-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Channel_with_Four_Cylinders_Static-4-procs --output-on-failure
```
