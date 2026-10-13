# Channel 2D static

Steady incompressible Navier–Stokes flow through a 2D channel immersed in a square background
mesh. Four level-set planes cut the channel out of the background grid; a parabolic velocity
profile enters through the inlet plane and the walls carry no-slip conditions, all imposed
weakly with Nitsche's method. The deck demonstrates immersed incompressible flow with SUPG/PSPG
stabilization and three ghost-penalty terms on the cut elements.

```{exa-summary} projects/EXA/fluid/laminar/Channel_2D_Static
```

## Problem

| | |
|---|---|
| Channel | x ∈ [−2.5, 2.5], y ∈ [−0.5, 0.5] (half-height r = 0.5) |
| Background domain | 10 × 10 square with offset (−4.63, −4.74) |
| Phases | phase 10 of the four-plane phase table is the fluid |
| Material | density 1.0, viscosity 1.0 |
| Inlet | parabolic velocity u_x = −(y − r)(y + r)/(2r²), peak 0.5 at y = 0, u_y = 0 |
| Walls | zero velocity |
| Time | steady; Newton iteration from VX = 0.5, VY = 0, P = 0 |

The inlet is the interface between phase 10 and phase 8 (`iside_b0_10_b1_8`); the walls are
the interfaces to phases 2 and 14. The fourth interface of phase 10 (to phase 11) is not
referenced by any IWG.

## Setup in MORIS

**Mesh (HMR).** 200 × 200 elements on the 10 × 10 domain (element edge 0.05), linear
Lagrange and B-spline meshes (`lagrange_to_bspline = "0"`), refinement and staircase buffers
of 3, no initial or adaptive refinement. `processor_decomposition_method` 1 with
`processor_dimensions` (2, 1).

**Geometry and phases (GEN).** Four `USER_DEFINED` level sets with the callback `Func_Plane`
φ = n_x (x − x_c) + n_y (y − y_c), parameters given as `variable_1`…`variable_4`:
bottom (normal (0, 1), y_c = −0.5), top (normal (0, 1), y_c = 0.5), left (normal (1, 0),
x_c = −2.5) and right (normal (1, 0), x_c = 2.5). No phase table is set, so the default table
over the four geometries applies; the deck takes phase 10 as the fluid
(see {doc}`../../../user-guide/modules`).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, `ghost_stab`
on, which produces the `ghost_p10` set used by the ghost IWGs
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
Integration-mesh output is off.

**Physics (FEM).** The block uses `hack_for_legacy_fem()` and explicit `mesh_set_names`.
Constitutive model `FLUID_INCOMPRESSIBLE` (`CMFluid`). Stabilization parameters:
`INCOMPRESSIBLE_FLOW` (`SPSUPGNS`, parameter 36.0), `VELOCITY_DIRICHLET_NITSCHE`
(`SPNitscheU`, "50.0/1.0"), `VISCOUS_GHOST` (0.05), `CONVECTIVE_GHOST` (0.05) and
`PRESSURE_GHOST` ("0.005/1.0"). IWGs:

- `INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK` on
  `HMR_dummy_c_p10, HMR_dummy_n_p10`;
- `INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE` and
  `INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE` on the inlet (property
  `PropInletU`, callback `Func_Inlet_U`) and on the walls (`PropZeroU`);
- three `GHOST_NORMAL_FIELD` IWGs (viscous, convective, pressure) on `ghost_p10`.

IQIs: `DOF` for VX, VY and P on the fluid sets.

**Solver (SOL).** `AMESOS_IMPL` direct solve; Newton with `NLA_rel_res_norm_drop` 10⁻⁶,
relaxation 1.0, at most 50 iterations, on `VX,VY;P`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). No
preconditioner.

**Output (VIS).** `Channel_2D_Static.exo` on the fluid sets with nodal VX, VY, P.

## What the test checks

The test file defines two Catch2 test cases, `Channel_2D_Static_Inlet_Velocity` and
`Channel_2D_Static_Inlet_Pressure`. Both set `gInletVelocityBCFlag = true` and
`gInletPressureBCFlag = false`, call the workflow with `./Channel_2D_Static.so`, and make no
assertion: the return code is not checked (`REQUIRE( tRet == 0 )` is commented out). A comment
in the test file records that the legacy test never toggled the flags between the two cases,
so the "Inlet_Pressure" case runs the same velocity-inlet configuration. The example is
registered for 1 and 2 processors.

## Run it

```bash
ctest -R Channel_2D_Static --output-on-failure   # runs Channel_2D_Static-1-procs and -2-procs
```

The deck reads two globals set by the test harness: `gInletVelocityBCFlag` (adds the inlet
Nitsche IWGs) and `gInletPressureBCFlag` (adds an `INCOMPRESSIBLE_NS_IMPOSED_PRESSURE` IWG with
inlet pressure 20.0 on the inlet set).

## Notes

- The `.dox` page (`EXA_Channel_2D_Static`) contains section headings only.
- The pressure-inlet path (`gInletPressureBCFlag = true`) is not exercised by any test case.
