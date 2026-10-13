# Bar with linear displacement and quadratic temperature interpolation

A 3D bar carries a structural body load and a heat source, and the two fields are discretized on
different B-spline meshes: displacements on a linear B-spline basis and temperature on a quadratic
one, both mapped to one quadratic Lagrange mesh. A planar level set splits the bar into two phases
with identical properties, so the deck also exercises Nitsche interface coupling of both fields
across an internal cut. The deck is named `Bar_Lin_Disp_Quad_Diff`; the example tests
per-DOF-type interpolation order selection through MSI.

```{exa-summary} projects/EXA/structure/linear/Thermoelastic_Varying_Order
```

## Problem

| | |
|---|---|
| Domain | 2.0 × 0.25 × 0.25 box |
| Phases | `p0` for x < 1.2, `p1` for x > 1.2 (plane level set x − 1.2), same material in both |
| Structure | E = 1.0, ν = 0.0, full 3D model; body load (1, 0, 0) |
| Thermal | conductivity 1.0; heat source 1.0 |
| Supports | `SideSet_4_n_p0` (x = 0 face): displacement (0, 0, 0) and temperature 1.0, both Nitsche |
| Interface | `dbl_iside_p0_0_p1_1`, Nitsche continuity of displacement and temperature |
| Analysis | one time step over a time frame of 1.0 |

`tCTE` is set to `"0.0e-5"` and no CTE property is passed to the structural constitutive models,
so in this deck the displacement solution does not depend on the temperature.

## Setup in MORIS

**Mesh (HMR).** 4 × 2 × 2 elements. One Lagrange mesh of order 2 (pattern 0) and two B-spline
meshes: order 2 on pattern 0 and order 1 on pattern 1 (`bspline_orders = "2,1"`,
`lagrange_to_bspline = "0,1"`). Refinement and staircase buffer 0.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Plane`, returning x − 1.2.

**Cut and enrichment (XTK).** Conformal decomposition with enrichment of both B-spline meshes
(`enrich_mesh_indices = "0,1"`); `ghost_stab = true`, although the deck adds no ghost IWGs
(`tUseGhost = false`)
([coupling C6](../../../user-guide/interlinks.md#c6-mesh-indices-hmr-xtk-gen)).

**DOF discretization (MSI).** `UX`, `UY`, `UZ` → B-spline mesh 1 (linear); `TEMP` → B-spline
mesh 0 (quadratic).

**Physics (FEM).** Per phase: a `STRUC_LIN_ISO` model (`FULL`) with `STRUC_LINEAR_BULK` and the
body load, and a `DIFF_LIN_ISO` model with `SPATIALDIFF_BULK` and the heat source. On the x = 0
face, `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` and `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE`
(`DIRICHLET_NITSCHE` 100.0 each). On the interface, `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE`
and `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE` 10.0 each). IQIs: `DOF` for
UX, UY, UZ and TEMP over both phases.

**Solver (SOL).** `AMESOS_IMPL` with `Solver_Type = "Amesos_Umfpack"`, combined residual/Jacobian
assembly, one monolithic nonlinear solver on `UX,UY,UZ,TEMP`, one time step
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Output (VIS).** `Bar_Lin_Disp_Quad_Diff.exo` on both phases with nodal UX, UY, UZ, TEMP.

## What the test checks

The test calls the workflow interface without checking its return code, then (serial only)
opens `Bar_Lin_Disp_Quad_Diff.exo` and requires:

| Quantity | Reference | Tolerance |
|---|---|---|
| dimensions / nodes / elements | 3 / 390 / 332 | exact |
| node 232 coordinates | (2.0, 0.125, 0.125) | relative 10⁻¹² |
| time value | 1.0 | relative 10⁻¹² |
| UX | 2.003084205519877 | relative 10⁻⁶ |
| UY | −4.907031775787233e−13 | absolute 10⁻⁸ |
| UZ | 1.806606057303358e−12 | absolute 10⁻⁸ |
| TEMP | 3.000000000000023 | relative 10⁻⁶ |

A `FIXME` in the test notes that the tolerances are loose because serial and parallel results
differ; only the serial case is configured.

## Run it

```bash
ctest -R Bar_Lin_Disp_Quad_Diff-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The test sets `gInterpolationOrder = 1`; the deck does not read it (orders are hard-coded in
  the HMR list).
- `PropDispNeumann`, `PropTempNeumann` (100.0) and a `SPATIALDIFF_NEUMANN` IWG appear only as
  unused properties or commented-out code.
- The `.dox` page for this directory has no content beyond its title.
