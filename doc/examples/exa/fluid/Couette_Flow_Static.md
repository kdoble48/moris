# Couette flow static

Steady incompressible flow in the annulus between two concentric circles (Taylor–Couette
geometry), immersed in a square background mesh. The inner and outer circles carry the
tangential velocity of the analytic Couette solution, imposed with Nitsche's method, and the
deck integrates the L2 and H1 errors of the computed velocity against that solution. It
demonstrates immersed incompressible flow on curved boundaries with ghost stabilization and
the analytic-error IQIs.

```{exa-summary} projects/EXA/fluid/laminar/Couette_Flow_Static
```

## Problem

| | |
|---|---|
| Domain | annulus 1.0 ≤ r ≤ 2.0 about the origin |
| Background domain | 5 × 5 square with offset (−2.25, −2.25) |
| Phases | phase 2 (between the circles) is the fluid |
| Material | density 1.0, viscosity 1.0 |
| Boundary conditions | v_θ(r) = A r + B/r prescribed on both circles, Nitsche-enforced |
| Rotation parameters | ω_in = −5, ω_out = 5, Kn = 0 |
| Time | steady; Newton iteration from VX = VY = 10⁻⁴, P = 0 |

The constants follow from the deck's formulas with Kn = 0: b = 1/R_in = 1, d = 1/R_out = 0.5,
B = (ω_out R_out − ω_in R_in R_out/R_in)/(d − b R_out/R_in) = −40/3 ≈ −13.33 and
A = (ω_in R_in − b B)/R_in = 25/3 ≈ 8.33, which give v_θ = −5 at r = 1 and v_θ = 10 at r = 2.
`ImposedVelocityFunc` evaluates (v_x, v_y) = v_θ (−y/r, x/r); `AnalyticdVelocitydxFunc` returns
the analytic velocity gradient.

## Setup in MORIS

**Mesh (HMR).** 22 × 22 elements on the 5 × 5 domain, linear Lagrange and B-spline meshes
(`lagrange_to_bspline = "0"`), refinement and staircase buffers of 3,
`adaptive_refinement_level` 2.

**Geometry and phases (GEN).** Two `CIRCLE` level sets with radii 1.0 and 2.0 (only the
radius is set; the centre keeps its default (0, 0)). The default phase table over two geometries is used; the fluid is phase 2,
the inner interface is `iside_b0_2_b1_0` and the outer interface `iside_b0_2_b1_3`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, `ghost_stab` on
(set `ghost_p2`); the enriched integration mesh is printed and written to Exodus.

**Physics (FEM).** The block uses `hack_for_legacy_fem()` and explicit `mesh_set_names`.
Constitutive model `FLUID_INCOMPRESSIBLE`. Stabilization parameters: `INCOMPRESSIBLE_FLOW`
(36.0), `VELOCITY_DIRICHLET_NITSCHE` ("100.0/1.0"), `VISCOUS_GHOST` (0.05), `CONVECTIVE_GHOST`
(0.05), `PRESSURE_GHOST` ("0.005/1.0"). IWGs:

- `INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK` on
  `HMR_dummy_n_p2, HMR_dummy_c_p2`;
- `INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_SYMMETRIC_NITSCHE` and
  `INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE` on each circle, with the
  properties `PropDirichletInU` and `PropDirichletOutU` (both `ImposedVelocityFunc`);
- three `GHOST_NORMAL_FIELD` IWGs (viscous, convective, pressure) on `ghost_p2`.

IQIs: `DOF` for VX, VY, P; `L2_ERROR_ANALYTIC` on VX, VY against `PropDirichletInU`
(`L2Check`); `H1_ERROR_ANALYTIC` against `PropdUdx` (`H1Check`).

**Solver (SOL).** `AMESOS_IMPL` direct solve; Newton with `NLA_rel_res_norm_drop` 10⁻⁶,
relaxation 1.0, at most 20 iterations, on `VX,VY;P`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). No
preconditioner.

**Output (VIS).** `Couette_Flow_Static.exo` on the fluid sets and both interface sets, with
nodal `VX, VY, P, L2_ERROR_ANALYTIC, H1_ERROR_ANALYTIC`.

## What the test checks

The test case sets `gPrintReferenceValues = false`, calls the workflow with
`./Couette_Flow_Static.so`, and makes no assertion: neither the return code nor any output
value is checked, although the file includes the Exodus reader headers. The example is
registered for 1 and 2 processors.

## Run it

```bash
ctest -R Couette_Flow_Static --output-on-failure   # runs Couette_Flow_Static-1-procs and -2-procs
```

## Notes

- The `.dox` page (`EXA_Couette_Flow_Static`) contains section headings only.
- The L2 and H1 errors are written to the Exodus file but are not compared by the test.
