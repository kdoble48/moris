# Stefan's problem: diffusion with phase change (cut mesh)

A transient 1D melting problem (Stefan's problem): a bar initially at 313 is heated from one end
to 350, and latent heat is absorbed over a melting range around 314.5. This deck is the cut-mesh
half of a pair: the bar is a 0.028 × 0.0015 rectangle bounded by four planar level sets and
immersed in a larger HMR background mesh, so XTK cuts every boundary of the bar. Its twin
{doc}`Stefans_Problem_conform <Stefans_Problem_conform>` uses a body-fitted HMR mesh; the
material model, stabilization and boundary conditions are the same.

```{exa-summary} projects/EXA/thermal/diffusion/Stefans_Problem_cut
```

## Problem

| | |
|---|---|
| Background domain | 0.030 × 0.0030, offset (−0.0015, −0.0008), 26 × 3 elements |
| Bar | 0.028 × 0.0015 rectangle, centre (0.014, 0.00075) |
| Density, heat capacity | ρ = 0.75, c_p = 2.4 |
| Conductivity | κ = 2.1 × 10⁻⁷ |
| Latent heat | 175 |
| Phase-change temperature | 314.5 |
| Phase-change constant / phase-state function | 3.0 / 2.0 |
| Initial temperature | 313.0 |
| Wall temperature | on the interface `iside_b0_15_b1_14`, ramped linearly from 313 to 350 over the first 9 time slabs (t = 288), then held at 350 |
| Other boundaries | no condition |
| Time | 15 time slabs over t_max = 480 |

The semi-analytical solution and the MATLAB scripts are documented with the conformal twin; this
directory holds no `.dox` or `.m` files of its own.

## Setup in MORIS

**Mesh (HMR).** 26 × 3 elements on 0.030 × 0.0030 with offset (−0.0015, −0.0008); Lagrange and
B-spline orders follow `gInterpolationOrder` (1 or 2). Refinement and staircase buffer 3, no
adaptive refinement.

**Geometry and phases (GEN).** Four `USER_DEFINED` level sets in this order: `Bottom_Boundary`,
`Top_Boundary`, `Right_Boundary`, `Left_Boundary`, each a plane that is positive inside the bar
[0, 0.028] × [0, 0.0015]. No phase function is set; FEM addresses the bar as phase 15.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets`, Exodus output of the integration mesh.

**Physics (FEM).** Legacy set-name input
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
`DIFF_LIN_ISO_PC` with the seven phase-change properties; `SPATIALDIFF_BULK` with
`GGLS_DIFFUSION` and `TIME_CONTINUITY_DOF` (weights 100.0, initial condition 313.0) on
`HMR_dummy_n_p15,HMR_dummy_c_p15`; `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE`
(`DIRICHLET_NITSCHE`, penalty 1000.0) on `iside_b0_15_b1_14`, the interface between the bar and
phase 14. As in the conformal deck, a `GHOST_DISPL` parameter (0.01) is defined but no ghost IWG
uses it. The only IQI is the `DOF` temperature.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type; Newton with relative residual drop
2 × 10⁻⁵, relaxation parameter 0.96, at most 20 iterations; 15 time steps over 480, two time
levels for `TEMP`.

**Output (VIS).** `Comsol_cut.exo` with nodal `TEMP` on the phase-15 bulk sets.

The differences from the conformal deck are tabulated on
{doc}`Stefans_Problem_conform <Stefans_Problem_conform>`.

## What the test checks

Two Catch2 cases each run the deck, require a return code of 0, abort unless run on one
processor, and check `Comsol_cut.exo` (written by both cases in turn):

| Case | Order | Nodes / elements | Reference node | Temperature (step index 14) | Tolerance |
|---|---|---|---|---|---|
| `Comsol_cut_Linear` | 1 | 386 / 235 | 26 | 328.9846289451254 | 10⁻⁴ relative |
| `Comsol_cut_Quadratic` | 2 | 479 / 235 | 53 | 328.7974558344104 | 10⁻⁴ relative |

Both cases also check 2 dimensions, the reference node coordinates
(3.115384615384614 × 10⁻³, 0.0012) to 10⁻⁸, and a time value of 32.0 to 10⁻⁸.

## Run it

```bash
ctest -R Comsol_cut --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is exactly `Comsol_cut` (`NO_PROCS_SUFFIX`).

## Notes

- The conformal twin's `.dox` describes the cut variant as using the same element size as the
  conformal mesh with no refinement; this deck's background element is 0.030/26 ≈ 1.15 × 10⁻³
  long in x, against 5.6 × 10⁻⁴ in the conformal deck.
