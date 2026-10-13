# Stefan's problem: diffusion with phase change (conformal mesh)

A transient 1D melting problem (Stefan's problem) on a thin 2D strip: a bar initially at 313
is heated from one end to 350, and latent heat is absorbed over a melting range around 314.5.
The example is adapted from a COMSOL benchmark and is one half of a pair. This deck meshes the
bar conformally with HMR; its twin {doc}`Stefans_Problem_cut <Stefans_Problem_cut>` immerses the
same bar in a larger background mesh and cuts it with XTK. Both exercise the phase-change
diffusion model `DIFF_LIN_ISO_PC` with GGLS stabilization and time-continuity (time-slab)
integration.

```{exa-summary} projects/EXA/thermal/diffusion/Stefans_Problem_conform
```

## Problem

| | |
|---|---|
| Domain | 0.028 × 0.0005 strip, 50 × 1 elements (h = 5.6 × 10⁻⁴) |
| Density, heat capacity | ρ = 0.75, c_p = 2.4 |
| Conductivity | κ = 2.1 × 10⁻⁷ |
| Latent heat | 175 |
| Phase-change temperature | 314.5 |
| Phase-change constant (`PhaseChangeConst`) | 3.0 (melting range ΔT_m = 3.0 in the `.dox`) |
| Phase-state function (`PhaseStateFunction`) | 2.0 |
| Initial temperature | 313.0 |
| Wall temperature | `SideSet_4` (x = 0), ramped linearly from 313 to 350 over the first 6 time slabs, then held at 350 |
| Other boundaries | adiabatic (no condition) |
| Time | 10 time slabs over t_max = 480 |

The `.dox` file gives the semi-analytical solution, T(x, t) = T_w + (T_m − T_w) erf(η)/erf(β)
with η = x / (2√(αt)), α = κ/(ρc_p), and β from β e^{β²} erf(β) = Ste/2 with
Ste = c_p (T_w − T_m)/L_H. `Stefans_Problem.m` and `solve_stefans_problem_1D.m` evaluate it; the
test does not compare against it.

## Setup in MORIS

**Mesh (HMR).** 50 × 1 elements on 0.028 × 0.0005. Lagrange and B-spline orders follow
`gInterpolationOrder` (1 or 2; any other value aborts with an error). Refinement and staircase
buffer 3, `adaptive_refinement_level` = 1.

**Geometry and phases (GEN).** One `USER_DEFINED` level set `Dummy_LS` = x + 1, positive on the
whole domain, so nothing is cut and every element belongs to phase 1.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, Exodus output of the integration mesh.

**Physics (FEM).** Legacy set-name input on `HMR_dummy_n_p1`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
The constitutive model `DIFF_LIN_ISO_PC` takes conductivity, density, heat capacity, latent heat,
phase-change temperature, phase-state function and phase-change constant. `SPATIALDIFF_BULK`
uses the stabilization parameter `GGLS_DIFFUSION` (same seven properties).
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, penalty 1000.0) imposes the
wall temperature on `SideSet_4_n_p1`. `TIME_CONTINUITY_DOF` (weights 100.0 current and previous,
initial condition 313.0) couples the time slabs. A `GHOST_DISPL` parameter (0.01) is defined, but
no ghost IWG uses it. IQIs: `DOF` temperature, `MAX_DOF` (function parameters `313.0/30.0`) and
`LATENT_HEAT_ABSORPTION`.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type. Newton with relative residual drop
10⁻⁷, at most 40 iterations, relaxation strategy `InvResNormAdaptive` (parameter 0.5, damping
0.5). Time solver: 10 steps over 480, two time levels for `TEMP`.

**Output (VIS).** `Comsol_conform_<gTestCaseIndex>.exo` on `HMR_dummy_n_p1`: nodal `TEMP`,
global `MAX_DOF` and `LATENT_HEAT_ABSORPTION`.

## Conformal vs cut

| | Stefans_Problem_conform | {doc}`Stefans_Problem_cut <Stefans_Problem_cut>` |
|---|---|---|
| Background mesh | 50 × 1 on 0.028 × 0.0005, no offset | 26 × 3 on 0.030 × 0.0030, offset (−0.0015, −0.0008) |
| Bar geometry | the background mesh itself (dummy level set) | four planar level sets bounding 0.028 × 0.0015 |
| Bulk / BC sets | `HMR_dummy_n_p1`, `SideSet_4_n_p1` | `HMR_dummy_{n,c}_p15`, `iside_b0_15_b1_14` |
| Adaptive refinement level | 1 | not set |
| Time slabs / ramp length | 10 / 6 slabs (t = 288) | 15 / 9 slabs (t = 288) |
| Newton | drop 10⁻⁷, 40 iterations, adaptive relaxation | drop 2 × 10⁻⁵, 20 iterations, relaxation 0.96 |
| IQIs in output | temperature, max temperature, latent heat | temperature |

## What the test checks

Two Catch2 cases each run the deck, require a return code of 0, and abort unless run on one
processor. Each checks the output file named after `gTestCaseIndex`:

| Case | Order | File | Nodes / elements | Reference node | Temperature (step index 9) | Tolerance |
|---|---|---|---|---|---|---|
| `Comsol_conform_Linear` | 1 | `Comsol_conform_0.exo` | 102 / 50 | 14 | 323.3263141372428 | 10⁻⁵ relative |
| `Comsol_conform_Quadratic` | 2 | `Comsol_conform_1.exo` | 303 / 50 | 39 | 323.3688441192814 | 10⁻³ relative |

Both cases also check 2 dimensions, the reference node coordinates (3.92 × 10⁻³, 0.0) to 10⁻⁸,
and a time value of 48.0 to 10⁻⁸. A comment in the source states that the quadratic case fails
on some machines, with results that depend on the BLAS implementation.

## Run it

```bash
ctest -R Comsol_conform --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is exactly `Comsol_conform` (`NO_PROCS_SUFFIX`).

## Notes

- The `.dox` file states a domain length of 0.28 and a 500 × 1 mesh; the deck uses 0.028 and
  50 × 1 (the element size 5.6 × 10⁻⁴ agrees). A commented-out line `tXlength = 0.28` remains in
  the deck.
- The `.dox` file attributes the Dirichlet ramp to instabilities when the full wall temperature
  is imposed from the first step, and notes that with this few time steps the results have
  limited meaning; the COMSOL benchmark uses a time frame of 57600.
