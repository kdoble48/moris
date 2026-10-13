# Parabolic indenter (nonconformal contact)

A 2D contact problem: an elastic body with a parabolic lower boundary is pressed into
a second elastic body under displacement control. Both bodies and the thin void gap between them
are immersed in one background mesh, and the contact condition is imposed on the nonconformal
interface between the two cut boundaries by an unbiased Nitsche contact IWG, with the gap
evaluated by ray tracing. The directory holds a small-strain deck, which the test runs, and a
finite-strain deck with the Mlika contact formulation.

```{exa-summary} projects/EXA/structure/contact/Parabolic_Indenter
```

## Problem

Values for `Parabolic_Indenter_Linear.cpp`, the deck the test runs:

| | |
|---|---|
| Domain | 1.0 × 1.0 square, 10 × 10 elements, linear interpolation (`gInterpolationOrder = 1`) |
| Top body (phase 1) | above the parabola y = 0.5035 + (x − 0.5)² |
| Bottom body (phase 2) | below the line y = 0.5025 |
| Gap (phase 0) | void band between the two; initial gap 0.001 at x = 0.5 |
| Material | both bodies linear elastic, E = 1000.0, ν = 0.0 |
| Supports | bottom face of the bottom body fixed; top face of the top body displaced by (0, −0.051) (gap 0.001 + penetration 0.05), scaled by the Newton load factor |
| Analysis | quasi-static, small strain, load-controlled Newton over one time step (time frame 10.0) |

## Setup in MORIS

**Mesh (HMR).** 10 × 10 elements on [0, 1]², one Lagrange and one B-spline mesh of order
`gInterpolationOrder`, `refinement_buffer = 0`, `staircase_buffer = 0`.

**Geometry and phases (GEN).** Two level sets: a `USER_DEFINED` field `Parabola`,
φ = (y − y₀) − c (x − x₀)² with (x₀, y₀, c) = (0.5, 0.5035, 1), and a `LINE` through (0.5, 0.5025)
with normal (0, 1). The phase function `Phase_Index_Split` returns 1 above both, 0 above the line
and below the parabola, and 2 below the line (`number_of_phases = 3`).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`, Exodus output of the integration mesh. The contact sets
are nonconformal side-set pairs such as `ncss|iside_b0_1_b1_0|iside_b0_2_b1_0`, formed from the
phase-1/void and phase-2/void interfaces
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Physics (FEM).** The FEM list is phase-based: `PhaseVoid` (0), `PhaseTop` (1), `PhaseBottom` (2),
`PhaseAll` (1, 2). Each body has a `STRUC_LIN_ISO` model, a `STRUC_LINEAR_BULK` IWG, and a
`GHOST_NORMAL_FIELD` IWG with a `GHOST_DISPL` parameter of 0.05. The supports are
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` IWGs on side ordinal 3 of the top body and side ordinal 1
of the bottom body (`DIRICHLET_NITSCHE`, 100.0); the top value function `Func_Dirichlet` multiplies
the prescribed displacement by the `LoadFactor` that the Newton solver logs. Contact uses two
`STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED` IWGs of element type
`NONCONFORMAL_SIDESET` — one with the top body as leader, one with the bottom body as leader, each
with `neighbor_phases = PhaseVoid` — and a `DIRICHLET_NITSCHE` parameter with
`function_parameters = "10.0;1.0"` in the `NitscheInterface` slot. The contact Jacobian is computed
by finite differences (`analytical_jacobian = false`; 3-point central scheme, absolute
perturbation 10⁻⁸). Nonconformal settings: integration order `BAR_6`, maximum ray lengths −0.005
and 0.01. IQIs: `DOF` (UX, UY), `NORMAL_STRESS` (x, y) per body, `CONTACT_PRESSURE_REFERENCE` per
body, and `GAP`.

**Solver (SOL).** Direct solve with `AMESOS_IMPL`, `Solver_Type = Amesos_Umfpack`, no
preconditioner. Newton with combined residual/Jacobian assembly, at most 200 iterations,
`NLA_rel_res_norm_drop = 1e-8`, relaxation 1.0, and linear load control (factor 0.1, 15 steps,
`NLA_load_control_relres = 1e-2`). The nonconformal mapping is recomputed with the
`MixedNthLoadStepAndResidualChange` strategy (every load step, every 10 iterations, residual-change
tolerance −10⁻²). One time step over a time frame of 10.0.

**Output (VIS).** `STANDARD` mesh on both bodies' bulk sets and the contact interface sets:
`GAP` (`FACETED_AVG`), nodal `UX, UY`, `STRESSBulk<Side><X|Y>` (`ELEMENTAL_AVG`),
`ContactPressure<Side>` (`FACETED_INT`) and `NORMALX, NORMALY` (`FACETED_AVG`), with
`Time_Offset = 0.1`, written as `Parabolic_Indenter_Linear_Case_<gCaseIndex>`.

## Decks

| | `Parabolic_Indenter_Linear.cpp` | `Parabolic_Indenter_NonLinear.cpp` |
|---|---|---|
| Domain | [0, 1]², 10 × 10 | [−0.5, 0.5]², 11 × 11, `truncate_bsplines = true`, buffers 2 |
| Bodies | parabola (c = 1) on a flat line | two parabolas, y = 0.005 + 0.5x² and y = −0.005 − 0.5x² |
| Gap, penetration | 0.001, 0.05 | 0.01, 0.02 |
| Loading | top face moves −0.051, bottom fixed | top and bottom faces each move 0.015 towards each other |
| Lateral faces | free | symmetry (x displacement fixed) on both side faces of both bodies |
| Material | `STRUC_LIN_ISO`, ν = 0.0 | `STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF`, ν = 0.3 |
| Bulk / support IWGs | `STRUC_LINEAR_BULK`, `STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` | `STRUC_NON_LINEAR_BULK_PF`, `STRUC_NON_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE_PF` |
| Contact IWG | `STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED`, γ = 10 | `STRUC_NONLINEAR_CONTACT_MLIKA_UNBIASED_UNSYMMETRIC`, γ = 30 |
| Newton | relaxation 1.0, load control 0.1 × 15 steps | relaxation 0.6, load control 0.01 × 20 steps, relres 1.0 |
| Remapping | mixed strategy, rays −0.005 / 0.01, `BAR_6` | every iteration, rays −0.08 / 0.5, `BAR_2` |
| Linear solver | `Amesos_Umfpack` | Amesos default (`Solver_Type` line commented out) |
| Output fields | GAP, UX/UY, stresses, contact pressure, normals | per-body displacements, `NORMAL_STRESS_CAUCHY`, `VON_MISES_STRESS_CAUCHY`, per-body GAP |

## What the test checks

The single test case `Parabolic_Indenter_Linear` runs only the small-strain deck. It is active
only when MORIS is built with ArborX (`MORIS_HAVE_ARBORX`); otherwise it logs that the example is
skipped. It aborts unless run on one processor. The return code of the workflow call is not
checked. The test opens `Parabolic_Indenter_Linear_Case_0.e-s.0000` and requires:

- 2 dimensions, 313 nodes, 192 elements;
- node 294 at (0.5, 0.5025), relative tolerance 10⁻⁵;
- time value 10.0, relative tolerance 10⁻⁸;
- displacement at node 294 equal to (−3.534434707883998e−06, −2.785067809916309e−02), relative
  tolerance 10⁻⁵.

The nonlinear deck is compiled as an input but no test case runs it.

## Run it

```bash
ctest -R Parabolic_Indentor-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name keeps the legacy spelling `Indentor`. A build without ArborX registers the test
but skips the run.

## Notes

- The small-strain VIS block requests `IQINormalX` and `IQINormalY`, which its FEM block does not
  define.
- The same two-body setup, with a flat top body, is the basis of
  {doc}`Struc_Contact_Patch_Test <Struc_Contact_Patch_Test>`; the parabolic geometry is reused in
  {doc}`Struc_Hertz_Contact <Struc_Hertz_Contact>`.
