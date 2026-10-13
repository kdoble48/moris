# Hertzian contact benchmark (nonconformal, small strain)

A 2D frictionless Hertzian-contact setup: an elastic indenter with a parabolic lower boundary
(local radius R = 0.5) is pressed under displacement control into a flat elastic body. Both
bodies and the void gap between them are immersed in one background mesh, and contact is imposed
on the nonconformal interface with the neutral unbiased Nitsche contact IWG. The deck sweeps
three penetration depths against three uniform mesh sizes and writes the full 2D stress tensor of
each body on a discontinuous output mesh, so that the interface normal traction can be compared
with the Hertz pressure distribution in post-processing.

```{exa-summary} projects/EXA/structure/contact/Struc_Hertz_Contact
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square, N × N elements, linear interpolation |
| Indenter (phase 1) | above the parabola y = 0.5035 + (x − 0.5)², R = 1/(2 · 1.0) = 0.5 |
| Flat body (phase 2) | below the line y = 0.5025 |
| Gap (phase 0) | void band between the two; initial gap 0.001 at x = 0.5 |
| Material | both bodies linear elastic, E = 1000.0, ν = 0.0 |
| Supports | bottom face fixed; top face displaced by (0, −(0.001 + δ)), scaled by the Newton load factor |
| Penetration δ | 0.0025, 0.005 or 0.01 |
| Mesh N | 40, 80 or 160 |
| Analysis | quasi-static, small strain, load-controlled Newton over one time step (time frame 10.0) |

The deck header describes the intended post-processing: recover the interface normal traction
σ_nn = n_x² σ_xx + 2 n_x n_y σ_xy + n_y² σ_yy with the analytic parabola normal
n = (2(x − 0.5), −1)/√(1 + 4(x − 0.5)²), and compare it to the Hertz distribution
p(x) = p₀ √(1 − (x/a)²). Its comments estimate a/R ≈ 0.10, 0.14 and 0.20 for the three
penetrations from the geometric relation a ≈ √(2Rδ). The test does not perform this comparison.

## Setup in MORIS

The deck is adapted from {doc}`Parabolic_Indenter <Parabolic_Indenter>` (geometry) and
{doc}`Struc_Contact_Patch_Test <Struc_Contact_Patch_Test>` (phase-based contact wiring, one
workflow per test case, bulk-only discontinuous output). Mesh, phases, supports, ghost
stabilization and nonconformal settings match those decks.

**Mesh (HMR).** N × N elements on [0, 1]² with N = 40, 80, 160 selected by `gCaseIndex % 3`; one
Lagrange and one B-spline mesh of order `gInterpolationOrder`; no refinement.

**Geometry and phases (GEN).** A `USER_DEFINED` level set `Parabola` with (x₀, y₀, c) =
(0.5, 0.5035, 1.0) and a `LINE` through (0.5, 0.5025) with normal (0, 1); `Phase_Index_Split`
gives phases 1 (indenter), 0 (gap) and 2 (flat body).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on.

**Physics (FEM).** Phases `PhaseVoid`, `PhaseTop`, `PhaseBottom`, `PhaseAll`. Per body:
`STRUC_LIN_ISO`, `STRUC_LINEAR_BULK`, a symmetric Nitsche support (`DIRICHLET_NITSCHE`, 100.0) on
side ordinal 3 (top) or 1 (bottom), and `GHOST_NORMAL_FIELD` with `GHOST_DISPL` 0.05. Contact uses
two `STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED` IWGs of element type
`NONCONFORMAL_SIDESET`, one per leader body, with a `DIRICHLET_NITSCHE` parameter
`"100.0;1.0"` (γ = `tContactGamma` = 100) in the `NitscheInterface` slot and a finite-difference
Jacobian (3-point central, absolute perturbation 10⁻⁸). Nonconformal settings: `BAR_6`, maximum
ray lengths −0.005 and 0.01. IQIs: `DOF` (UX, UY); per body `NORMAL_STRESS` with
`vectorial_field_index` 0 (σ_xx) and 1 (σ_yy) and `SHEAR_STRESS` with index 2 (σ_xy);
`CONTACT_PRESSURE_REFERENCE` per body (defined, not written); and `GAP` with
`IQI_bulk_type = NONCONFORMAL_SIDESET`, leader `PhaseBottom`, follower `PhaseTop`.

**Solver (SOL).** `Amesos_Umfpack`, no preconditioner. Newton with combined residual/Jacobian
assembly, at most 200 iterations, `NLA_rel_res_norm_drop = 1e-8`, relaxation 1.0, linear load
control (factor 0.1, 10 steps, `NLA_load_control_relres = 1e-2`), remapping with
`MixedNthLoadStepAndResidualChange` (every load step, every 10 iterations, tolerance −10⁻²).
`TSA_Output_Indices = "0,1"` triggers both output meshes.

**Output (VIS).** Two output meshes with distinct `Output_Index`
([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)):

| Index | File | Mesh type | Sets | Fields |
|---|---|---|---|---|
| 0 | `Struc_Hertz_Contact_Case_<i>` | `FULL_DISCONTINUOUS` | bulk sets of both bodies | nodal `UX, UY`; `STRESSBulk<Side><Comp>` (`ELEMENTAL_AVG`) and `STRESSNodal<Side><Comp>` (`NODAL`) for Comp = XX, YY, XY |
| 1 | `Struc_Hertz_Contact_Case_<i>_gap` | `STANDARD` | bulk sets and the two nonconformal contact sets | `GAP` (`FACETED_AVG`) |

The deck comments explain the split: a nonconformal side set in a `FULL_DISCONTINUOUS` mesh
crashes `VIS_Factory::populate_nodal_point_pairs`, so the gap is written on a separate mesh.

## Cases

`gCaseIndex = load · 3 + mesh`:

| δ \ N | 40 | 80 | 160 |
|---|---|---|---|
| 0.0025 | 0 | 1 | 2 |
| 0.005 | 3 | 4 | 5 |
| 0.01 | 6 | 7 | 8 |

Each case is its own Catch2 test case with tag `[EXA_Struc_Hertz_Contact_<i>]` and its own ctest,
because, per the test header, ArborX/Kokkos cannot run a second workflow in one process.

## What the test checks

Each case sets `gInterpolationOrder = 1` and `gCaseIndex`, runs the workflow once, and opens
`Struc_Hertz_Contact_Case_<i>.e-s.0000`. The return code of the workflow call is not checked.
Without `MORIS_HAVE_ARBORX` the case only logs that it is skipped; with it, the case aborts unless
run on one processor. The checks are structural and finiteness-only:

- 2 dimensions, more than 0 nodes and elements;
- σ_yy sampled from `STRESSNodalTopYY` at nodes with y > 0.563 and from `STRESSNodalBottomYY` at
  nodes with y < 0.443 (excluding ±0.06 around y = 0.503); non-finite values on duplicated
  vertices are skipped and counted;
- at least one sample per body, and finite minimum and maximum of σ_yy in each body.

Per-body mean, minimum and maximum σ_yy and a `HERTZ_SWEEP` summary line are logged but not
compared to references; a strict nodal check is disabled (`gHaveReference = false`).

## Run it

```bash
ctest -R Struc_Hertz_Contact_ --output-on-failure            # all nine cases
ctest -R Struc_Hertz_Contact_8-1-procs --output-on-failure   # one case
```

in a build with `MORIS_USE_EXAMPLES=ON`. The ctest names are `Struc_Hertz_Contact_<i>-1-procs`
for i = 0…8. The `CMakeLists.txt` comment records case 8 (δ = 0.01, 160 × 160) at about 78 s.
