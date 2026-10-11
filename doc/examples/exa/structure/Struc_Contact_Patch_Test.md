# Contact patch test (nonconformal, small strain)

A constant-stress contact patch test: two flat linear-elastic blocks, separated by a thin void
band and immersed in one background mesh, are compressed vertically until the gap closes. Per the
deck header, with ν = 0 the exact solution is a spatially uniform uniaxial normal stress σ_yy
across the nonconformal contact interface. The deck sweeps four contact formulations — symmetric,
unsymmetric and neutral unbiased Nitsche, and penalty — against four values of the contact
stabilization γ, one case per process.

```{exa-summary} projects/EXA/structure/contact/Struc_Contact_Patch_Test
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square, 10 × 10 elements, linear interpolation |
| Top block (phase 1) | above the line y = 0.5035 |
| Bottom block (phase 2) | below the line y = 0.5025 |
| Gap (phase 0) | void band of width 0.001 between the two lines |
| Material | both blocks linear elastic, E = 1000.0, ν = 0.0 |
| Supports | bottom face fixed; top face displaced by (0, −0.006) (gap 0.001 + penetration 0.005), scaled by the Newton load factor |
| Analysis | quasi-static, small strain, load-controlled Newton over one time step (time frame 10.0) |

## Setup in MORIS

The deck is the flat-interface variant of
{doc}`Parabolic_Indenter <Parabolic_Indenter>` (`Parabolic_Indenter_Linear.cpp`): same mesh,
phase function, phase-based FEM list, Dirichlet supports, ghost stabilization and nonconformal
settings. The differences are listed here.

**Geometry and phases (GEN).** Both level sets are `LINE` geometries with normal (0, 1), through
(0.5, 0.5035) and (0.5, 0.5025). `Phase_Index_Split` maps them to phases 1 (top), 0 (gap) and
2 (bottom); the `Parabola` field function is kept in the deck but not used.

**Physics (FEM).** The global `gCaseIndex` selects the contact IWG and γ as
`method = gCaseIndex / 4`, `gamma = gCaseIndex % 4`:

| method | Contact IWG | Stabilization parameter (slot) | Jacobian |
|---|---|---|---|
| 0 | `STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED` | `DIRICHLET_NITSCHE`, `"γ;1.0"` (`NitscheInterface`) | finite difference |
| 1 | `STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED` | same | finite difference |
| 2 | `STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED` | same | finite difference |
| 3 | `STRUC_LINEAR_CONTACT_PENALTY` | `PENALTY_CONTACT`, `"γ"`, leader and follower Young's modulus (`PenaltyContact`) | analytic |

with γ ∈ {10, 100, 1000, 10000} for `gamma` = 0, 1, 2, 3. As in the indenter, two contact IWGs of
element type `NONCONFORMAL_SIDESET` are defined, one per leader body, with
`neighbor_phases = PhaseVoid`. Each body has a `STRUC_LIN_ISO` model, a `STRUC_LINEAR_BULK` IWG,
symmetric Nitsche supports (`DIRICHLET_NITSCHE`, 100.0) and a `GHOST_NORMAL_FIELD` IWG
(`GHOST_DISPL`, 0.05). IQIs: `DOF` (UX, UY), `NORMAL_STRESS` (x, y) per body,
`CONTACT_PRESSURE_REFERENCE` per body, and `GAP`.

**Solver (SOL).** `Amesos_Umfpack`, no preconditioner. Newton with at most 200 iterations,
`NLA_rel_res_norm_drop = 1e-8`, relaxation 1.0, linear load control (factor 0.2, 5 steps,
`NLA_load_control_relres = 1e-2`), remapping with `MixedNthLoadStepAndResidualChange`.
`NLA_combined_res_jac_assembly` is `true` for the Nitsche methods and `false` for the penalty
method; the deck comment states that the penalty IWG does not implement the combined
residual/Jacobian evaluation.

**Output (VIS).** A `FULL_DISCONTINUOUS` mesh on the bulk sets of both blocks only
(`HMR_dummy_{n,c}_p1`, `HMR_dummy_{n,c}_p2`), so stresses are not averaged across the interface.
Fields: nodal `UX, UY`; `STRESSBulk<Side><X|Y>` (`ELEMENTAL_AVG`); and `STRESSNodal<Side>Y`, the
same σ_yy IQIs projected to nodes. `GAP` and the contact pressure are defined but not written; the
deck comments state that a nonconformal side set in `Set_Names` of a `FULL_DISCONTINUOUS` mesh
crashes `VIS_Factory::populate_nodal_point_pairs`. Output: `Struc_Contact_Patch_Test_Case_<i>`.

## Cases

| `gCaseIndex` | method | γ | Registered as ctest |
|---|---|---|---|
| 0, 1, 2, 3 | symmetric Nitsche | 10, 100, 1000, 10000 | 0, 1, 2 |
| 4, 5, 6, 7 | unsymmetric Nitsche | 10, 100, 1000, 10000 | 4, 5, 6 |
| 8, 9, 10, 11 | neutral Nitsche | 10, 100, 1000, 10000 | 8, 9, 10 |
| 12, 13, 14, 15 | penalty | 10, 100, 1000, 10000 | none |

All 16 test cases are compiled and carry the tag `[EXA_Struc_Contact_Patch_Test_<i>]`. The
`CMakeLists.txt` comment gives the reasons for the unregistered cases: for γ = 10⁴ (cases 3, 7, 11)
Newton diverges (residual above 10²⁰), and the penalty cases 12–15 end in a SIGSEGV in
`fem::Property::val()` via `SP_Penalty_Contact` on the nonconformal contact set.

## What the test checks

Each test case sets `gInterpolationOrder = 1` and `gCaseIndex`, runs the workflow once (one
workflow per process, because the test header states ArborX/Kokkos is not re-entrant), and then
checks `Struc_Contact_Patch_Test_Case_<i>.e-s.0000`. The return code of the workflow call is not
checked. Without `MORIS_HAVE_ARBORX` the case only logs that it is skipped; with it, the case
aborts unless run on one processor.

The checks are structural and finiteness-only:

- 2 dimensions, more than 0 nodes and elements;
- σ_yy is sampled from `STRESSNodalTopY` at nodes with y > 0.563 and from `STRESSNodalBottomY` at
  nodes with y < 0.443 (a band of ±0.06 around y = 0.503 excludes the cut rows); non-finite values
  on duplicated vertices of the discontinuous mesh are skipped and counted;
- at least one sample per block, and the per-block spread (max − min) and the cross-block
  difference of the means must be finite.

The spread, the means and the jump are logged (a `GAMMA_SWEEP` line per case) but not compared to
a reference. A nodal-displacement regression against node 294 exists in the code but is disabled
(`gHaveReference = false`; the reference values are placeholders).

## Run it

```bash
ctest -R Struc_Contact_Patch_Test_ --output-on-failure      # all registered cases
ctest -R Struc_Contact_Patch_Test_9-1-procs --output-on-failure     # one case
```

in a build with `MORIS_USE_EXAMPLES=ON`. The registered ctest names are
`Struc_Contact_Patch_Test_<i>-1-procs` for i ∈ {0, 1, 2, 4, 5, 6, 8, 9, 10}. The unregistered
cases can be run by invoking the EXA test binary with the tag `[EXA_Struc_Contact_Patch_Test_<i>]`.
