# MORIS bug ledger (kdoble48 fork)

One entry per MORIS defect found while using the code: crashes, wrong results, silent failures.
Missing features belong in a study's capability list, not here (exception: a missing piece that fails
*silently*, which is a bug). The ledger is versioned with the code, so a fix commit and its entry move together.

**Rules**
- New find → add an entry with the next `MB-NNN`, status `open`, before or while fixing it.
- Every fix commit message cites its ID (`fix(XTK): ... [MB-001]`) and the entry records the commit and the regression test.
- Status: `open` → `fix-written` (on a branch, not yet verified) → `fixed` (built, regression test + original
  reproducer pass) → `merged` (on `develop` of the fork). `wontfix`/`not-a-bug` with a reason.
- Reproducers point at durable paths (`runs/...` decks, unit tests). Workspace agents read this file before triage
  (`.claude/rules/moris-bug-ledger.md`).

## Index

| ID | Status | Module | Summary | Found by |
|---|---|---|---|---|
| MB-001 | fix-written | XTK | Segfault building the enriched integration mesh when a cell's vertices all lie on one geometry's interface → now a clear error | lunch_and_learn_xiga B13 |
| MB-002 | partly fixed | FEM/INT | Double-sided FD geometry sensitivity uses the leader's parametric coords for the follower → dR/dp ~1e8 on interfaces along cell edges | lunch_and_learn_xiga B12 |
| MB-003 | open | GEN/FEM | Static user-defined geometry with `discretization_mesh_index = -2` (analytic) makes design gradients NaN, silently | moment_fitting (2026-07-16) |
| MB-004 | open | FEM | Perimeter (interface VOLUME IQI) adjoint off by 12–150% on a few small-gradient ADVs | lunch_and_learn_xiga B11 |
| MB-005 | open | OPT/HMR | Adaptive refinement during optimization: ADV count drifts between remeshes → abort in `compute_design_criteria_gradients` | seeding_impact (2026-07-13) |
| MB-006 | open (suspect) | FEM/INT | Nonconformal branch of `select_dRdp_FD_geometry_double` also copies the leader ξ into the follower row | review of MB-002 |
| MB-007 | open | FEM/VIS | Requested IQI evaluated on zero sets writes zeros to Exodus with no warning | Hertz contact (2026-07) |
| MB-008 | open (limitation) | FEM/INT | `IQI_Max_Stress::compute_dQIdu` stops with "not implemented"; decks must switch the whole forward analysis to FD | moris_stress_demo (2026-10-06) |
| MB-009 | fix-written | GEN | Design factory rejects a level-set geometry listed before the fields it depends on ("does not recognize geometry type"); properties built on the wrong field | MB-001 regression test |
| MB-010 | open | XTK | Phase boundary along background grid lines creates no interface facets → phases silently decoupled | MB-001 verification |

---

## MB-001 · XTK segfault: cell with every vertex on an interface
- **Status:** fix-written on `fix/xiga-deck-bugs` (2026-10-09, second version); rebuild + test pending.
- **Symptom:** SIGSEGV "Address not mapped" in `xtk::Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()` at GCMMA iteration 1, right after "XTK - Enriched Integration Mesh - Construction". Opt build only.
- **Reproducer:** `studies/lunch_and_learn_xiga/decks/mbb_levelset.py`, MBB with `"elements": [60, 20]` and on-grid pads (no `pad_offgrid`), config `runs/studies/lunch_and_learn_xiga/configs/MB001_repro.json`. The two fixed pads are exactly 2×1 cells with every face on a grid line.
- **Root cause:**
  1. GEN `get_element_phase_index` returns `MORIS_INDEX_MAX` when every vertex of a cell has |φ| ≤ tolerance for one geometry.
  2. `XTK/src/cl_XTK_Integration_Mesh_Generator.cpp` `deduce_ig_cell_bulk_phase_from_facets` returned the sentinel `MORIS_INDEX_MAX - 1` for QUAD/HEX cells, which counted as assigned.
  3. It was then used unchecked as an index at `cl_XTK_Enriched_Integration_Mesh.cpp` (`tNoChildBlockSetOrds( tBulkPhaseIndex )`).
- **Fix, v1 (rejected 2026-10-09):** a GEN tie-break treating such cells as positive. The run then completed, but the pad was **silently decoupled**: design strain energy exactly 0, and no `dbl_iside_p0_1_p1_2` set processed. A phase boundary along grid lines gets no interface facets (MB-010), so v1 traded a crash for a wrong answer.
- **Fix, v2:** XTK stops with a `MORIS_ERROR` naming the cell and the remedy (shift the geometry off the grid lines). A `MORIS_ERROR` guard also sits before the block-set index. GEN is unchanged.
- **Regression test:** `XTK/test/xtk/UT_Snapping.cpp` "Band One Cell Thick on Grid Lines 2d" `[Snapping_2d]`, expecting `CHECK_THROWS( perform_decomposition() )`.
- **Workaround:** keep fixed-geometry faces off grid lines.

## MB-002 · FEM: follower parametric coordinates in double-sided FD geometry sensitivities
- **Status:** partly fixed on `fix/xiga-deck-bugs` (2026-10-09). After the fix the reproducer gives ADVs 262/342/343/349 = −6.58e-3/−6.64e-3/−2.64e-2/−6.93e-4 vs FD −6.93e-3/−7.00e-3/−2.68e-2/−6.21e-4: the 1e8 spikes are gone, signs and magnitudes are right, but 1.7–12% error remains on the grid-aligned interface row. The residual is under investigation (candidates: follower point leaving its element during the FD step, ghost facets on cell edges, the step size chosen from the leader only).
- **Symptom:** analytic d(strain energy)/dADV of order 1e7–1e8, sometimes with the wrong sign, for ADVs whose intersection nodes sit on a double-sided interface facet lying along a background-cell edge. Volume and perimeter gradients are fine, because they use explicit FD and don't go through λᵀ∂R/∂p. No NaN and no warning; GCMMA climbs and topology degrades.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/triage_bridge/sweep_pad20_idx1/`, a bridge with the deck pad's top face on node row J=2. SWEEP ADVs 262/342/343/349: analytic −1.25e8/+1.31e8/+6.9e7/−6.4e6 vs central FD −6.93e-3/−7.00e-3/−2.69e-2/−6.2e-4. The control with the pad face at 2.5h (`sweep_pad25_an`) is clean.
- **Root cause:** `FEM/INT/src/IWG/cl_FEM_IWG.cpp` `IWG::select_dRdp_FD_geometry_double` computed the perturbed node's ξ in the *leader* IP element and wrote it into the *follower* IG cell's parametric coefficients. That is only valid when both sides share one IP element. Across a cell edge (η = −1 vs +1) the follower point jumped a whole cell, giving an O(1) residual jump divided by δ ≈ 5.6e-9. GEN's dx/dφ was verified correct (the ADV ratio test matches, 1.90).
- **Fix:** for conformal sets, the follower ξ comes from the follower IP geometry interpolator (`update_parametric_coordinates`). The nonconformal path is unchanged (see MB-006).
- **Verification to run:** the reproducer above; pass = analytic within ~1e-4 relative of FD, and Criteria norms drop from 1e7–1e8 to O(1).
- **Regression test:** TODO, `FEM/INT/test/UT_FEM_IWG_Struc_Linear_Interface.cpp`, two-IP-element interface facet vs hand FD.
- **Was misdiagnosed first** as "fixed geometry on the grid" and as the −2 static path (MB-003). The −1 re-test ruled the latter out.

## MB-003 · GEN: analytic static user geometry poisons design gradients with NaN
- **Status:** open (deck workaround known).
- **Symptom:** once a static (ADV-free) user geometry's boundary becomes an exposed material interface, every `GCMMA - Criteria(N): norm = -nan`. The design freezes, with no error or warning.
- **Reproducer:** moment_fitting p2 padded warm restart (`runs/studies/moment_fitting/phase1_p2_opt/deck/boxbeam3d_frameless_lean_p2.py`, pads with index −2). A/B: −1 is clean.
- **Root cause (partial):** the analytic user-field path (`cl_GEN_User_Defined_Field`: FD `get_dfield_dcoordinates`, no `get_dfield_dadvs`) produces NaN in the intersection-node sensitivity chain.
- **Notes:** pymoris `add_user_geometry` defaults to −2, so every pymoris deck with fixed pads is exposed. Either fix GEN or make −2 error out for geometries without ADVs in optimization runs.
- **Workaround:** `discretization_mesh_index = -1`. Memory note: gen-static-geometry-stored-field.

## MB-004 · FEM: perimeter adjoint mismatch
- **Status:** open.
- **Symptom:** X6 SWEEP at the MBB seed. Strain-energy and volume adjoints match central FD (≤ 6.5e-5, ~1e-8). The perimeter criterion (`IQI VOLUME` on `iside_*` sets) is off by 12%, 18% and 150% on 3 of 12 ADVs. FD is stable to 4 digits over steps 1e-5 to 1e-7. The absolute size is ~1e-4 of the largest perimeter gradient.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X6_sweep/` (runB_fd, ADVs 417, 4015, 1103).
- **Root cause:** not investigated.

## MB-005 · OPT/HMR: adaptive refinement + optimization ADV drift
- **Status:** open.
- **Symptom:** abort with "incompatible matrix dimensions 1x4370 and 1x4374" in `opt::Problem::compute_design_criteria_gradients`. The B-spline ADV count changes between remeshes because the refinement pattern follows the evolving topology.
- **Reproducer:** seeding_impact deck family at 60×20 + refine-1 with remeshing (`studies/seeding_impact/FINDINGS.md`, S3-W2 it.70).
- **Workaround:** a static uniform mesh.

## MB-006 · FEM: nonconformal FD geometry sensitivity copies the leader ξ (suspect)
- **Status:** open (suspect, unverified).
- **Symptom:** none observed yet. Found by reading `select_dRdp_FD_geometry_double` while fixing MB-002. In the nonconformal branch the follower node is not perturbed, but its parametric row is overwritten with the leader's ξ before the ray remap.
- **Owner:** whoever owns the nonconformal contact work (`feature/nonconformal-*`). Check whether the remap later overrides it; if not, this is MB-002's twin.

## MB-007 · FEM/VIS: IQI on zero sets writes silent zeros
- **Status:** open (usability).
- **Symptom:** a contact IQI (`GAP`, `CONTACT_PRESSURE_*`) missing `IQI_bulk_type = NONCONFORMAL_SIDESET` is matched to 0 sets (`nIQIonSet=0`). VIS writes uninitialized zeros and nothing warns.
- **Suggested fix:** warn or error when a VIS-requested IQI is evaluated on no set. Memory note: nonconformal-contact-iqis-zero-at-vis.

## MB-008 · FEM: no analytic dQ/du for the p-norm stress IQI
- **Status:** open (limitation; it fails loudly, so it is kept here only because it forces whole-model FD).
- **Symptom:** `IQI_Max_Stress::compute_dQIdu` calls `MORIS_ERROR("not implemented")` (`cl_FEM_IQI_Max_Stress.hpp:99-114`). Stress-constrained decks must set `is_analytical_forward = false`, which is slow in 3D.

## MB-009 · GEN: design factory ignores unbuilt dependencies
- **Status:** fix-written on `fix/xiga-deck-bugs` (2026-10-09).
- **Symptom:** `GEN does not recognize geometry type: level_set` (`cl_GEN_Design_Factory.cpp`) when a geometry with
  `dependencies` comes before its fields, e.g. a COMBINED_FIELDS geometry in GEOMETRIES with its fields in the FIELDS
  sublist, which `Geometry_Engine` appends *after* the geometries.
- **Reproducer:** the MB-001 regression test as first written (fields in `tGENParams( 2 )`).
- **Root cause:** the design-building block ran even when `tCanBuild` was false. The level-set branch tested
  `level_set and tCanBuild`, so an unbuildable level set fell through to the error. The property branch had no check
  and would bind `mFields( tFieldIndex - 1 )` (the wrong field). Both erased `design_type` even when nothing was built.
- **Fix:** build designs only `if ( exists( "design_type" ) and tCanBuild )`; the existing outer loop retries.
- **Regression test:** `GEN/test/ut_GEN_Field.cpp` "Combined Field Declared Before Its Fields".

## MB-010 · XTK: phase boundary along background grid lines is not coupled
- **Status:** open (found 2026-10-09 while verifying MB-001 v1).
- **Symptom:** two phases meet along background-cell edges (e.g. a fixed pad whose faces lie exactly on grid lines,
  next to non-cut design cells). No double-sided interface facets are created there, so an interface IWG (Nitsche) has
  no set to act on. The phases separate silently: in the MB-001 reproducer with the tie-break, design strain energy was
  exactly 0, with no warning.
- **Reproducer:** MB-001 reproducer with fix v1 (the GEN tie-break), or any deck where a static geometry's boundary
  coincides with grid lines over whole cells adjacent to another non-cut phase.
- **Root cause:** not investigated. Interface facets appear to come only from cut cells. Check whether XTK should emit
  facets on background edges shared by non-cut cells of different bulk phases.
- **Workaround:** keep static geometry boundaries off grid lines. MB-001 v2 now errors in the worst case (all vertices
  on the interface); partially grid-aligned boundaries are not detected yet.

---

## Entry template
```
## MB-NNN · <module>: <one-line summary>
- **Status:** open | fix-written (<branch>, <date>) | fixed (<commit>) | merged | wontfix (<why>)
- **Symptom:** what the user sees; exact error text or log signature
- **Reproducer:** durable path (deck/run dir/unit test) + the minimal config that triggers it
- **Root cause:** file:line and the data condition (or "not investigated")
- **Fix:** what changed (commit), and why that is the right place
- **Regression test:** test file + tag
- **Workaround:** deck-level avoidance, if any
```
