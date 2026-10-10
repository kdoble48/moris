# MORIS bug ledger (kdoble48 fork)

One entry per MORIS defect found while using the code: crashes, wrong results, silent failures.
Missing features belong in a study's capability list, not here (exception: a missing piece that fails
*silently*, which is a bug). The ledger is versioned with the code, so a fix commit and its entry move together.

**Rules**
- New find → add an entry with the next `MB-NNN`, status `open`, before or while fixing it.
- Every fix commit message cites its ID (`fix(XTK): ... [MB-001]`) and the entry records the commit and the regression test.
- Status: `open` → `fix-written` (on a branch, not yet verified) → `fixed` (built, regression test + original
  reproducer pass) → `merged` (on `develop` of the fork). `wontfix`/`not-a-bug` with a reason.
- Reproducers point at durable paths (`runs/...` decks, unit tests). **Evidence bundles** (pathology figure from the real output, log excerpt, copies of the small raw artifacts) live in `runs/moris_bugs/MB-NNN/` and are rebuilt by `studies/lunch_and_learn_xiga/scripts/build_bug_evidence.py`; add a builder function for each new bug. Workspace agents read this file before triage
  (`.claude/rules/moris-bug-ledger.md`).

## Index

| ID | Status | Module | Summary | Found by |
|---|---|---|---|---|
| MB-001 | fixed (992639917) | XTK | Segfault building the enriched integration mesh when a cell's vertices all lie on one geometry's interface → now a clear error | lunch_and_learn_xiga B13 |
| MB-002 | fixed (22c3c2799; residual = MB-011) | FEM/INT | Double-sided FD geometry sensitivity uses the leader's parametric coords for the follower → dR/dp ~1e8 on interfaces along cell edges | lunch_and_learn_xiga B12 |
| MB-003 | open | GEN/FEM | Static user-defined geometry with `discretization_mesh_index = -2` (analytic) makes design gradients NaN, silently | moment_fitting (2026-07-16) |
| MB-004 | open | FEM | Perimeter (interface VOLUME IQI) adjoint off by 12–150% on a few small-gradient ADVs | lunch_and_learn_xiga B11 |
| MB-005 | open (root-caused; native-GCMMA workaround) | OPT/HMR | Adaptive refinement during optimization: ADV count drifts between remeshes → abort in `compute_design_criteria_gradients` | seeding_impact (2026-07-13) |
| MB-006 | open (suspect) | FEM/INT | Nonconformal branch of `select_dRdp_FD_geometry_double` also copies the leader ξ into the follower row | review of MB-002 |
| MB-007 | open | FEM/VIS | Requested IQI evaluated on zero sets writes zeros to Exodus with no warning | Hertz contact (2026-07) |
| MB-008 | open (limitation) | FEM/INT | `IQI_Max_Stress::compute_dQIdu` stops with "not implemented"; decks must switch the whole forward analysis to FD | moris_stress_demo (2026-10-06) |
| MB-009 | fixed (a856e9c13) | GEN | Design factory rejects a level-set geometry listed before the fields it depends on ("does not recognize geometry type"); properties built on the wrong field | MB-001 regression test |
| MB-010 | open | XTK | Phase boundary along background grid lines creates no interface facets → phases silently decoupled | MB-001 verification |
| MB-011 | fixed (06515ef67) | FEM/INT | Cluster-measure FD sensitivities (Nitsche-interface SP) use a base residual from CM/SP caches of the last perturbed geometry → 0.5–15% gradient error next to interfaces | MB-002 verification |
| MB-012 | open | VIS | STANDARD VIS mesh writes non-finite nodal values (UX/UY/VM) at 11–15% of nodes, mostly void-side but some with φ ≥ 0 | X18 peak extraction |
| MB-013 | open (limitation) | WRK/OPT | Every-step in-place SDF redistancing edits the ADV vector between criteria and gradient evaluation → ADV-consistency `MORIS_ERROR` after O(100) iterations (native MMA only; GCMMA cannot use it) | sdf_reinit_sweep (2026-06) |
| MB-014 | open | WRK | Scheduled-restart SDF apply in `Workflow_HMR_XTK::initialize()` replaces the whole ADV vector with one field's coefficients → silently skipped ("does not match current ADV count") on composite designs (level set + PDVs) | sdf_reinit_sweep (2026-06), review 2026-10-10 |
| MB-015 | open | WRK | `sdf_reinit_mode = "field_remap"` throws `std::out_of_range` in `Reinitialize_Performer::perform()` | sdf_reinit_sweep (2026-06) |
| MB-016 | fixed (f6d21bfb7) | OPT | Native in-tree GCMMA (`algorithm = "mma"`) wrote no `ADV_Alg_0_Iter_N.hdf5` restart files: no warm restart, no post-run tail-velocity check, silently | native-GCMMA validation 2026-10-10 |
| MB-017 | fixed (7a3115ee7) | GEN/WRK | Every-step in-place SDF redistance wrote coefficients at their global ADV IDs instead of ADV-vector positions → design shifted by one coefficient per redistance ("wrote 4960 of 4961"), MBB cut off its load pad | warm-restart work 2026-10-10 (D-015) |
| MB-018 | fixed (a0fe9502e) | GEN | `Stored_Field::get_mtk_field()` built an MTK field on a null mesh pair → segfault in `mtk::Field::Field` whenever SDF reinit (or the remeshing restart) asked GEN for its MTK fields with a static stored geometry | L-bracket SDF port 2026-10-10 |
| MB-019 | fixed (a0fe9502e) | GEN/WRK | `SIGNED_DISTANCE_IMAGE` field ignored the `name` parameter (named after the image file), so SDF reinit found no design: redistance skipped silently every iteration | L-bracket SDF port 2026-10-10 |
| MB-016 | open | XTK | Per-geometry cut-cell selection uses only background corner nodes; an analytic field evaluated at Box-cut child vertices splits a shared edge on one side only → nonconforming IG facets → SIGSEGV in `construct_subphase_neighborhood` (opt) | gaussian_splat_ls G-N24 (D-012) |
| MB-017 | fix-written (feature/gen-additive @ee83435b7) | GEN/MRS | `Design_Variable` ID is a `char` from a static counter; `ADV_Manager::create_adv` aliases any non-constant variable whose ID repeats → >~250 analytic Design_Variables silently share ADVs | gaussian_splat_ls G-N60 (D-011) |
| MB-020 | fixed (abb99ee66) | GEN | Multilinear intersection divides by the norm of the four corner φ values; a cell whose corners are all exactly 0 (order-2 field nonzero only at mid-side/interior nodes) gives NaN residuals → abort "Newton did not converge: Reference residual nan" | lunch_and_learn_xiga X40_cheese_o2 (D-020) |
| MB-021 | fixed (4ed1e2393) | OPT | Native GCMMA in-place SDF mode: the "redistanced checkpoint" rewrite (3ebdd8050) lands in the NEXT iteration's file and is then overwritten by the raw trial point, so Iter_0..N-1 hold ψ, never φ = R(ψ); only the last file holds φ | lunch_and_learn_xiga theory manual audit (X9r_rcp, X81_o2_cellorder) |
| MB-022 | fix-written (fix/mb022-025-audit @307236d81) | GEN | Nested intersection recursion drops the incoming sensitivity factor | implementation audit 2026-10-10 (D-028) |
| MB-023 | fix-written (fix/mb022-025-audit @307236d81) | OPT | Native GCMMA accepts NaN criteria in the conservativeness inner loop | implementation audit 2026-10-10 (D-029) |
| MB-024 | fix-written (fix/mb022-025-audit @307236d81) | OPT | Exhausted MMA backtracking leaves Problem on the rejected trial | implementation audit 2026-10-10 (D-030) |
| MB-025 | fix-written (fix/mb022-025-audit @307236d81) | OPT | SWEEP single finite-difference modes silently skip all perturbations | implementation audit 2026-10-10 (D-031) |

---

## MB-001 · XTK segfault: cell with every vertex on an interface
- **Status:** fixed, commit 992639917 on `fix/xiga-deck-bugs` (2026-10-09). Verified: the regression test passes and the original reproducer now stops with the new error instead of a SIGSEGV (`runs/studies/lunch_and_learn_xiga/MB001_repro_v2`); the full XTK/GEN/INT suites pass.
- **Evidence:** `runs/moris_bugs/MB-001/` — `v2_error_excerpt.txt` (the clean error), `why_not_a_tiebreak_decoupled_pads.png`.
- **Symptom:** SIGSEGV "Address not mapped" in `xtk::Enriched_Integration_Mesh::setup_blockset_with_cell_clusters()` at GCMMA iteration 1, right after "XTK - Enriched Integration Mesh - Construction". Opt build only.
- **Reproducer:** `studies/lunch_and_learn_xiga/decks/mbb_levelset.py`, MBB with `"elements": [60, 20]` and on-grid pads (no `pad_offgrid`), config `runs/studies/lunch_and_learn_xiga/configs/MB001_repro.json`. The two fixed pads are exactly 2×1 cells with every face on a grid line.
- **Root cause:**
  1. GEN `get_element_phase_index` returns `MORIS_INDEX_MAX` when every vertex of a cell has |φ| ≤ tolerance for one geometry.
  2. `XTK/src/cl_XTK_Integration_Mesh_Generator.cpp` `deduce_ig_cell_bulk_phase_from_facets` returned the sentinel `MORIS_INDEX_MAX - 1` for QUAD/HEX cells, which counted as assigned.
  3. It was then used unchecked as an index at `cl_XTK_Enriched_Integration_Mesh.cpp` (`tNoChildBlockSetOrds( tBulkPhaseIndex )`).
- **Fix, v1 (rejected 2026-10-09):** a GEN tie-break treating such cells as positive. The run then completed, but the pad was **silently decoupled**: design strain energy exactly 0, and no `dbl_iside_p0_1_p1_2` set processed. A phase boundary along grid lines gets no interface facets (MB-010), so v1 traded a crash for a wrong answer.
- **Fix, v2:** XTK stops with a `MORIS_ERROR` naming the cell and the remedy (shift the geometry off the grid lines). A `MORIS_ERROR` guard also sits before the block-set index. GEN is unchanged.
- **Code** (commit 992639917):
  ```diff
  // projects/XTK/src/cl_XTK_Integration_Mesh_Generator.cpp — deduce_ig_cell_bulk_phase_from_facets()
   if ( (uint)aCell->get_geometry_type() == 4 || (uint)aCell->get_geometry_type() == 2 )
   {
  -    return MORIS_INDEX_MAX - 1;      // sentinel != MORIS_INDEX_MAX, so the cell counted as "assigned"
  +    MORIS_ERROR( false,
  +            "Integration_Mesh_Generator::deduce_ig_cell_bulk_phase_from_facets() - every vertex of background cell %d "
  +            "(id %d) lies on a geometry interface, so its bulk phase cannot be determined. [...] "
  +            "Shift the geometry off the grid lines (e.g. by half a cell).",
  +            aCell->get_index(), aCell->get_id() );
   }

  // projects/XTK/src/cl_XTK_Enriched_Integration_Mesh.cpp — setup_blockset_with_cell_clusters()
   if ( tPrimarySubPhase > -1 )
   {
  +    MORIS_ERROR( tBulkPhaseIndex >= 0 && tBulkPhaseIndex < (moris_index)tNoChildBlockSetOrds.size(),
  +            "... enriched IP cell %d has invalid bulk phase index %d (number of bulk phases: %d).", ... );
       if ( tCluster.is_trivial() )
       {
           tSetOrd = tNoChildBlockSetOrds( tBulkPhaseIndex );   // was the segfault (index 2147483646)
  ```
  Upstream of it, GEN reports the undetermined phase (unchanged; `projects/GEN/src/cl_GEN_Geometry_Engine.cpp:546`):
  ```cpp
              // All vertices are (somehow) on the interface; return that this is a problem
              return MORIS_INDEX_MAX;
  ```
- **Regression test:** `XTK/test/xtk/UT_Snapping.cpp` "Band One Cell Thick on Grid Lines 2d" `[Snapping_2d]`, expecting `CHECK_THROWS( perform_decomposition() )`.
- **Workaround:** keep fixed-geometry faces off grid lines.

## MB-002 · FEM: follower parametric coordinates in double-sided FD geometry sensitivities
- **Status:** fixed, commit 22c3c2799 (2026-10-09); the 1e8 defect is gone. The remaining 0.5–15% error near the interface is a separate defect (MB-011), also present in the off-grid pad25 control. After the fix the reproducer gives ADVs 262/342/343/349 = −6.58e-3/−6.64e-3/−2.64e-2/−6.93e-4 vs FD −6.93e-3/−7.00e-3/−2.68e-2/−6.21e-4: the 1e8 spikes are gone, signs and magnitudes are right, but 1.7–12% error remains on the grid-aligned interface row. The residual is under investigation (candidates: follower point leaving its element during the FD step, ghost facets on cell edges, the step size chosen from the leader only).
- **Evidence:** `runs/moris_bugs/MB-002/` — `gradient_vs_fd.png` (before/after bars), `bridge_strain_energy_history.png` + `bridge_v1_designs_*.png` (the optimization pathology), raw SWEEP hdf5 ×3, deck + config.
- **Symptom:** analytic d(strain energy)/dADV of order 1e7–1e8, sometimes with the wrong sign, for ADVs whose intersection nodes sit on a double-sided interface facet lying along a background-cell edge. Volume and perimeter gradients are fine, because they use explicit FD and don't go through λᵀ∂R/∂p. No NaN and no warning; GCMMA climbs and topology degrades.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/triage_bridge/sweep_pad20_idx1/`, a bridge with the deck pad's top face on node row J=2. SWEEP ADVs 262/342/343/349: analytic −1.25e8/+1.31e8/+6.9e7/−6.4e6 vs central FD −6.93e-3/−7.00e-3/−2.69e-2/−6.2e-4. The control with the pad face at 2.5h (`sweep_pad25_an`) is clean.
- **Root cause:** `FEM/INT/src/IWG/cl_FEM_IWG.cpp` `IWG::select_dRdp_FD_geometry_double` computed the perturbed node's ξ in the *leader* IP element and wrote it into the *follower* IG cell's parametric coefficients. That is only valid when both sides share one IP element. Across a cell edge (η = −1 vs +1) the follower point jumped a whole cell, giving an O(1) residual jump divided by δ ≈ 5.6e-9. GEN's dx/dφ was verified correct (the ADV ratio test matches, 1.90).
- **Fix:** for conformal sets, the follower ξ comes from the follower IP geometry interpolator (`update_parametric_coordinates`). The nonconformal path is unchanged (see MB-006).
- **Code** (commit 22c3c2799, `projects/FEM/INT/src/IWG/cl_FEM_IWG.cpp`, `IWG::select_dRdp_FD_geometry_double`):
  ```diff
   // xi of the perturbed node, computed in the LEADER IP element
   tLeaderIPGI->update_parametric_coordinates( tXCoords, tXiCoords );

   Matrix< DDRMat > tFollowerParamCoeffPert = tFollowerParamCoeff;
  -tFollowerParamCoeffPert.get_row( iFollowerNode ) = tXiCoords.matrix_data();     // leader xi reused for the follower
  +if ( tIsNonconformal )
  +{
  +    tFollowerParamCoeffPert.get_row( iFollowerNode ) = tXiCoords.matrix_data();  // unchanged, see MB-006
  +}
  +else
  +{
  +    // follower xi from the FOLLOWER IP element (differs when the facet lies on a background-cell edge)
  +    Matrix< DDRMat > tFollowerXCoords  = tFollowerCoeffPert.get_row( iFollowerNode );
  +    Matrix< DDRMat > tFollowerXiCoords = tFollowerParamCoeff.get_row( iFollowerNode );
  +    tFollowerIPGI->update_parametric_coordinates( tFollowerXCoords, tFollowerXiCoords );
  +    tFollowerParamCoeffPert.get_row( iFollowerNode ) = tFollowerXiCoords.matrix_data();
  +}
  ```
  With leader η = −1 and follower η = +1 on a horizontal cell edge, the old line moved the follower Gauss point one
  whole cell per FD perturbation (δ ≈ 5.6e-9), so dR/dp ≈ O(1)/δ ≈ 1e8.
- **Verification to run:** the reproducer above; pass = analytic within ~1e-4 relative of FD, and Criteria norms drop from 1e7–1e8 to O(1).
- **Regression test:** TODO, `FEM/INT/test/UT_FEM_IWG_Struc_Linear_Interface.cpp`, two-IP-element interface facet vs hand FD.
- **Was misdiagnosed first** as "fixed geometry on the grid" and as the −2 static path (MB-003). The −1 re-test ruled the latter out.

## MB-003 · GEN: analytic static user geometry poisons design gradients with NaN
- **Status:** open (deck workaround known).
- **Evidence:** `runs/moris_bugs/MB-003/` — none yet (logs not found; README says how to rebuild it).
- **Symptom:** once a static (ADV-free) user geometry's boundary becomes an exposed material interface, every `GCMMA - Criteria(N): norm = -nan`. The design freezes, with no error or warning.
- **Reproducer:** moment_fitting p2 padded warm restart (`runs/studies/moment_fitting/phase1_p2_opt/deck/boxbeam3d_frameless_lean_p2.py`, pads with index −2). A/B: −1 is clean.
- **Root cause (partial):** the analytic user-field path (`cl_GEN_User_Defined_Field`: FD `get_dfield_dcoordinates`, no `get_dfield_dadvs`) produces NaN in the intersection-node sensitivity chain.
- **Code** (`projects/GEN/src/ADV/fields/cl_GEN_User_Defined_Field.cpp`): the analytic (−2) path differentiates the user
  function by an internal FD with a fixed 1e-8 step and has no ADV sensitivity for ADV-free geometries. The NaN appears
  somewhere in this chain (exact site not yet located):
  ```cpp
  void User_Defined_Field::get_dfield_dcoordinates( const Matrix< DDRMat >& aCoordinates, Matrix< DDRMat >& aSensitivities )
  {
      const real tPertubation = 1e-8;
      Matrix< DDRMat > tCoordinates( aCoordinates );
      for ( uint idim = 0; idim < aCoordinates.numel(); ++idim )
      {
          tCoordinates( idim ) += tPertubation;
          real tFieldValue = this->get_field_value_user_defined( tCoordinates, mFieldVariables );
          tCoordinates( idim ) -= 2.0 * tPertubation;
          ...
  ```
- **Notes:** pymoris `add_user_geometry` defaults to −2, so every pymoris deck with fixed pads is exposed. Either fix GEN or make −2 error out for geometries without ADVs in optimization runs.
- **Workaround:** `discretization_mesh_index = -1`. Memory note: gen-static-geometry-stored-field.

## MB-004 · FEM: perimeter adjoint mismatch
- **Status:** open. **Re-test after MB-011 (2026-10-10, queue job 091, `runs/studies/lunch_and_learn_xiga/X6_sweep/mb004_retest/`): not fixed by it.** Large-gradient ADVs match FD to 1e-5..2e-3 (416, 1104, 1112, 1113, 3812, 4014); small-gradient ADVs (|g| <= 6e-4) are off by 14-552%, three with the wrong sign (417 +6.02e-4 vs FD +5.27e-4; 1103 +5.42e-5 vs -2.73e-5; 1102 +7.07e-5 vs -1.57e-5; 4016 -5.24e-5 vs +3.13e-5). Volume matches to 1e-8 everywhere.
- **Evidence:** `runs/moris_bugs/MB-004/` — `perimeter_gradient_table.md` (pre-MB-011 mismatches), raw SWEEP, post-MB-011 diagnostic logs.
- **Symptom:** X6 SWEEP at the MBB seed. Strain-energy and volume adjoints match central FD (≤ 6.5e-5, ~1e-8). The perimeter criterion (`IQI VOLUME` on `iside_*` sets) is off by 12%, 18% and 150% on 3 of 12 ADVs. FD is stable to 4 digits over steps 1e-5 to 1e-7. The absolute size is ~1e-4 of the largest perimeter gradient.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X6_sweep/` (runB_fd, ADVs 417, 4015, 1103).
- **Root cause:** not investigated.

- **Code:** defect not located yet; add the excerpt when the root cause is found.
## MB-005 · OPT/HMR: adaptive refinement + optimization ADV drift
- **Status:** open.
- **Evidence:** `runs/moris_bugs/MB-005/` — `crash_excerpt.txt` (the abort with matrix sizes).
- **Symptom:** abort with "incompatible matrix dimensions 1x4370 and 1x4374" in `opt::Problem::compute_design_criteria_gradients`. The B-spline ADV count changes between remeshes because the refinement pattern follows the evolving topology.
- **Code** (`projects/OPT/src/cl_OPT_Problem.cpp:276-287`): the product fails when the ADV count of `get_dcriteria_dadv()`
  changed after a remesh, while `mADVs` and the user callbacks still use the old count:
  ```cpp
  mObjectiveGradient =
          this->compute_dobjective_dadv() +
          this->compute_dobjective_dcriteria() * mInterface->get_dcriteria_dadv();   // 1x4370 * ... vs 1x4374
  ```
- **Reproducer:** seeding_impact deck family at 60×20 + refine-1 with remeshing (`studies/seeding_impact/FINDINGS.md`, S3-W2 it.70).
- **Root cause (2026-10-10, read from code + the archived S4-W3 log, `runs/studies/seeding_impact/_archive_v1_remeshing/runs/S4-W3/logs/03_running_moris.log` lines 37526-39887):**
  a double failure. A remesh restart fires, then the first evaluation on the new mesh fails the XTK same-level check
  again (`check_intersected_background_cell_levels`, `cl_XTK_Integration_Mesh_Generator.cpp:407-430` → NaN + restart flag).
  The TPL GCMMA wrapper calls the gradient anyway, and `Criteria_Interface::get_dcriteria_dadv` (`cl_OPT_Criteria_Interface.cpp`)
  returns the cached pre-remesh `mSensitivities` while the restart flag is set (1x4260 vs 1x4475). The native GCMMA
  (`cl_OPT_Algorithm_MMA.cpp:917`) checks `restart_optimization()` right after the first `func` and exits before `grad`.
- **Workaround:** a static uniform mesh, or the native GCMMA (`Optimization_Algorithm_Type::MMA`): verified in
  `runs/studies/lunch_and_learn_xiga/X50_adaptive_smoke` (adaptive refinement + in-place interface SDF + multilinear, one
  reactive remesh at eval 8, ADVs 13073 → 13246, no abort).
- **Proposed fix:** `get_dcriteria_dadv` raises `MORIS_ERROR` when the restart flag is set and the cached column count differs
  from the current ADV count; the TPL GCMMA wrapper checks `restart_optimization()` after its first evaluation, like
  `Algorithm_MMA::mma_solve`; `Manager::reinitialize` counts consecutive restarts with no successful evaluation (error after 3)
  and loops instead of recursing.

## MB-006 · FEM: nonconformal FD geometry sensitivity copies the leader ξ (suspect)
- **Status:** open (suspect, unverified).
- **Symptom:** none observed yet. Found by reading `select_dRdp_FD_geometry_double` while fixing MB-002. In the nonconformal branch the follower node is not perturbed, but its parametric row is overwritten with the leader's ξ before the ray remap.
- **Code** (`projects/FEM/INT/src/IWG/cl_FEM_IWG.cpp`, nonconformal branch kept as-is by the MB-002 fix): the follower
  coordinates are not perturbed for nonconformal sets, but its parametric row is still overwritten with the LEADER's xi:
  ```cpp
  if ( tIsNonconformal )
  {
      // unchanged behavior for nonconformal sets (follower not perturbed; see the remapping below)
      tFollowerParamCoeffPert.get_row( iFollowerNode ) = tXiCoords.matrix_data();
  }
  ```
- **Owner:** whoever owns the nonconformal contact work (`feature/nonconformal-*`). Check whether the remap later overrides it; if not, this is MB-002's twin.

## MB-007 · FEM/VIS: IQI on zero sets writes silent zeros
- **Status:** open (usability).
- **Symptom:** a contact IQI (`GAP`, `CONTACT_PRESSURE_*`) missing `IQI_bulk_type = NONCONFORMAL_SIDESET` is matched to 0 sets (`nIQIonSet=0`). VIS writes uninitialized zeros and nothing warns.
- **Code:** defect not located yet; add the excerpt when the root cause is found.
- **Suggested fix:** warn or error when a VIS-requested IQI is evaluated on no set. Memory note: nonconformal-contact-iqis-zero-at-vis.

## MB-008 · FEM: no analytic dQ/du for the p-norm stress IQI
- **Status:** open (limitation; it fails loudly, so it is kept here only because it forces whole-model FD).
- **Symptom:** `IQI_Max_Stress::compute_dQIdu` calls `MORIS_ERROR("not implemented")` (`cl_FEM_IQI_Max_Stress.hpp:99-114`). Stress-constrained decks must set `is_analytical_forward = false`, which is slow in 3D.

- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Max_Stress.hpp:99-114`):
  ```cpp
  void compute_dQIdu( real aWStar ) override
  {
      MORIS_ERROR( false, "IQI_Max_Stress::compute_dQIdu - not implemented." );
  }
  ```
## MB-009 · GEN: design factory ignores unbuilt dependencies
- **Status:** fixed, commit a856e9c13 on `fix/xiga-deck-bugs` (2026-10-09); GEN suite passes (80101 assertions).
- **Symptom:** `GEN does not recognize geometry type: level_set` (`cl_GEN_Design_Factory.cpp`) when a geometry with
  `dependencies` comes before its fields, e.g. a COMBINED_FIELDS geometry in GEOMETRIES with its fields in the FIELDS
  sublist, which `Geometry_Engine` appends *after* the geometries.
- **Reproducer:** the MB-001 regression test as first written (fields in `tGENParams( 2 )`).
- **Root cause:** the design-building block ran even when `tCanBuild` was false. The level-set branch tested
  `level_set and tCanBuild`, so an unbuildable level set fell through to the error. The property branch had no check
  and would bind `mFields( tFieldIndex - 1 )` (the wrong field). Both erased `design_type` even when nothing was built.
- **Fix:** build designs only `if ( exists( "design_type" ) and tCanBuild )`; the existing outer loop retries.
- **Code** (commit a856e9c13, `projects/GEN/src/ADV/cl_GEN_Design_Factory.cpp`):
  ```diff
  -// Check if design needs building
  -if ( iParameterList.exists( "design_type" ) )
  +// Check if design needs building (only once its field dependencies exist; otherwise the outer loop retries it)
  +if ( iParameterList.exists( "design_type" ) and tCanBuild )
   {
       ...
  -    if ( tGeometryType == "level_set" and tCanBuild )      // false -> fell through to the error below
  +    if ( tGeometryType == "level_set" )
       { ... }
       ...
       else
       {
           MORIS_ERROR( false, "GEN does not recognize geometry type: %s", tGeometryType.c_str() );
       }
  ```
- **Regression test:** `GEN/test/ut_GEN_Field.cpp` "Combined Field Declared Before Its Fields".

## MB-010 · XTK: phase boundary along background grid lines is not coupled
- **Status:** open (found 2026-10-09 while verifying MB-001 v1).
- **Evidence:** `runs/moris_bugs/MB-010/` — `decoupled_pads_ux.png` (design carries nothing), `v1_log_excerpt.txt` (SE design = 0, no dbl_iside set), raw exodus.
- **Symptom:** two phases meet along background-cell edges (e.g. a fixed pad whose faces lie exactly on grid lines,
  next to non-cut design cells). No double-sided interface facets are created there, so an interface IWG (Nitsche) has
  no set to act on. The phases separate silently: in the MB-001 reproducer with the tie-break, design strain energy was
  exactly 0, with no warning.
- **Reproducer:** MB-001 reproducer with fix v1 (the GEN tie-break), or any deck where a static geometry's boundary
  coincides with grid lines over whole cells adjacent to another non-cut phase.
- **Root cause:** not investigated. Interface facets appear to come only from cut cells. Check whether XTK should emit
  facets on background edges shared by non-cut cells of different bulk phases.
- **Code:** defect not located yet; add the excerpt when the root cause is found.
- **Workaround:** keep static geometry boundaries off grid lines. MB-001 v2 now errors in the worst case (all vertices
  on the interface); partially grid-aligned boundaries are not detected yet.

## MB-011 · FEM: stale CM/SP caches in cluster-measure FD geometry sensitivities
- **Status:** fixed, commit 06515ef67 (2026-10-10). Falsifier (unfixed binary, central FD): 262/342/343/349 within 8.5e-5/1.4e-4/4.4e-6/9.6e-4 of FD, as
- **Evidence:** `runs/moris_bugs/MB-011/` — `gradient_vs_fd.png`, falsifier + verification logs, pad25 control SWEEP before/after.
  predicted. After the fix (default forward FD): the same, and the pad25 control is within 1e-5..1e-3 (was 15%). INT suite passes.
  Logs: `runs/studies/lunch_and_learn_xiga/queue/logs/005_mb011_falsify.log`, `007_build_mb011.log`.
- **Symptom:** the analytic objective gradient is off by 0.5–15% on ADVs next to a Nitsche interface whose stabilization
  depends on cluster measures; the far field is fine. Deterministic, independent of the FD step size, symmetric.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/triage_bridge/sweep_pad25_fd` (ADVs 343/349/397: 0.5%/15%/1.0%) and
  `sweep_verify_mb002` (262/342/343/349: 5%/5%/1.7%/12%).
- **Root cause (static analysis):** `IWG::select_dRdp_FD_geometry_double` / `_geometry` restore the coordinates, then
  call `add_cluster_measure_dRdp_FD_geometry[_double]`. That routine computes its "unperturbed" residual without
  `reset_eval_flags()`, so the CM flux/traction and SP caches (`cl_FEM_Constitutive_Model.cpp` mFluxEval/mTractionEval)
  still hold the last perturbed state. With the default one-sided (forward) FD scheme, that stale base enters
  (∂R/∂x)·(h/3)/MEA·dMEA/dp for every PDV of the cluster. `IQI::add_cluster_measure_dQIdp_FD_geometry` has the same pattern.
  The elementwise bulk paths already reset (`cl_FEM_IWG.cpp` ~3422, `cl_FEM_IQI.cpp` ~2669).
- **Fix:** `this->reset_eval_flags()` before the base evaluation in the two IWG routines and the IQI routine.
- **Code** (commit 06515ef67): the same insertion in `IWG::add_cluster_measure_dRdp_FD_geometry`,
  `IWG::add_cluster_measure_dRdp_FD_geometry_double` (`cl_FEM_IWG.cpp`) and `IQI::add_cluster_measure_dQIdp_FD_geometry`
  (`cl_FEM_IQI.cpp`). The old comment already said "reset", but nothing reset:
  ```diff
   // reset, evaluate and store the residual for unperturbed case
  +this->reset_eval_flags();     // CM/SP caches still held the last perturbed geometry
   mSet->get_residual()( 0 ).fill( 0.0 );
   this->compute_residual( aWStar );
  ```
- **Workaround:** FEM computation `finite_difference_scheme = POINT_3_CENTRAL` (no base point).
- **Regression test:** TODO, forward vs central dR/dp on a Nitsche-interface facet with an active PDV (agree to 1e-6).

## MB-012 · VIS: non-finite nodal output values
- **Status:** open (found 2026-10-10).
- **Evidence:** `runs/moris_bugs/MB-012/` — `nonfinite_nodes_x9.png` (non-finite nodes trace every hole boundary), raw VIS exodus + GEN field.
- **Symptom:** nodal IQI fields in the STANDARD VIS exodus (UX, UY, VM) are non-finite at 11–15% of nodes: X9 MBB
  917/5946, stress demo 409/3735, X18 1605/14667. Of X9's 917, 793 lie in the void (φ < 0) but 124 have φ ≥ 0.
  Post-processing that takes a max without masking returns NaN; this broke the first X18 σ* calibration
  (σ* = 0.6 × NaN made design 2 abort with an infinite P-norm).
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X9_baseline/results/exodus/mbb_levelset.exo` (last step) with the
  GEN level set `GEN_mbb.e-s.0151`.
- **Root cause:** not investigated. Possibly nodal averaging with zero weight for nodes of material clusters that no
  material integration cell touches. Check whether non-finite output is intended as a "no value" marker; if it is, it
  should be documented and kept out of material nodes.
- **Code:** defect not located yet; add the excerpt when the root cause is found.
- **Workaround:** mask non-finite values in post-processing (`scripts/lbracket_peak.py`).

---

## MB-013 · WRK/OPT: every-step in-place SDF redistancing trips the ADV-consistency check
- **Status:** open (limitation; recorded 2026-10-10 from the 2026-06 sweep).
- **Symptom:** with `sdf_reinit_every_step = true` (native MMA, `reinit_in_place = true`) the run aborts after O(100)
  iterations (iteration ~829 in IP_f10): `Problem::compute_design_criteria_gradients - given ADV vector does not match
  ADVs used to compute design criteria.` (`cl_OPT_Problem.cpp:265`). Before the abort the objective oscillates
  (0.398–0.55 vs baseline 0.529).
- **Reproducer:** `runs/studies/sdf_reinit_sweep/IP_f10/` (H18 60×20 box beam, every-step interface reinit, cadence 10).
- **Evidence:** `runs/moris_bugs/MB-013/` — `assert_excerpt_IP_f10.txt`, `baseline_vs_reinit10.png`, the sweep's
  `FINDINGS.md`, the deck `IP_f10.cpp`.
- **Root cause:** the in-place block in `cl_WRK_Workflow_HMR_XTK.cpp` (`perform()`, "Every-step in-place redistance",
  ~lines 480–530) overwrites `aNewADVs` after the criteria evaluation, so the gradient call sees a different ADV vector
  than the one the criteria were computed for. Documented in-code as KNOWN LIMITATION (commit d50ede932). GCMMA (TPL)
  cannot use the mode at all.
- **Code:**
  ```cpp
  // projects/OPT/src/cl_OPT_Problem.cpp:263-267 — the check that fires
  MORIS_ERROR( std::abs( mADVs( iADVIndex ) - aNewADVs( iADVIndex ) ) <= mADVNormTolerance,
          "Problem::compute_design_criteria_gradients - given ADV vector does not match ADVs used to compute design criteria.\n" );
  // projects/WRK/src/cl_WRK_Workflow_HMR_XTK.cpp (perform(), every-step block) — the write between the two calls
  aNewADVs( tAdvId ) = tCoeffs( iC );
  ```
- **Workaround:** use scheduled-restart mode (`sdf_reinit_every_step = false`, redistance every N iterations through a
  GCMMA restart). Valid for single-field designs; see MB-014 for composite ones.

## MB-014 · WRK: scheduled-restart SDF apply assumes the level set is the whole ADV vector
- **Status:** open (recorded 2026-10-10).
- **Symptom:** on a restart after a scheduled redistance, `initialize()` logs `SDF coefficient count (%d) does not match
  current ADV count (%d) ... Skipping SDF application.` and the redistance is silently dropped whenever the ADV vector
  also holds other designs or PDVs (H18: 524 level-set coefficients vs 3767 ADVs, so 100% skipped).
- **Reproducer:** the H18 composite deck of `runs/studies/sdf_reinit_sweep/` with `sdf_reinit_every_step = false`.
- **Evidence:** `runs/moris_bugs/MB-013/sdf_reinit_sweep_FINDINGS.md` ("Layer 6", "What a complete study needs" 1).
- **Root cause:** the in-place path was converted to a slice write (d50ede932: coefficients written at
  `get_design_determining_adv_ids( adv_field )`), but the restart path still builds a full-length ADV vector from the
  coefficients.
- **Code:**
  ```cpp
  // projects/WRK/src/cl_WRK_Workflow_HMR_XTK.cpp, Workflow_HMR_XTK::initialize() (~lines 131-152)
  if ( tStoredSize != tExpectedSize ) { MORIS_LOG_WARNING( "SDF coefficient count (%d) does not match current ADV count ..." ); }
  else
  {
      Vector< real > tNewADVs( tCoeffs.n_rows() );          // whole vector = one field's coefficients
      for ( uint i = 0; i < tCoeffs.n_rows(); ++i ) { tNewADVs( i ) = tCoeffs( i ); }
      mPerformerManager->mGENPerformer( 0 )->set_advs( tNewADVs );
  }
  ```
- **Fix (proposed):** the same slice write as the in-place path, onto a copy of the current GEN ADVs.
- **Workaround:** keep the level set as the only design field (no PDVs, static geometries with
  `discretization_mesh_index = -1`), as `studies/lunch_and_learn_xiga/decks/mbb_levelset.py` does.

## MB-015 · WRK: `field_remap` reinitialization throws `std::out_of_range`
- **Status:** open (recorded 2026-10-10 from the 2026-06 sweep; not re-reproduced).
- **Symptom:** `std::out_of_range` from `Reinitialize_Performer::perform()` with `sdf_reinit_mode = "field_remap"` (the default mode).
- **Reproducer:** the `FR_f*` arms of `runs/studies/sdf_reinit_sweep/` (H18 deck, field_remap mode).
- **Evidence:** `runs/moris_bugs/MB-013/sdf_reinit_sweep_FINDINGS.md` ("Layer 2").
- **Root cause:** not investigated.
- **Code:** defect not located yet; mode selection at `projects/WRK/src/cl_WRK_Reinitialize_Performer.cpp:75` and :125.
- **Workaround:** use `sdf_reinit_mode = "interface"`.

## MB-016 · XTK: analytic-field cut-cell selection from background corners → nonconforming child meshes → SIGSEGV
- **Status:** open (2026-10-10)
- **Symptom:** SIGSEGV "Address not mapped", failing address = <base>+0x1fffffffc ((2^31−1)×4), in `xtk::Integration_Mesh_Generator::construct_subphase_neighborhood` (+0xe8d) ← `Integration_Mesh_Generator::perform` ← `Model::decompose`, right after "XTK - Integration_Mesh_Generator - Creating Facets". Opt build; a debug build would stop at the `MORIS_ASSERT` "Adjacent IG Cell not found" (predicted, build_dbg predates the branch). Seen at GCMMA eval 7 of `runs/studies/gaussian_splat_ls/runs/G-N24`.
- **Reproducer:** `runs/moris_defects/D-012_G-N24/triage/repro_eval7.json` (crash θ `theta_G-N24_eval7.npy`, max_iter 1, ~10 s); minimal: `triage/synth_N1_crash.json` (one Gaussian, SEGV) vs `triage/synth_N1_control.json` (μ_y + 0.003, runs). Predictor `triage/predict_nonconformal.py` matches 24/24 drop-one runs.
- **Evidence:** `runs/moris_defects/D-012_G-N24/triage/` (logs, `drop/RESULTS.txt`, θ files); bank to `runs/moris_bugs/MB-016/`.
- **Root cause:** with ≥ 2 geometries, geometry k > 0 subdivides only background cells flagged by GEN `is_intersected` on the 4 background corner nodes (`projects/GEN/src/cl_GEN_Geometry_Engine.cpp:300-323`; list built at `projects/XTK/src/cl_XTK_Integration_Mesh_Generator.cpp:476`, consumed by `projects/XTK/src/cl_XTK_Node_Hierarchy_Interface.cpp:411-430`). An analytic field (`discretization_mesh_index = -2`) is evaluated exactly at derived vertices (`projects/GEN/src/ADV/fields/cl_GEN_Field_Analytic.hpp:104-108`), so a geometry-0 cut vertex on a shared background edge can carry the opposite geometry-k sign to both endpoints. Only the flagged neighbour splits that edge → an IG facet with one cell → neighbour index left at `MORIS_INDEX_MAX` → unchecked read. B-spline / linear fields cannot produce this (linear along the edge).
- **Code** (open; `projects/XTK/src/cl_XTK_Integration_Mesh_Generator.cpp:1409-1429`):
  ```cpp
  moris_index tNeighborSubphaseIndex = MORIS_INDEX_MAX;
  for ( const auto& iCell : aFacetConnectivity->mFacetToCell( tIgFacetIndex ) )
  {
      if ( iCell->get_index() != tMyIgCellIndex )
      {
          tNeighborSubphaseIndex = aCutIntegrationMesh->get_ig_cell_subphase_index( iCell->get_index() );
      }
  }
  MORIS_ASSERT( tNeighborSubphaseIndex != MORIS_INDEX_MAX,
          "Integration_Mesh_Generator::construct_subphase_neighborhood() - "
          "Adjacent IG Cell not found." );
  ...
  moris_index tNeighborBulkPhaseIndex = aCutIntegrationMesh->get_subphase_bulk_phase( tNeighborSubphaseIndex );   // SIGSEGV: index 2^31-1
  ```
- **Fix:** (a) select cells for geometry k > 0 from the current child-mesh vertices (or derive child-vertex regions consistently with the corner test); (b) promote the `MORIS_ASSERT` to `MORIS_ERROR` naming the background cell.
- **Regression test:** to add — XTK 2D, Box + analytic sliver through a Box-cut edge (synthetic θ above).
- **Workaround:** `use_multilinear_interpolation = true` on the analytic design geometry (tested: the crash design runs; changes strain energy +11.7% at a sliver-rich design → FD-check gradients before use). σ_min ≥ 2h only lowers the frequency.

## MB-017 · GEN/MRS: `Design_Variable` char ID wraps → analytic ADVs silently aliased
- **Status:** fix-written (feature/gen-additive @ee83435b7, 2026-10-10)
- **Symptom:** none at run time — a deck with > ~250 non-constant `Design_Variable`s gets fewer ADVs than variables; variables whose IDs repeat mod 256 share an ADV (value and bounds). Found by static review of the gaussian_splat_ls G-N60 deck (360 variables); D-011.
- **Reproducer:** `runs/studies/gaussian_splat_ls/runs/compile_G-N60/` (360 Design_Variables; compare the ADV count MORIS reports). Unit test `[mb017]` in GEN-test.exe (300 variables → 300 distinct ADVs).
- **Evidence:** `runs/moris_defects/D-011_compile_G-N60/`.
- **Root cause:** `char mID` assigned from `inline static char mCount` (`projects/MRS/CNT/src/cl_Design_Variable.hpp:23,30`, ctors `cl_Design_Variable.cpp:20,35`); `GEN::ADV_Manager::create_adv` returns an alias to an earlier ADV when the ID matches (`projects/GEN/src/ADV/cl_GEN_ADV_Manager.cpp:77-98`).
- **Code** (fix, `git show ee83435b7`, trimmed):
  ```diff
  diff --git a/projects/GEN/src/ADV/cl_GEN_ADV_Manager.hpp b/projects/GEN/src/ADV/cl_GEN_ADV_Manager.hpp
  --- a/projects/GEN/src/ADV/cl_GEN_ADV_Manager.hpp
  +++ b/projects/GEN/src/ADV/cl_GEN_ADV_Manager.hpp
  @@ -21,3 +21,3 @@ namespace moris::gen
         private:
  -        Vector< char > mParameterIDs;
  +        Vector< sint > mParameterIDs;
           bool mParameterIDsFinalized = false;
  @@ -47,3 +47,3 @@ namespace moris::gen
            */
  -        void register_parameter_ids( const Vector< char >& aParameterIDs );
  +        void register_parameter_ids( const Vector< sint >& aParameterIDs );
   
  diff --git a/projects/GEN/src/ADV/fields/fn_GEN_create_field.cpp b/projects/GEN/src/ADV/fields/fn_GEN_create_field.cpp
  --- a/projects/GEN/src/ADV/fields/fn_GEN_create_field.cpp
  +++ b/projects/GEN/src/ADV/fields/fn_GEN_create_field.cpp
  @@ -240,6 +240,6 @@ namespace moris::gen
   
  -    Vector< char > get_active_parameter_ids( const Parameter_List& aFieldParameterList )
  +    Vector< sint > get_active_parameter_ids( const Parameter_List& aFieldParameterList )
       {
           // Vector of parameter IDs
  -        Vector< char > tParameterIDs;
  +        Vector< sint > tParameterIDs;
   
  diff --git a/projects/MRS/CNT/src/cl_Design_Variable.cpp b/projects/MRS/CNT/src/cl_Design_Variable.cpp
  --- a/projects/MRS/CNT/src/cl_Design_Variable.cpp
  +++ b/projects/MRS/CNT/src/cl_Design_Variable.cpp
  @@ -51,3 +51,3 @@ namespace moris
   
  -    char Design_Variable::get_id() const
  +    sint Design_Variable::get_id() const
       {
  diff --git a/projects/MRS/CNT/src/cl_Design_Variable.hpp b/projects/MRS/CNT/src/cl_Design_Variable.hpp
  --- a/projects/MRS/CNT/src/cl_Design_Variable.hpp
  +++ b/projects/MRS/CNT/src/cl_Design_Variable.hpp
  @@ -22,3 +22,3 @@ namespace moris
         private:
  -        char mID;
  +        sint mID;
           real mValue;
  @@ -29,3 +29,3 @@ namespace moris
           // How many design variables have been created, used for having a unique ID
  -        inline static char mCount = 0;
  +        inline static sint mCount = 0;
   
  @@ -65,3 +65,3 @@ namespace moris
            */
  -        char get_id() const;
  +        sint get_id() const;
   
  ```
- **Fix:** `Design_Variable::mID`/`mCount`/`get_id()` widened `char` → `sint`, and the `Vector< char >` parameter-ID vectors in `ADV_Manager` / `get_active_parameter_ids` → `Vector< sint >`. Types only, no behaviour change; no serialized/MPI/HDF5 use of the ID. Commit ee83435b7 on `feature/gen-additive`.
- **Regression test:** `TEST_CASE "Design variable IDs beyond 256" [mb017]` in `projects/GEN/test/ut_GEN_Field.cpp` — 300 non-constant variables → 300 distinct ADVs with correct values/bounds (pre-fix: 256), plus same-object aliasing preserved.
- **Workaround:** keep ≤ 250 non-constant Design_Variables per deck (pymoris `GaussianSplatLevelSet` enforces `MAX_MORIS_ANALYTIC_ADVS = 250`).

---

## MB-016 · OPT: native GCMMA wrote no restart files
- **Status:** fixed, commit f6d21bfb7 on `fix/xiga-deck-bugs` (2026-10-10). Build green; OPT-test passes (32 assertions) with the new check.
- **Symptom:** a run with `Optimization_Algorithm_Type::MMA` (the native in-tree GCMMA port) leaves no `ADV_Alg_0_Iter_N.hdf5`
  at all. Nothing warns; warm restarts (`restart_file`) and post-run tail-velocity checks have nothing to read.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X27_native_cheese/` (30 its, `"algorithm": "mma"` in the MBB deck): zero
  checkpoint files; the TPL twin `X25_cheese` writes one per iteration.
- **Root cause:** `Algorithm_MMA::func()` (`projects/OPT/src/cl_OPT_Algorithm_MMA.cpp`) calls `compute_design_criteria` without
  `write_advs_to_file`, which the TPL wrapper (`cl_OPT_Algorithm_GCMMA.cpp:235`), ROL, LBFGS and SQP all call first.
- **Fix / Code** (commit f6d21bfb7):
  ```diff
           Vector< real > tADVs( mNumVar );
           for ( int i = 0; i < mNumVar; ++i ) tADVs( i ) = aX[ i ];
   
  +        // Write the restart file for this design, as the TPL GCMMA wrapper does before every
  +        // criteria evaluation (MB-016: the native port wrote none, so no warm restart was possible)
  +        this->write_advs_to_file( tADVs );
  +
           this->compute_design_criteria( tADVs );
  ```
- **Regression test:** `projects/OPT/test/opt_test.cpp`, section "MMA native": `REQUIRE( std::ifstream( "ADV_Alg_0_Iter_1.hdf5" ).good() )`.
- **Related tooling gap (fixed the same day, not MORIS):** native GCMMA logs `MMA - Objective:` where the TPL logs
  `GCMMA - Objective:`; the workspace merit gate, guard and post-run report matched only the TPL form, so native runs ran unguarded.

---

## MB-017 · GEN/WRK: in-place SDF redistance wrote at ADV IDs, not positions
- **Status:** fixed, commit 7a3115ee7 on `fix/xiga-deck-bugs` (2026-10-10).
- **Symptom:** with `sdf_reinit_every_step = true` the log says `wrote 4960 of 4961 level-set coeffs`; every redistance
  iterate sweeps 9–14% of the nodes (ordinary iterates ~3%); the MBB design is cut off its load pad after ~12 iterations
  and the strain energy freezes at the pad-only value (0.077352). Checkpoints across a redistance align best when shifted
  by one coefficient (sign flips 5.8% at shift −1 vs 11.1% unshifted; ordinary iterates align at shift 0).
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X28_inplace_recut/` and `X28_recut_f10/` (MBB cheese, `algorithm`
  `mma`, `sdf_reinit` `every_step`).
- **Evidence:** `runs/moris_bugs/MB-017/` — `detached_evals12_15.png`, `log_before.txt` / `log_after.txt`.
- **Root cause:** the slice write in `Workflow_HMR_XTK::perform()` (every-step block) used
  `Geometry_Engine::get_design_determining_adv_ids()`, which returns the B-spline field's GLOBAL ADV IDs (shared
  `ADV_Handler`, `cl_GEN_ADV_Handler.cpp:62`), as indices into `aNewADVs`. `set_advs()` places `aNewADVs( k )` at global ID
  `mFullADVIds( k )` (`cl_GEN_Geometry_Engine.cpp:149`); serially the IDs are one past the positions. The scheduled-restart
  path writes the coefficients in order and was not affected. Likely also behind the June sweep's "surface roughening /
  frayed members" in its in-place arms.
- **Fix / Code** (commit 7a3115ee7):
  ```diff
  // projects/WRK/src/cl_WRK_Workflow_HMR_XTK.cpp, every-step in-place block
  -                        Vector< sint >          tLSAdvIds = mPerformerManager->mGENPerformer( 0 )
  -                                ->get_design_determining_adv_ids( tReinitPerformer->get_adv_field_name() );
  +                        // positions in aNewADVs, not global ADV IDs (MB-017: using the IDs shifted the design by one index)
  +                        Vector< sint >          tLSAdvIds = mPerformerManager->mGENPerformer( 0 )
  +                                ->get_design_adv_positions( tReinitPerformer->get_adv_field_name() );
  ```
  plus `Geometry_Engine::get_design_adv_positions()` (ID → position through `mFullADVIds`), and a re-cut on the
  redistanced design before FEM.
- **Verification:** `X29_inplace_f1` (same deck, redistance every iteration): `wrote 4961 of 4961`, stays attached,
  strain energy 77.18 at iteration 30 (scheduled-restart arm 77.16, no upkeep 78.42).
- **Regression test:** `projects/GEN/test/ut_GEN_Field.cpp`, "Design ADV positions" [MB-017]: two named B-spline level sets; positions valid, unique, covering; IDs differ from positions (the test discriminates); a slice write at one design's positions sets only that design's field (partition of unity: 1 vs 0).

---

## MB-018 · GEN: stored-field MTK export dereferenced a null mesh
- **Status:** fixed, commit a0fe9502e on `fix/xiga-deck-bugs` (2026-10-10).
- **Symptom:** `Segmentation fault ... moris::mtk::Field::Field(Mesh_Pair const&, ...)` from
  `gen::Field_Discrete_Integration::get_mtk_field()` ← `Geometry_Engine::get_mtk_fields()` ←
  `Reinitialize_Performer::compute_sdf_from_interface()`, at the first redistance.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/queue/failed/050_x18sdf_smoke.sh` (L-bracket stress deck, `XIGA_SDF=1`;
  corner block and load pad are static user-defined geometries with `discretization_mesh_index -1`).
- **Evidence:** `runs/moris_bugs/MB-018/segfault_excerpt.txt`.
- **Root cause:** `Field_Discrete_Integration`'s ADV constructor sets `mMeshPair( nullptr, nullptr )` (upstream `// FIXME`);
  `Stored_Field` derives from it and never sets a pair, but inherited `get_mtk_field()` → `create_mtk_field( mMeshPair )`.
- **Fix / Code:**
  ```diff
  // projects/GEN/src/ADV/fields/cl_GEN_Stored_Field.hpp
  +      std::shared_ptr< mtk::Field > get_mtk_field() override
  +      {
  +          return nullptr;
  +      }
  ```
  (`Level_Set_Geometry::get_mtk_fields()` already returns no field for a null one.)
- **Verification:** `052_x18sdf_smoke_r2` runs 5 in-place redistances with no crash.

## MB-019 · GEN/WRK: image-seeded designs could not be redistanced (name ignored)
- **Status:** fixed, commit a0fe9502e (2026-10-10), with guards so the failure can no longer be silent.
- **Symptom:** `SDF redistance skipped this step: 0 level-set ADV ids vs 10404 coeffs (mesh changed)` every iteration; the
  start-of-run initialization re-fires every evaluation; the run otherwise proceeds as if no redistancing were configured.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/queue/done/051_x18sdf_smoke_r.sh` (design geometry is a
  `SIGNED_DISTANCE_IMAGE` with `name "Design"`).
- **Evidence:** `runs/moris_bugs/MB-019/skip_excerpt.txt`, `after_fix.txt`.
- **Root cause:** `fn_GEN_create_field.cpp` built `Image_Signed_Distance_Field` without the `name` parameter, and its
  constructor named the field after the image file; `get_design_adv_positions( "Design" )` found nothing. The two
  label lookups in `cl_WRK_Reinitialize_Performer.cpp` additionally fell back silently (interface mode: mesh pair 0 /
  index 0, the June sweep's "layer 3") or dereferenced the end iterator (field_remap).
- **Fix / Code:**
  ```diff
  // projects/GEN/src/ADV/fields/cl_GEN_Image_Signed_Distance_Field.hpp
  -                bool                  aSdfInterpolate )
  -                : Field_Analytic< 0 >( {}, aImageFileName )
  +                bool                  aSdfInterpolate,
  +                const std::string&    aName = "" )    // MB-019: the "name" parameter was ignored (field named after the image file)
  +                : Field_Analytic< 0 >( {}, aName.empty() ? aImageFileName : aName )
  ```
  plus `tName` passed in `fn_GEN_create_field.cpp`, and `MORIS_ERROR`s for a missing `adv_field` label (both lookups) and for
  a design with no ADV positions (workflow slice write).
- **Verification:** `052_x18sdf_smoke_r2`: `wrote 10404 of 10404` on all 5 redistances.
- **Regression test:** none yet.

## MB-020 · GEN: multilinear intersection NaN when every corner φ of the cell is zero
- **Status:** fixed, commit abb99ee66 on `fix/xiga-deck-bugs` (2026-10-10), with regression test. First written in worktree `moris-mb020` (branch `fix/mb020-bilinear-nan`).
- **Symptom:** abort in the XTK cut (`Node_Hierarchy_Interface::determine_intersected_edges_and_make_requests` → `Level_Set_Geometry::create_intersection_node` → `Intersection_Node_Bilinear`):
  `Intersection_Node_Bilinear::compute_intersection - Newton did not converge: Reference residual nan Current residual nan Edge coordinate 1.000000e+00`
  at `projects/GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:424`. Only with `use_multilinear_interpolation = true`; the same run with it off (`X40_cheese_o2_nomlin`) completes 151 evals.
- **Reproducer:** `runs/moris_bugs/MB-020/run_repro.sh <outdir> [MORISROOT]` (config `runs/moris_bugs/MB-020/repro_full.json` = X40_cheese_o2: MBB, order-2 B-spline + Lagrange, multilinear, swiss-cheese seed, native MMA, in-place SDF redistance every step, 13 its). Deterministic: aborts in iteration 12 during "Re-cut on redistanced design", ~80 s, ~0.6 GB.
- **Evidence:** `runs/moris_bugs/MB-020/` — `signature_excerpt.txt`, `diag_excerpt.txt` (instrumented build: QUAD9 cell x∈[0.325,0.35], y∈[0.15,0.175], corner φ = 0, 0, 0, 0 exactly, mid-side φ = ±5e-18, −0.075, −0.075, centre −0.075; threshold 0; scaling = inf), `after_fix.txt`, `fix.diff`.
- **Root cause:** `projects/GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:227` scales the corner values by `1/norm(tPhiBCNodes)` before any degeneracy check. With all four (eight in 3D) corner values exactly 0, `tPhiScaling = inf`, so `tPhiBCNodes = 0·inf = NaN` and the scaled threshold `0·inf = NaN`. The sign-change test (`NaN > 0` is false) and the coincident-parent test (`|NaN| < eps` is false) both pass through, the Newton residual is NaN for all three initial guesses, and the final guess `1.0` is kept by the `min/max` clamp. Such cells occur at order 2: the corner (knot) values of a quadratic B-spline are averages of neighbouring coefficients, so an alternating coefficient pattern at the ±bound (here ±0.15 after the SDF clamp) is 0 at every corner but nonzero inside. The multilinear model sees φ ≡ 0 on the cell.
- **Fix:** skip the scaling when the corner norm is zero. The existing checks then classify the edge: both parents on the zero isocontour → coincident-parent branch returns 0.0, and XTK discards edges whose parents are on the interface.
- **Code:** (trimmed from `git show abb99ee66`)
  ```diff
  --- a/projects/GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp
  +++ b/projects/GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp
  @@ -223,8 +223,12 @@ namespace moris::gen
  -            // Scale element level set field such that norm equals 1.0
  -            real tPhiScaling = 1.0 / norm( tPhiBCNodes );
  +            // Scale element level set field such that norm equals 1.0. If every corner value is zero, the multilinear
  +            // field vanishes on the whole cell [...]: 1/0 would turn every scaled value into NaN (MB-020), so
  +            // leave the values unscaled; the checks below then classify the edge (no sign change / coincident parents).
  +            real tPhiNorm    = norm( tPhiBCNodes );
  +            real tPhiScaling = tPhiNorm > MORIS_REAL_MIN ? 1.0 / tPhiNorm : 1.0;
               tPhiBCNodes      = tPhiScaling * tPhiBCNodes;
               tIsocontourThreshold *= tPhiScaling;
  ```
- **Verification:** the reproducer against the worktree build runs all 13 iterations to rc 0 (`runs/moris_bugs/MB-020/repro_fix`), passing iteration 12 with no MB-001 "every vertex on an interface" error.
- **Regression test:** `Bilinear Intersections with Zero Corner Values` (tag `[MB-020]`) in `projects/GEN/test/ut_GEN_Intersection_Level_Set.cpp`: one QUAD4 with an all-zero nodal field, multilinear on, threshold 0; queues every edge and diagonal through `Geometry_Engine::queue_intersection` and checks the local coordinate is finite and in [-1, 1]. Without the fix it fails with the MB-020 signature (unexpected exception, `Reference residual nan`); with it, it passes (36 assertions), and the full `GEN-test.exe` passes (26 cases). Logs: `runs/moris_bugs/MB-020/unit_test/`.
- **Workaround:** `use_multilinear_interpolation = false` at order 2 (X40_cheese_o2_nomlin ran clean).

---

## Entry template
```
## MB-NNN · <module>: <one-line summary>
- **Status:** open | fix-written (<branch>, <date>) | fixed (<commit>) | merged | wontfix (<why>)
- **Symptom:** what the user sees; exact error text or log signature
- **Reproducer:** durable path (deck/run dir/unit test) + the minimal config that triggers it
- **Evidence:** `runs/moris_bugs/MB-NNN/`: pathology figure, log excerpt, raw artifacts (copied, so cleanups can't delete them)
- **Root cause:** file:line and the data condition (or "not investigated")
- **Fix:** what changed (commit), and why that is the right place
- **Code:** the offending lines (open) or a trimmed `diff` of the fix hunk (fixed), with file:line / commit
- **Regression test:** test file + tag
- **Workaround:** deck-level avoidance, if any
```

## MB-021 · OPT: in-place SDF checkpoints hold the raw trial point, not the redistanced design
- **Status:** fixed, commit 4ed1e2393 on `fix/xiga-deck-bugs` (2026-10-10), with regression test. Written in worktree `moris-mb021`
  (branch `fix/mb021-advs-iter`). Verified end to end: `runs/studies/lunch_and_learn_xiga/X82_mb021_check` (X9r_rcp, `max_iter` 4,
  worktree binary) has 936/950/1032/1091 coefficients at |s| = 0.15 in Iter_1..4 (X9r_rcp: 0 in each), no stray Iter_5, and the
  objective trace of its 5 evaluations equals X9r_rcp's first 5 to all 15 printed digits. Iter_0 holds the seed in both runs: the
  in-place redistance first fires in evaluation 1. Main `build_opt` rebuild queued as XIGA job `1950_build_mb021.sh`.
- **Evidence:** `runs/moris_bugs/MB-021/` — `checkpoint_counts.txt` (per-file clamp counts before/after + objective traces),
  `ut_prefix.log` / `ut_postfix.log` (the [MB-021] test on the unfixed and fixed source), `mb021_count.py`.
- **Symptom:** in `reinit_in_place` runs, `ADV_Alg_0_Iter_0..N-1.hdf5` contain no coefficient at the SDF clamp ±0.15; only the last file does
  (X9r_rcp: 0 of 4961 at |s| = 0.15 in Iter_0/1/50/100/150, 1194 in Iter_151; X81_o2_cellorder on 3bc5a1ca3: 0 in Iter_1..150, 725 in
  Iter_151). The redistanced φ = R(ψ) clamps at ±0.15, so the intermediate files are the raw GCMMA trial points ψ.
- **Root cause (confirmed by reading the code):** `write_advs_to_file` names the file from `gLogger.get_opt_iteration()`. The criteria call
  advances it: `Algorithm::compute_design_criteria` → `communicate_running_status()` → `MORIS_LOG_ITERATION()` for
  `Task::compute_criteria_forward_analysis` (`projects/OPT/src/cl_OPT_Algorithm.cpp:390-394`). So in `Algorithm_MMA::func` the ψ write
  went to Iter_k, the φ rewrite to Iter_{k+1}, and the next evaluation's ψ write overwrote Iter_{k+1}. Only the last φ survived, in a
  stray Iter_N past the last evaluation. The TPL GCMMA wrapper writes once, before the criteria call, so it is not affected.
- **Fix:** an overload `Algorithm::write_advs_to_file( aADVs, aOptIter )`; `func` captures the index before `compute_design_criteria`
  and names both writes from it. The one-argument signature delegates with the logger's iteration (other callers unchanged).
- **Code** (commit 4ed1e2393):
  ```diff
  // projects/OPT/src/cl_OPT_Algorithm_MMA.cpp — Algorithm_MMA::func()
  -        this->write_advs_to_file( tADVs );
  +        // Capture the index first: compute_design_criteria() advances the optimization iteration (MB-021)
  +        const uint tOptIter = gLogger.get_opt_iteration();
  +        this->write_advs_to_file( tADVs, tOptIter );
  
           this->compute_design_criteria( tADVs );
  ...
               // checkpoints (and restarts) hold phi = R(psi) rather than the raw trial point psi (written above).
  -            this->write_advs_to_file( tADVs );
  +            this->write_advs_to_file( tADVs, tOptIter );

  // projects/OPT/src/cl_OPT_Algorithm.cpp
       Algorithm::write_advs_to_file( const Vector< real >& aADVs )
       {
  -        // Get iteration from global clock
  -        uint tOptIter = gLogger.get_opt_iteration();
  +        // Name the file from the iteration of the global clock
  +        this->write_advs_to_file( aADVs, gLogger.get_opt_iteration() );
  +    }
  +
  +    void
  +    Algorithm::write_advs_to_file( const Vector< real >& aADVs, uint aOptIter )
  +    {
           // Create file name
           std::string tRestartFileName =
  -                "ADV_Alg_" + std::to_string( mCurrentOptAlgInd ) + "_Iter_" + std::to_string( tOptIter ) + ".hdf5";
  +                "ADV_Alg_" + std::to_string( mCurrentOptAlgInd ) + "_Iter_" + std::to_string( aOptIter ) + ".hdf5";
  ```
- **Regression test:** `OPT/test/opt_test.cpp` "Native GCMMA in-place redistance checkpoints hold the redistanced design"
  `[optimization][MB-021]`: native MMA (3 its) with the reinit hook on the Rosenbrock + far-field problem of the MB-013 test; each
  Iter_k must equal evaluation k's redistanced ADVs and no Iter_N may exist. Unfixed source: 8 of 26 assertions fail; fixed: all pass.
- **Consequences (before the fix):** warm restarts from Iter_k < N started from ψ (one optimizer step off the redistanced design);
  per-iteration |∇φ| measurements (`sdf_quality.py`) and coefficient renders from intermediate checkpoints showed ψ, not φ. The final
  design was correct. After the fix the last evaluation's φ is in Iter_{N-1}, not Iter_N: tooling that picks the newest file
  (`sort -V | tail -1`) is unaffected, tooling that hard-codes `Iter_<max_iter+1>` must change.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/X81_o2_cellorder/results/checkpoints/` (count |s| = 0.15 per file with
  `runs/moris_bugs/MB-021/mb021_count.py`).
- **Workaround (runs on builds before 4ed1e2393):** use only the final checkpoint as a φ restart; treat intermediate checkpoints as ψ.


## MB-022 · GEN: Nested intersection recursion drops the incoming sensitivity factor
- **Status:** fix-written (`fix/mb022-025-audit @307236d81`; found by the implementation audit 2026-10-10 at `a3b602862`). Not yet in the main `build_opt`.
- **Symptom:** Scaling the incoming chain factor from I to 2I leaves parent contributions unscaled: nodes 21 and 22 have maximum absolute errors 0.410097050636513 and 0.41009705063651286; nodes 9 through 20 pass.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/implementation_audit_20261010/gen_chain_audit.cpp`, Catch test `AUDIT nested intersection chain factor`. Build and run instructions in that directory's `README.md`.
- **Evidence:** `runs/moris_bugs/MB-022/`; native log `gen_chain_test.log`; defect intake D-028.
- **Root cause:** `projects/GEN/src/PDV/cl_GEN_Intersection_Node_Level_Set.cpp:85`. The direct field term includes aSensitivityFactor, but each recursive parent call constructs only its local Jacobian and discards the incoming factor. For incoming A and parent Jacobian J, recursion must receive A*J. A top-level call with I masks the issue for a single layer of moving parents; deeper nesting loses ancestor factors.
- **Fix:** parents receive `aSensitivityFactor * J` (incoming factor times the local parent Jacobian) in both parent branches.
- **Code** (fix, `git show 307236d81`, trimmed):
  ```diff
  -            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_first_parent() );
  +            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_first_parent() );
               tFirstParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
  @@ -95,3 +99,3 @@ namespace moris::gen
  -            Matrix< DDRMat > tSensitivityFactor = 0.5 * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_second_parent() );
  +            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_second_parent() );
               tSecondParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
  ```
- **Regression test:** `GEN-test.exe "[linear intersection]"` (`ut_GEN_Intersection_Level_Set.cpp`) checks D(A) = A·D(I) for every derived node with the non-commuting A = [[1, 2], [-3, 0.5]]; on `a3b602862` it fails for nodes 21 and 22 (2 of 510 assertions), with the fix it passes. The three-level fixed-branch coordinate FD test proposed by the audit is not written yet.
- **Workaround / scope:** Until the fix is in the production build: The manual recipe cuts ADV-free pads first and the moving geometry afterward, avoiding the affected moving-parent recursion. Restrict nesting until fixed; this audit does not establish a failing recipe gradient.


## MB-023 · OPT: Native GCMMA accepts NaN criteria in the conservativeness inner loop
- **Status:** fix-written (`fix/mb022-025-audit @307236d81`; found by the implementation audit 2026-10-10 at `a3b602862`). Not yet in the main `build_opt`.
- **Symptom:** A fault-injection test with max_inner_its=3 returns NaN at the third criteria evaluation. MMA requests a gradient at that failed point and returns objective=nan at x=0.48074069641608658.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/implementation_audit_20261010/opt_invalid_audit.cpp`, Catch test `AUDIT native MMA rejects invalid inner trial`. Build and run instructions in that directory's `README.md`.
- **Evidence:** `runs/moris_bugs/MB-023/`; native log `opt_invalid_test.log`; defect intake D-029.
- **Root cause:** `projects/OPT/src/cl_OPT_Algorithm_MMA.cpp:1037`. The outer trial is checked by criteria_finite(), but the trial inside the conservativeness loop is not. Every comparison of the form fapp+eps < NaN is false, leaving conserv=1; the failed point is committed and grad() is called.
- **Fix:** every conservatism inner trial goes through the same NaN-guard backtrack as the first trial (`Algorithm_MMA::backtrack_non_finite`: halve the move limit, re-solve, re-evaluate, at most 15 times). A trial that stays non-finite rejects the iteration (MB-024 path); it is never tested for conservatism, committed, or differentiated.
- **Code** (fix, `git show 307236d81`, trimmed):
  ```diff
  -                    conserv = 1;
  -                    if ( ( f0app + epsimin ) < f0valnew ) conserv = 0;
  +                // a failed inner trial must not reach the conservatism test: NaN comparisons are false, so it
  +                // would pass as conservative and be committed and differentiated (MB-023)
  +                tTrialStatus = this->backtrack_non_finite( iter, nfunc );
  +                if ( tTrialStatus == Trial_Status::restart )
  +                {
  +                    return istop;
  +                }
  +                if ( tTrialStatus == Trial_Status::failed )
  +                {
  +                    break;
  +                }
  +                conserv = 1;
  +                if ( ( f0app + epsimin ) < f0valnew ) conserv = 0;
  ```
- **Regression test:** `OPT-test.exe "[MB-023]"` (same quadratic, `max_inner_its = 3`, NaN injected at evaluation 3): asserts no gradient at a failed evaluation and a finite, self-consistent returned objective. Fails on `a3b602862`, passes with the fix.
- **Workaround / scope:** max_inner_its=0 avoids this particular branch, including in the manual recipe. This is a test of failure handling using injected NaN, not proof that the recipe naturally generates this failure.


## MB-024 · OPT: Exhausted MMA backtracking leaves Problem on the rejected trial
- **Status:** fix-written (`fix/mb022-025-audit @307236d81`; found by the implementation audit 2026-10-10 at `a3b602862`). Not yet in the main `build_opt`.
- **Symptom:** With a finite initial x=0.5 and injected NaN for every subsequent trial, all 15 backtracks fail. Despite logging holding previous design, solve returns Problem objective=nan and x=0.49999084573096514.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/implementation_audit_20261010/opt_invalid_audit.cpp`, Catch test `AUDIT native MMA restores accepted state after failed backtracking`. Build and run instructions in that directory's `README.md`.
- **Evidence:** `runs/moris_bugs/MB-024/`; native log `opt_backtrack_test.log`; defect intake D-030.
- **Root cause:** `projects/OPT/src/cl_OPT_Algorithm_MMA.cpp:1007`. continue preserves private xval/f0val but does not restore mProblem, which every func evaluation overwrites through Problem::compute_design_criteria (lines 218-236). Manager getters return this polluted Problem state. If the iteration budget ends here, the returned result is the rejected trial.
- **Fix:** after exhausted backtracking, `Algorithm_MMA::restore_accepted_state` re-evaluates the held design `xval` (one extra evaluation) so the shared Problem state (ADVs, geometry, analysis, criteria) matches it before `continue`. In reinit mode the restore evaluates R(xval), so the Problem holds that redistance with consistent criteria. If the restore evaluation is itself non-finite the solve stops (returns 1) instead of iterating on an unevaluable design.
- **Code** (fix, `git show 307236d81`, trimmed):
  ```diff
  -            if ( !this->criteria_finite() )
  -            {
  -                // could not recover: hold the previous (finite) design and continue
  -                MORIS_LOG_INFO( "Native MMA iter %d: criteria still non-finite after %d backtracks; "
  -                                "holding previous design and continuing.",
  -                        (int)iter, tMaxBacktrack );
  -                continue;
  -            }
  +            if ( tTrialStatus == Trial_Status::failed )
  +            {
  +                // could not recover: hold the previous (finite) design. The failed evaluations overwrote the
  +                // shared Problem state, so re-evaluate xval before continuing; the Manager reads that state (MB-024)
  +                if ( !this->restore_accepted_state( iter, nfunc ) )
  +                {
  +                    return mProblem->restart_optimization() ? istop : 1;
  +                }
  +                MORIS_LOG_INFO( "Native MMA iter %d: holding previous design and continuing.", (int)iter );
  +                continue;
  ```
- **Regression test:** `OPT-test.exe "[MB-024]"`: NaN at every design other than the initial x = 0.5 (a design-dependent failure, unlike the audit harness, which also failed the accepted design and so could not be restored). Asserts 18 evaluations (initial, trial, 15 backtracks, restore), no gradient on a failed evaluation, and the Problem returning x = 0.5 with objective 0.25. Fails on `a3b602862`, passes with the fix.
- **Workaround / scope:** Treat exhausted-backtracking output as a failed solve and recover from a known accepted state. max_inner_its=0 does not avoid this branch; the default recipe can reach it if all retries fail.


## MB-025 · OPT: SWEEP single finite-difference modes silently skip all perturbations
- **Status:** fix-written (`fix/mb022-025-audit @307236d81`; found by the implementation audit 2026-10-10 at `a3b602862`). Not yet in the main `build_opt`.
- **Symptom:** forward, backward, and central each perform only one baseline criteria evaluation and produce HDF5 files with no gradient datasets. The all control performs three evaluations and writes analytic plus all three FD gradients.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/implementation_audit_20261010/opt_invalid_audit.cpp`, Catch test `AUDIT SWEEP individual FD modes perform perturbations`. Build and run instructions in that directory's `README.md`.
- **Evidence:** `runs/moris_bugs/MB-025/`; native log `opt_sweep_test.log`; defect intake D-031.
- **Root cause:** `projects/OPT/src/cl_OPT_Algorithm_Sweep.cpp:291`. The parameter API documents all three individual modes, but sweep_solve executes finite differences only under mFiniteDifferenceType == all. No alternative dispatch handles individual modes. The manual already notes this defect; the audit independently reproduces it.
- **Fix:** `sweep_solve` dispatches `forward`/`backward`/`central` to `SA_Type::forward/backward/central` and saves `fd_<mode>` datasets; the analytic gradient is always saved as the reference; unknown `finite_difference_type` strings abort in the constructor; `all` writes `fd_central` only for the gradient types requested. `none`/`all` output is otherwise unchanged (no in-repo deck uses an individual mode).
- **Code** (fix, `git show 307236d81`, trimmed):
  ```diff
  +                    else
  +                    {
  +                        SA_Type tType = SA_Type::central;
  +                        if ( mFiniteDifferenceType == "forward" )
  +                        {
  +                            tType = SA_Type::forward;
  +                        }
  +                        else if ( mFiniteDifferenceType == "backward" )
  +                        {
  +                            tType = SA_Type::backward;
  +                        }
  +
  +                        this->set_sensitivity_analysis_type( tType );
  +                        this->compute_design_criteria_gradients( tEvaluationPoint );
  +                        this->evaluate_objective_gradients( tEvaluationName + " fd_" + mFiniteDifferenceType );
  +                        this->evaluate_constraint_gradients( tEvaluationName + " fd_" + mFiniteDifferenceType );
  +                    }
  ```
- **Regression test:** `OPT-test.exe "[MB-025]"`: for none/forward/backward/central/all, asserts the criteria-evaluation count (1/2/2/3/3) and that every requested HDF5 dataset exists and matches the analytic gradient within 1.1e-2. Fails on `a3b602862`, passes with the fix.
- **Workaround / scope:** Use finite_difference_type=all and inspect the central-difference datasets. The manual reported sweeps use this workaround, so their results are not invalidated by this finding.
