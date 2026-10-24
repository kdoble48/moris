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
| MB-001 | merged (#25, 82b186bf7; fixed 992639917) | XTK | Segfault building the enriched integration mesh when a cell's vertices all lie on one geometry's interface → now a clear error | lunch_and_learn_xiga B13 |
| MB-002 | merged (#25, 82b186bf7; fixed 22c3c2799; residual = MB-011) | FEM/INT | Double-sided FD geometry sensitivity uses the leader's parametric coords for the follower → dR/dp ~1e8 on interfaces along cell edges | lunch_and_learn_xiga B12 |
| MB-003 | open | GEN/FEM | Static user-defined geometry with `discretization_mesh_index = -2` (analytic) makes design gradients NaN, silently | moment_fitting (2026-07-16) |
| MB-004 | open | FEM | Perimeter (interface VOLUME IQI) adjoint off by 12–150% on a few small-gradient ADVs | lunch_and_learn_xiga B11 |
| MB-005 | open (root-caused; native-GCMMA workaround) | OPT/HMR | Adaptive refinement during optimization: ADV count drifts between remeshes → abort in `compute_design_criteria_gradients` | seeding_impact (2026-07-13) |
| MB-006 | open (suspect) | FEM/INT | Nonconformal branch of `select_dRdp_FD_geometry_double` also copies the leader ξ into the follower row | review of MB-002 |
| MB-007 | open | FEM/VIS | Requested IQI evaluated on zero sets writes zeros to Exodus with no warning | Hertz contact (2026-07) |
| MB-008 | merged (#63, d355343dd) | FEM/INT | `IQI_Max_Stress::compute_dQIdu` stops with "not implemented"; decks must switch the whole forward analysis to FD | moris_stress_demo (2026-10-06) |
| MB-009 | merged (#25, 82b186bf7; fixed a856e9c13) | GEN | Design factory rejects a level-set geometry listed before the fields it depends on ("does not recognize geometry type"); properties built on the wrong field | MB-001 regression test |
| MB-010 | open | XTK | Phase boundary along background grid lines creates no interface facets → phases silently decoupled | MB-001 verification |
| MB-011 | merged (#25, 82b186bf7; fixed 06515ef67) | FEM/INT | Cluster-measure FD sensitivities (Nitsche-interface SP) use a base residual from CM/SP caches of the last perturbed geometry → 0.5–15% gradient error next to interfaces | MB-002 verification |
| MB-012 | open | VIS | STANDARD VIS mesh writes non-finite nodal values (UX/UY/VM) at 11–15% of nodes, mostly void-side but some with φ ≥ 0 | X18 peak extraction |
| MB-013 | open (limitation) | WRK/OPT | Every-step in-place SDF redistancing edits the ADV vector between criteria and gradient evaluation → ADV-consistency `MORIS_ERROR` after O(100) iterations (native MMA only; GCMMA cannot use it) | sdf_reinit_sweep (2026-06) |
| MB-014 | open | WRK | Scheduled-restart SDF apply in `Workflow_HMR_XTK::initialize()` replaces the whole ADV vector with one field's coefficients → silently skipped ("does not match current ADV count") on composite designs (level set + PDVs) | sdf_reinit_sweep (2026-06), review 2026-10-10 |
| MB-015 | open | WRK | `sdf_reinit_mode = "field_remap"` throws `std::out_of_range` in `Reinitialize_Performer::perform()` | sdf_reinit_sweep (2026-06) |
| MB-017 | merged (#34, 4900ba1a8) | GEN/WRK | Every-step in-place SDF redistance wrote coefficients at their global ADV IDs instead of ADV-vector positions → design shifted by one coefficient per redistance ("wrote 4960 of 4961"), MBB cut off its load pad | warm-restart work 2026-10-10 (D-015) |
| MB-018 | merged (#34, 4900ba1a8) | GEN | `Stored_Field::get_mtk_field()` built an MTK field on a null mesh pair → segfault in `mtk::Field::Field` whenever SDF reinit (or the remeshing restart) asked GEN for its MTK fields with a static stored geometry | L-bracket SDF port 2026-10-10 |
| MB-019 | merged (#34, 4900ba1a8) | GEN/WRK | `SIGNED_DISTANCE_IMAGE` field ignored the `name` parameter (named after the image file), so SDF reinit found no design: redistance skipped silently every iteration | L-bracket SDF port 2026-10-10 |
| MB-016 | merged (#40, 97c863cb7; fixed 2ae1dfb3f) | XTK | Per-geometry cut-cell selection uses only background corner nodes; an analytic field evaluated at Box-cut child vertices splits a shared edge on one side only → nonconforming IG facets → SIGSEGV in `construct_subphase_neighborhood` (opt) | gaussian_splat_ls G-N24 (D-012) |
| MB-020 | merged (#25, 82b186bf7; fixed abb99ee66) | GEN | Multilinear intersection divides by the norm of the four corner φ values; a cell whose corners are all exactly 0 (order-2 field nonzero only at mid-side/interior nodes) gives NaN residuals → abort "Newton did not converge: Reference residual nan" | lunch_and_learn_xiga X40_cheese_o2 (D-020) |
| MB-021 | merged (#26, 33be1e5c6; fixed 4ed1e2393) | OPT | Native GCMMA in-place SDF mode: the "redistanced checkpoint" rewrite (3ebdd8050) lands in the NEXT iteration's file and is then overwritten by the raw trial point, so Iter_0..N-1 hold ψ, never φ = R(ψ); only the last file holds φ | lunch_and_learn_xiga theory manual audit (X9r_rcp, X81_o2_cellorder) |
| MB-022 | merged (#29, de7990fc2; fixed ef9f2a0c0, dc1935cad) | GEN | Nested intersection recursion drops the incoming sensitivity factor | implementation audit 2026-10-10 (D-028) |
| MB-023 | merged (#29, de7990fc2; fixed ef9f2a0c0, dc1935cad) | OPT | Native GCMMA accepts NaN criteria in the conservativeness inner loop | implementation audit 2026-10-10 (D-029) |
| MB-024 | merged (#29, de7990fc2; fixed ef9f2a0c0, dc1935cad) | OPT | Exhausted MMA backtracking leaves Problem on the rejected trial | implementation audit 2026-10-10 (D-030) |
| MB-025 | merged (#29, de7990fc2; fixed ef9f2a0c0, dc1935cad) | OPT | SWEEP single finite-difference modes silently skip all perturbations | implementation audit 2026-10-10 (D-031) |
| MB-026 | merged (#29, de7990fc2; fixed ef9f2a0c0, dc1935cad) | OPT | FD perturbation sign toggles on every backward call: SWEEP fd_forward/fd_backward swapped for every second (point, epsilon) call | code review of fix/mb022-025-audit, 2026-10-10 |
| MB-027 | merged (#34, 4900ba1a8) | WRK/OPT | In-place SDF redistance R is not idempotent (each application raises the KS stress aggregate ~0.6 % and the volume ~0.075 %) and is applied to every inner-loop trial, so GCMMA conservativeness compares f(R(ψ)) with an approximation of f(ψ): trials converge to f(R(x_k)) > approximation, the inner loop hits its cap, and the accepted R(x_k) ratchets the constraint up | X18 v8 KS mimic h1 (D-056), perf/x18v8 |
| MB-028 | merged (#26, 33be1e5c6; was cited as MB-016) | OPT | Native in-tree GCMMA (`algorithm = "mma"`) wrote no `ADV_Alg_0_Iter_N.hdf5` restart files: no warm restart, no post-run tail-velocity check, silently | native-GCMMA validation 2026-10-10 |
| MB-029 | merged (#46, 52b489c33; fixed ee83435b7 on feature/gen-additive, cited as MB-017) | GEN/MRS | `Design_Variable` ID is a `char` from a static counter; `ADV_Manager::create_adv` aliases any non-constant variable whose ID repeats → >~250 analytic Design_Variables silently share ADVs | gaussian_splat_ls G-N60 (D-011) |
| MB-030 | merged (#47, 3a5950ed8; fixed 496ac24d2, d3f17589c on feature/gen-additive) | OPT | ROL wrapper: `al_initial_penalty` silently ignored (ROL "Use Default Initial Penalty Parameter" defaults true), and ROL's Type-G AL updates multipliers only when a capped inner solve CONVERGES, so with `subproblem_iteration_limit` the AL silently degenerates to a fixed quadratic penalty with lambda = 0 | gaussian_splat_ls task 34 (code reading; explains the 2026-07 rol_applicability "penalty keys inert" result) |
| MB-036 | open (was cited as MB-030; test tag still `[MB-030]` on develop) | GEN | Nested intersection sensitivities use the field gradient projected onto an intersection parent's own edge; wrong once that parent moves off its edge (nesting level ≥ 3): 38–64 % relative error vs FD | MB-022 FD regression test 2026-10-11 |
| MB-031 | merged (#40, 97c863cb7; fixed 2ae1dfb3f) | XTK/GEN | Cut-cell detection tests only background-corner signs, but an analytic (`discretization_mesh_index = -2`) field is evaluated exactly at subdivision child vertices: a sub-cell void (corners φ > 0, cell centre φ < 0) is invisible until one corner flips sign, then appears whole → objective and mass jump discontinuously (same root as MB-016) | gaussian_splat_ls task 34 FD jumpcheck (D-074) |
| MB-032 | open | SOL/NLA (+deck) | NLA's absolute `tResNorm > 1e20` guard aborts the linear adjoint solve (`solve_implicit_DqDs`) when any criterion's dQ/du is large, including a criterion the optimizer weights by zero (KS stress IQI, fixed shift 1, in compliance mode) | triage of D-077/D-107, 2026-10-12 |
| MB-033 | merged (#38, 2c99051e1; fixed 63c42e266) | WRK | A remesh requested right after a failed XTK decomposition refines on MDL fields (e.g. `FieldZienkiewiczZhu`) still living on the pre-remesh mesh → SIGSEGV (null parent) in `Refinement_Mini_Performer::perform_refinement_2` | triage of D-150 (ZZ showcase), 2026-10-12 |
| MB-034 | merged (#38, 2c99051e1; fixed f0fb28d47) | WRK | `minimum_refinement_level` (level, until-iteration) pairs are read with stride 1: from the second pair on, the wrong entries are used, silently (`0, 1, 30, 2, 100` gives level 0, not 2, for iterations 30-99) | triage of D-150, 2026-10-12 |
| MB-035 | open | MTK/GEN | B-spline design field with Lagrange order < B-spline order (p = 2 on Lagrange 1): the node-to-B-spline L2 projection is singular; above 25000 coefficients the mapper switches to Belos, which does not converge and GEN uses the garbage field silently (φ up to ±2373, cut area +1049 %) | sphere3d_precision (2026-10-13, D-158) |
| MB-038 | fixed (86ce7f9b2, PR #52 open) | MTK/HMR | GCC 13.3 miscompiles the base-object ctor of classes with a virtual `mtk::Mesh` base (16-aligned via armadillo `Matrix` members): aligned `movaps` stores at `this+off` while the ABI places the subobject at 8 mod 16 → SIGSEGV in `Integration_Mesh::Integration_Mesh()`; latent, exposed by any change that makes GCC call the out-of-line C2 | Vector move-semantics gate 2026-10-20 (XTK-test "XTK Cell Clusters") |
| MB-039 | fixed (4b6a87b31, PR #52 open) | GEN | `gen::ADV` keeps a raw `real*` into `ADV_Manager::mADVs`; any reassignment that replaces the buffer leaves dangling ADVs, and `Geometry_Engine` clears that vector before `distribute_advs()` repoints them (reads past `size()`, worked by luck) | Vector move-semantics gate 2026-10-20 (GEN-test Circle/Sphere/User-defined) |
| MB-037 | open | XTK/GEN | An interface just outside the 1e-12 value tolerance of an existing vertex (a corner, a subdivision vertex or another geometry's cut vertex) produces IG tets with volume ≤ 0: 46 % of near-zero corner patterns at δ = 1e-11 (min −1.9e-18), and two parallel planes 1e-11 apart. At δ = 1e-12 a two-plane layer is partly snapped (keeps 7.6e-14 of 1.06e-12) yet survives as a subphase. No distance-based snapping; `intersection_tolerance` is never read | XTK walk-through pathology sweep `[XTK_PATHOLOGY_HEX_PARALLEL]`, 2026-10-20 |
| MB-040 | merged (#53, a6c1905d7; fixed 132e55764) | XTK/GEN | Analytic or multilinear level sets, ≥ 2 geometries: two-locator derived nodes are labelled from their locators while the node hierarchy cuts edges by the field values; inconsistent cut patterns (3 cut edges forming a face, 5–6 cut edges) give undefined templates, an uninitialized `mNumCells` and an out-of-bounds IG cell write (heap corruption). 828 (analytic) / 808 (multilinear) of 104 976 two-geometry single-hex cases | XTK walk-through pathology sweep `[XTK_PATHOLOGY_HEX_TWO_GEOMETRIES]`, 2026-10-20 |
| MB-041 | resolved by #55 (0f191448f: no IG cell is dropped); PR #54 closed, VIS fix kept on branch (ce59f590f) | FEM/VIS | IG elements dropped from the residual (cluster relative volume < 1e-6) are also skipped by VIS output, so their nodal/elemental fields keep the initial 0: an exact uniform stress shows as 0 in 240 of 1032 cells of an immersed patch test | promoted from D-182 (MB-037 patch test), 2026-10-21 |
| MB-042 | open | WRK/HMR/EXA | EXA `Levelset_Boxbeam_Adaptive_Refinement_Restart-1-procs` fails on develop: the create run writes its reinit refinement pattern as `HMR_Background_Refinement_Iter_10.hdf5`, the restart deck opens `..._Iter_11.hdf5` (the older moris/ line wrote Iter_11) | promoted from D-151, 2026-10-21 |
| MB-043 | open | HMR | Parallel deadlock: `Database::perform_refinement` calls the collective `Background_Mesh_Base::perform_refinement` as `tFlag \|\| perform_refinement()`, so ranks whose rank-local `mHaveRefinedAtLeastOneElement` is true (set by `flag_element` / `flag_parent`) skip it while the others enter it. Manual flagging on some ranks only hangs forever. No production caller of the manual flag API today | found by the XTK multi-rank pathology tests, D-186, 2026-10-21 |
| MB-044 | open | XTK | `check_intersected_background_cell_levels` compares each rank's intersected cells with that rank's first one and combines with `all_land`: cut cells on different refinement levels on different ranks pass in parallel and are decomposed, while 1 rank refuses the same problem ("Decomposition Failed. Not all child meshes on the same level", the remesh signal) | XTK walk-through finding 6, D-187, 2026-10-21 |
| MB-045 | open (limitation, mechanism not isolated) | XTK/FEM | A level-set neck of near-zero width at a grid node (phi = +1.3e-6) stays connected through cut cells of area ~1e-11 and carries load: strain energy 54.78 vs 76.34 (+39 %) when the node value is -7.0e-6; the discontinuity stalls MMA (step collapse) | ls_smoothing study, helm_L2 seed 1 (D-188) |
| MB-046 | merged (#57, 0117bd3ad; verified, runs/studies/xtk_pathology/patch_test/VERIFICATION.md) | FEM/INT | 3D principal stress uses σyz instead of σzz in the σxy² term of the third invariant (`IQI_Stress`, `IQI_Max_Stress`, `IQI_Zienkiewicz_Zhu`): σ1 off by 1.2 % in a uniform-stress patch test, every cell, every output type | stress-output-path study 2026-10-21 (MB-037 patch test) |
| MB-047 | merged (#57, 0117bd3ad; verified, runs/studies/xtk_pathology/patch_test/VERIFICATION.md) | FEM/INT | One CM named on both the leader and the follower side of an IWG/IQI is a single shared object: its traction cache ignores the normal argument and its FI manager is overwritten by the follower's. Interface Nitsche residual inconsistent with its Jacobian: Newton converges linearly to a different solution (1.6 % in nodal σxx, 2x2x2 non-patch case), 1-step runs report a residual drop of 1e-5..3e-4 (this is D-183); the patch test cannot see it | stress-output-path study 2026-10-21 |
| MB-048 | merged (#57, 0117bd3ad; verified, runs/studies/xtk_pathology/patch_test/VERIFICATION.md) | FEM/INT | `JUMP_TRACTION` returns t_L(n) − t_F(−n) = (σ_L + σ_F)n, the sum, not the jump: 2|σn| in exact equilibrium (separate CMs); 0 identically with one shared CM (MB-047) | stress-output-path study 2026-10-21 |
| MB-049 | merged (#57, 0117bd3ad; verified, runs/studies/xtk_pathology/patch_test/VERIFICATION.md) | FEM/INT | `IQI_Max_Stress` bounds check inverted (`MORIS_ERROR( mIQITypeIndex > 2, ...)`): MAX_NORMAL/SHEAR/PRINCIPAL_STRESS abort for every valid index 0..2; only MAX_VON_MISES_STRESS is usable | stress-output-path study 2026-10-21 |
| MB-050 | merged (#57, 0117bd3ad), NOT exercised: no deck or test gives an SP follower CMs (suites unchanged) | FEM/INT | `Model_Initializer_Legacy::create_stabilization_parameters` passes an SP's follower CMs without the side argument (default LEADER): the follower slots stay null and the leader slots are overwritten. Only `SP_Spalart_Allmaras_Nitsche_Interface` reads follower CMs | code review during the MB-047 fix, 2026-10-21 |
| MB-051 | merged (#57, 0117bd3ad; verified, runs/studies/xtk_pathology/patch_test/VERIFICATION.md) | FEM/INT | 2D principal stresses use entry 2 of the standardized stress vector (σzz) instead of entry 5 (σxy) in the max-shear radius (`IQI_Stress`, `IQI_Max_Stress`, `IQI_Zienkiewicz_Zhu`): in-plane shear ignored in every 2D PRINCIPAL_STRESS output | code review while documenting MB-046, 2026-10-21 |
| MB-052 | open | SOL/DLA | Belos: setting `"Output Frequency"` aborts the linear solve with `std::bad_variant_access` (`Linear_Solver_Belos` reads the integer parameter as `real` when passing it to Belos) | found by the MB-037 patch-test linear-solver ladder, 2026-10-21 |
| MB-053 | merged (#60) | FEM/INT | `IQI_J_Integral` wrong in every use: W = σᵀε (2W), `gradx(1)(1,0)` = ∂ux/∂y used for ∂uy/∂x, `aWStar` multiplies only W in the criterion overload, normal hard-coded (1,0): uniform patch J = +1.98e-7 vs exact −4.92e-7 | IQI audit tier 3 (2026-10-21) |
| MB-054 | merged (#59) | FEM/INT | `IQI_Eigen_Value::compute_QI(real)` forces `mIQITypeIndex = 0` when `dof_quantity` has ≤ 1 entry (the EXA Column_Buckling declaration): every EIGEN_VALUE criterion returns λ₁, and the reset persists into VIS | IQI audit tier 2 (2026-10-21) |
| MB-055 | merged (#59) | FEM/INT | `IQI_Stabilization` forces `mIQITypeIndex = 0` when `dof_quantity` has ≤ 1 entry: a request for tauC (any SP value but the first) silently returns tauM, value and gradient | IQI audit tier 3 (2026-10-21) |
| MB-056 | merged (#60) | FEM/INT | `IQI_Power_Dissipation_Bulk` computes σᵀε with the fluid CM's tensor-shear Voigt strain, dropping the ×2 on shear (the NS IWG weights diag(1,1,2)), and adds ρ\|v\|²: exact Poiseuille gives 6 + 1.2ρ instead of the dissipation 12 (numerically confirmed 2026-10-22) | IQI audit tier 3 (2026-10-21) |
| MB-059 | merged (#60) | FEM/INT | `DRAG_COEFF` / `LIFT_COEFF` return the traction on the fluid (fluid outward normal), i.e. −C_D and −C_L: Poiseuille bottom wall drag −6 vs the wall force +6, lift +6 vs −6 | IQI audit suspect, confirmed 2026-10-22 |
| MB-060 | merged (#59) | FEM/INT | `STRONG_RESIDUAL_INCOMPRESSIBLE_NS` (and `SP_CROSSWIND_INCOMPRESSIBLE_NS`, same code) reset `vectorial_field_index` to 0 when `dof_quantity` has ≤ 1 entry: the x-momentum residual request (`"VX"`, index 1) returns div v (0 vs 12; `"VX,VY"` gives 12) | IQI audit suspect, confirmed 2026-10-22 |
| MB-061 | merged (#60) | FEM/INT | `CONTACT_PRESSURE_CURRENT` returns n·(P n) (1st Piola–Kirchhoff with the current normal), 1.10–1.19× the Cauchy pressure n·σn at F = diag(0.9,1.1) / general F; `CONTACT_PRESSURE_REFERENCE` uses the current normal and aborts with a nonlinear CM | IQI audit suspect, confirmed 2026-10-22 |
| MB-062 | merged (#61) | FEM/INT | analytic geometric sensitivities (PR #35, `is_analytical_sensitivity` true) of STRAIN_ENERGY and STRUC_LINEAR_BULK omit the spatial variation of the CM properties: dSE/dr 2.7% off with E = 1 + 2x (FD path 7e-9, constant E 1e-8), nothing guards against it | IQI audit suspect, confirmed 2026-10-22 |
| MB-063 | merged (#60) | FEM/INT | `HEAT_METHOD_PENALTY` returns NaN in release when the projection exponent −2 g s PHID / B exceeds 709.8 (exp overflows, inf/inf); the guard is a `MORIS_ASSERT`. With the deck rule g/B = ln199/(2·band) this is \|PHID\| > 134·band on the penalized side | IQI audit suspect, confirmed 2026-10-22 |
| MB-064 | merged (#59) | FEM/INT | `IQI_Dof` analytic dQIdu aborts for a scalar `dof_quantity`; scalar derivatives on the output path (and in `IQI_Jump_Dof`) read component −1 | IQI audit B-1/B-3, unit-tested 2026-10-22 |
| MB-065 | merged (#59) | FEM/INT | DOF IQI with `dof_quantity "UY"` (one member of a vector group) returns ∫UX and a zero dQIdu | IQI audit B-2, unit-tested 2026-10-22 |
| MB-066 | open | FEM/INT | `CM_Struc_Linear::eval_geometric_stiffness` in 2D takes `mFlux( 2 )` as σxy; for PLANE_STRAIN / AXISYMMETRIC the flux has 4 entries and entry 2 is σzz / σθθ (shear is entry 3): wrong geometric stiffness and buckling loads (same class as MB-051). EXA Column_Buckling decks use plane stress (unaffected) | IQI audit A-12, confirmed by code reading 2026-10-22 |
| MB-067 | open | FEM/INT | `STRESS_VECTOR(_CAUCHY)` in 2D: `compute_QI( real )` adds the standardized 6×1 vector into a QI sized 3×1 (plane stress) or 4×1 (plane strain) by `get_matrix_dim`: Armadillo size error in debug, unchecked `+=` in opt (`ARMA_NO_DEBUG`); AXISYMMETRIC aborts in `get_matrix_dim` | IQI audit A-1, confirmed by code reading 2026-10-22 |
| MB-068 | open | FEM/INT | `SHEAR_STRESS` (and `MAX_SHEAR_STRESS`) in 2D: the header says the index is ignored, but the code reads standardized entry index+3, so index 0 or 1 returns σyz/σxz = 0 and only index 2 gives σxy. Decks in use all pass index 2 | IQI audit A-2, confirmed by code reading 2026-10-22 |
| MB-069 | open | FEM/INT | `*_STRESS_CAUCHY` IQIs with a 2D nonlinear CM and `model_type = PLANE_STRAIN`: the nonlinear CM returns a 3-entry Cauchy vector (no σzz), and `IQI_Stress::get_stress_vector` treats 3 entries as plane stress (σzz = 0): von Mises and principal values wrong, NORMAL_STRESS index 2 = 0 | IQI audit A-3 (trigger corrected), confirmed by code reading 2026-10-22 |
| MB-070 | open | FEM/INT | `MAX_*_STRESS` and `MAX_DOF` integrand pow( q/ref − shift, p ) with an undocumented default shift 1.0: NaN for a non-integer p wherever q/ref < shift (MAX_DOF `sign = -1`: everywhere), and the default penalizes q below ref; no guard in release (`MORIS_ASSERT`), ref = 0 unguarded | IQI audit A-4 + B-14, confirmed by code reading 2026-10-22 |
| MB-071 | open | FEM/INT | IQI `normalization` has no effect: `mReferenceValue` is stored but never read, so a numeric reference is ignored; `"design"` sets the rank-local value to 1.0 at iteration 0 before `sum_all` (criterion = nproc), returns raw values from iteration 1, never scales gradients, and sets the reference on `mIQIs( requested index )` | IQI audit C-1, confirmed by code reading 2026-10-22 |
| MB-072 | open | FEM/INT | NODAL IQI evaluation sets no normal (and, on standard/double-sided clusters, no IG point): normal-dependent IQIs (TRACTION, JUMP_TRACTION, THERMAL_ENERGY_DIFFUSIVE_FLUX, …) as NODAL VIS fields use the last facet's normal or an empty one (out-of-bounds read in release); the guard is a debug `MORIS_ASSERT` | IQI audit C-2 (with A-8, B-13), confirmed by code reading 2026-10-22 |
| MB-073 | open | FEM/INT | `EIGEN_VECTOR` evaluated at a stale point outside `Element_Bulk::compute_QI` and the NODAL routines: the eigen FI manager is never positioned in the GLOBAL, ELEMENTAL_INT/AVG, field-output or dQIdu loops, so those values are the eigenvector at one leftover point | IQI audit C-3, confirmed by code reading 2026-10-22 |
| MB-074 | merged (#61) | FEM/INT | with `is_analytical_sensitivity = true` the element never reads the deck's `finite_difference_scheme` / perturbation for sensitivities: material pdvs, the FD geometry fallback and cluster measures run POINT_1_FORWARD with the default step | found while fixing MB-062, 2026-10-22 |
| MB-075 | merged (#71, 30c2cc00a) | XTK/HMR | a cut coarse cell next to an uncut finer cell: from the fine side HMR reports no transition and the fine cell's own facet, so `construct_subphase_neighborhood` reads a null child-facet list → SIGSEGV (release) / "Null ptr on facet that should have children facets" (debug), 1/2/4 ranks | connectivity robustness tests, 2026-10-21 |
| MB-078 | merged (#66, warning only) | FEM/INT | a set name repeated in an IWG's (or IQI's) `mesh_set_names` adds the same IWG to that FEM set again (`add_IWG` is a plain push_back) → assembled once per listing, silently: the EXA SIMP hole-seeding deck's `"SideSet_4_n_p1,SideSet_4_n_p1"` applies its traction twice (SE ×4) | serial/parallel study (D-193 follow-up), 2026-10-24 |
| MB-077 | merged (#70, 2b35a72b1) | OPT | `Problem::read_restart_file()` reorders restart ADVs by IjklIds when the file was written on a different rank count, then overwrites the result with the file-order vectors (`mADVs = tRestartADVs`, same for the bounds) → the design (and bounds) load scrambled on any other rank count, silently | serial/parallel study, same-design restarts, 2026-10-24 |
| MB-076 | fix-written (PR #64) | FEM/VIS | nonconformal output IQIs (GAP, CONTACT_PRESSURE, traction) written as ~0 in Exodus: VIS evaluates them at the frozen reference-configuration leader/follower pairing, never remapped to the deformed follower | nonconformal contact VIS (2026-07), filed in the branch cleanup 2026-10-23 |
| MB-079 | open | FEM/INT | stress IQIs (`IQI_Stress`, `IQI_Max_Stress`, `IQI_Zienkiewicz_Zhu`) check `vectorial_field_index` only when the first leader dof group has > 1 member: with a single-member first group (e.g. `"TEMP;UX,UY"`) and no index, NORMAL_STRESS reads standardized entry (uint)−1 (out of bounds in release), SHEAR_STRESS silently returns σzz, PRINCIPAL_STRESS aborts with a misleading "Only 1st and 2nd principal stresses known" error. 0 EXA decks meet the trigger | code review of PR #63, 2026-10-24 |
| MB-080 | merged (#69, c455c797b) | GEN | debug builds probe every user-defined field once with a 3-component point `{ 0, 0, 0 }` (`User_Defined_Field::validate_user_defined_functions`); a 2D field that loops over the point's size and indexes a 2-entry centre throws `vector::_M_range_check` (MB-031 2D XTK tests fail in debug on develop) | PR #67 debug check, 2026-10-24 |
| MB-081 | merged (#70, 2b35a72b1) | OPT | `Problem::read_restart_file` opens the restart file with write permission and no parallel extension (`open_hdf5_file( mRestartFile, false )`): any parallel restart aborts in debug on the HDF5_Tools.hpp:107 assert | PR #67 debug check, 2026-10-24 |
| MB-082 | merged (#68, 78e70a57d) | XTK | after #67 (merge before ids), a merged enriched IP vertex is shared by several unzipped cells; `sort_unzipped_vertices_into_owned_and_not_owned` lists and numbers it once per cell: debug stops on `Assertion tVertexTracker( tVertexIndex ) failed`, release wastes ids and keeps the last of k unchecked answers | review of #67, 2026-10-24 |
| MB-057 | merged (#60) | FEM/INT | `IQI_L2_Error_Analytic` field reference: the field FI is fetched with the dof type (error ≡ 0 in output), and the criterion overload has no field branch (null property dereference) | IQI audit (= studies/iqi_audit FINDINGS B-7), 2026-10-21 |
| MB-058 | open | FEM/INT/VIS | `IQI_Strain_Energy::compute_QI(Matrix)` ignores `Thickness`: the Exodus GLOBAL strain energy with Thickness 2 equals the thickness-free value (6.2857e-7 vs 1.2571e-6); the criterion path applies it | IQI audit tier 2 (2026-10-21) |

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
- **Status:** merged (#63, squash `d355343dd` on develop, 2026-10-24). The PR also moved the stress measures into
  `CORE/fn_FEM_Stress_Measures` (shared by `IQI_Stress`, `IQI_Max_Stress`, `IQI_Zienkiewicz_Zhu`, `IWG_Struc_Stress`), added the
  analytic geometric dQI/dp with the dC/dx term and the MB-062 guard, and the tests `UT_FEM_Stress_Measures`,
  `UT_FEM_IWG_Struc_Stress`, `Analytic_Geo_Sens_Varying_IQI_Max_Stress`; full INT-test.exe passed (310 cases). Earlier:
  fix-written, PR #63 (branch `feature/max-stress-ks-dqidu`, worktree `moris-pr-max-stress-ks`, 2026-10-23):
  analytic von Mises dQI/du plus a `KSParameter` (KS) integrand option, ported from `feature/analytic-geo-sens-ks`;
  `UT_FEM_IQI_Max_Stress` checks dQI/du against FD; INT/GEN/OPT/XTK suites pass. EXA sensitivity ctests not run.
  Previously: open (limitation; it fails loudly, so it is kept here only because it forces whole-model FD).
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
- **Status:** merged (#40, landed on develop as 97c863cb7, squash, 2026-10-20; fixed 2ae1dfb3f, 2026-10-13): same commit as MB-031. Fix (a) is implemented as `add_background_cells_intersected_at_child_vertices`, run before each node-hierarchy geometry pass. The new MB-031 subdivision-vertex test also covers the 3D form, where a negative hex face centre is split by only one neighbour. Fix (b), promoting the assert to an error, was not done.
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
- **Regression test:** `UT_XTK_Interior_Cut_Detection.cpp` `[MB-016_2D]` (line y = 0.4 + disc on the shared cut vertex: SIGSEGV on develop, passes with the fix) and `[MB-016_3D]` (ball on a shared face centre). See MB-031.
- **Verified (3D):** `runs/moris_bugs/MB-031/repro3d/`, the `AM_crash_Iter13` design pads-first, 1-eval SWEEP version of the `3d/probes/mb016_old` deck. It SIGSEGV'd before; on `local/mb031-repro` it runs with rc 0 (85 s, 0.8 GB, obj 0.6672469; 10 interior-only cells, 1 cell added at child vertices for geometry 1). See `stdout.txt` there.
- **Workaround:** `use_multilinear_interpolation = true` on the analytic design geometry (tested: the crash design runs; changes strain energy +11.7% at a sliver-rich design → FD-check gradients before use). σ_min ≥ 2h only lowers the frequency.
- **Workaround 2 (geometry order, 2026-10-12, gaussian_splat_ls task 37, 3D):** list the analytic design geometry FIRST (geometry 0) and the stored/linear fixed geometry (pads) second, remapping the phase table (2 geometries: index = 2 b(g0) + b(g1); pads-first `0;1;2;2` -> design-first `0;2;1;2`, same phases and set names). Consistent with the root cause: only geometry k > 0 is evaluated at geometry-0 child vertices, and a stored (linear) field cannot split a shared edge one-sided. Tested: the 3D crash design (`runs/studies/gaussian_splat_ls/3d/probes/AM_crash_Iter13.hdf5`) SIGSEGVs pads-first (`3d/probes/mb016_old`) and runs design-first (`3d/probes/mb016_df`); the seed's eval-1 criteria agree to 1.5e-5 between orders (`3d/probes/DF` vs `S3D_N120w25`). Does not address MB-031 (corner-sign detection of geometry 0 itself). **Partial only:** a later design (`3d/runs/G3D-N120df-s1`, GCMMA step 0.005 from `G3D-N120r-s1_Iter_2.hdf5`) crashes with the same signature design-first, so the order lowers but does not remove the crash (plausibly because the stored pad field is re-evaluated by its box min/max function at child vertices, which is not linear along a cut edge).
- **3D occurrences (task 37, 48x16x16, N=120 GS field + pads, s_min = 1.5 h):** 7 crashes in 154 outers of TPL GCMMA step 0.01/0.005 (D-080, D-081, D-083, D-084, D-088, D-093, D-094) and 9 crash records within the first ~13 outers of native MMA step 0.015 on the geometric ±40 box (D-085..D-096 AM3D rows; one deterministic replay), all pads-first. Defect rows link back here. The 1.5 h width floor does not stop it in 3D.

## MB-029 · GEN/MRS: `Design_Variable` char ID wraps → analytic ADVs silently aliased
- **Renumbered** 2026-10-11 from MB-017 (duplicate ID; MB-017 is the SDF redistance bug). Commit ee83435b7 on `feature/gen-additive` still says "MB-017".
- **Status:** merged (#46), landed on develop as 52b489c33 (2026-10-21). Fix-written on feature/gen-additive @ee83435b7 (2026-10-10); the PR branch `fix/mb029-design-variable-id` @4790db742 was a cherry-pick with the test tag renamed `[MB-029]`, gated 2026-10-20. In production after the next `work/develop-local` rebuild
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

## MB-028 · OPT: native GCMMA wrote no restart files
- **Renumbered** 2026-10-11 from MB-016 (duplicate ID; MB-016 is the XTK cut-cell bug). Commit f6d21bfb7, PR #26 and a comment in `cl_OPT_Algorithm_MMA.cpp` still say "MB-016".
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

## MB-030 · OPT: ROL augmented-Lagrangian keys silently inert (initial penalty ignored; multipliers never update under a capped inner solve)
- **Status:** merged (#47), landed on develop as 3a5950ed8 (2026-10-21). Fix-written on feature/gen-additive (2026-10-11: 496ac24d2 + d3f17589c); the PR branch `fix/mb030-rol-persistent-al` @6757c1e17, with 920e602d9, was gated 2026-10-20
- **Symptom:** in `al_mode = "rol"` (the only mode before the fix) changing `al_initial_penalty` does not change the run, and with
  `subproblem_iteration_limit` set (default 25) the AL table's multiplier never moves: the inner Lin-More solve exits
  `Iteration Limit Exceeded` every outer, so the run is a fixed quadratic penalty with lambda = 0 and the final design misses the
  constraint by O(1/rho). Matches the 2026-07 `runs/studies/rol_applicability/param/` sweep: 15/16 penalty / growth / radius
  combinations gave the identical design, and the first A/B ended infeasible (+0.035) or under-filled.
- **Reproducer:** any ROL deck with a constraint and `subproblem_iteration_limit` < the inner iterations needed to reach ROL's inner
  tolerance (min(1e-2 x initial projected gradient, ...)); e.g. `runs/studies/rol_applicability/param/Lbracket_sm_ROL.cpp`.
- **Evidence:** `runs/moris_bugs/MB-030/`: R1 (persistent mode, first version) outer log `R1_rolal_outers.log` showing the second
  defect below (rho 40 -> 2560), config `t34_R1-geo-N80.json`.
- **Root cause:** Trilinos 15.1.1 `ROL_TypeG_AugmentedLagrangianAlgorithm_Def.hpp:61-62` (`useDefaultInitPen_` true by default, so
  "Initial Penalty Parameter" is overwritten in `initialize`), and `:255` (`if ( algo->getState()->statusFlag == EXITSTATUS_CONVERGED )`
  gates the multiplier AND penalty update). The MORIS wrapper (`cl_OPT_Algorithm_ROL.cpp`, rol_solve) set neither flag.
  Second defect, in the first version of the new persistent mode: the penalty-growth test compared against a previous violation
  of 0 after a feasible outer, so rho grew 4x after every overshoot (40 -> 2560 in 11 outers) and the L-BFGS model became stiff
  (steps ~0.07 h; the GS-LS R1 run stalled at SE 5.17).
- **Fix:** 496ac24d2 sets "Use Default Initial Penalty Parameter" = false so the key takes effect, and adds `al_mode = "persistent"`:
  a PHR AL driven by MORIS with one secant and TR radius kept across outers and an unconditional multiplier update after each
  capped inner solve. d3f17589c grows rho only after an infeasible outer that did not cut the violation 4x.
- **Code:** (fix, `cl_OPT_Algorithm_ROL.cpp`)
```diff
+                tAugLag.set( "Use Default Initial Penalty Parameter", false );
                 tAugLag.set( "Initial Penalty Parameter", mALInitialPenalty );
...
-            if ( tViol > mALFeasibilityReduction * tPrevViolation && tViol > mConstraintTol )
+            if ( tPrevViolation > mConstraintTol && tViol > mConstraintTol
+                    && tViol > mALFeasibilityReduction * tPrevViolation )
```
- **Regression test:** `projects/OPT/test/opt_test.cpp` SECTIONs "ROL persistent augmented Lagrangian (binding constraint)" (the
  unconstrained minimum is infeasible; converges to (1,1), lambda 2) and "ROL forward-solve budget".
- **Workaround:** before the fix, pass a full `xml_file` with "Use Default Initial Penalty Parameter" false and an uncapped inner solve.

## MB-031 · XTK/GEN: corner-sign cut detection hides sub-cell features of an analytic field → objective/mass jump when a corner flips
- **Status:** merged (#40, landed on develop as 97c863cb7, squash, 2026-10-20; fixed 2ae1dfb3f on `fix/mb031-interior-cut-detection`, 2026-10-13, dev worktree `moris-mb031`). Cherry-picked to `feature/gen-additive` @d765f1a15 (not yet built there). Same root cause as MB-016 (the crash form; fixed by the same commit); this is the silent wrong-physics form, and it occurs with a single design geometry.
- **Symptom:** finite-difference SWEEP on a Gaussian-splat analytic level set (2D box beam 120x40, linear, `discretization_mesh_index = -2`): a 1e-4 ADV step (2.5e-6 length) on 4 of 60 ADVs (152, 240, 256, 304) changes the strain-energy objective by +6.43e-4 (0.11 %, identical for all four and for both FD directions that trigger it) and the mass constraint by -5.1e-4; the other 56 ADVs match the adjoint to 0.5 %. No error or warning. Log signature (base → jumped): `Enrichment - Num Enriched Bases: 6617` → `6621`, `Sweep - All_IP_cells: 5563` → `5565`, `All_IG_cells: 14704` → `14724`, `Ghost - Total Ghost Facets for Bulk Phase 0: 1451` → `1454`, `Total number of DOFs: 6650` → `6652`.
- **Reproducer:** `runs/moris_bugs/MB-031/probe.py` + `decks/probe_default.cpp` (compile with `create_shared_object.sh . build_opt probe`, MORISROOT=moris-gen-additive; ~4 s per eval): design = `runs/studies/gaussian_splat_ls/fp/runs/R1b-geo-N80/results/checkpoints/ADV_Alg_0_Iter_390.hdf5`; `probe.py base` (obj 0.570867999591) vs `probe.py m152 152:-1e-4` (obj 0.571510763554). Bisection: the jump happens exactly when the background node N = (1.0, 0.525) changes sign (ADV152 -3.6e-5: φ(N) = +1.1e-10, no jump; -3.9e-5: φ(N) = -1.3e-10, jump). Working copies: `runs/studies/gaussian_splat_ls/d074/`.
- **Evidence:** `runs/moris_bugs/MB-031/`: `MB-031_ig_mesh_near_node.png` (IG mesh from the real `xtk_temp.exo` of the four states), `log_signature_base_vs_m152.diff`, `cases/{base,m152,ml_base,ml_m152}/` (xtk_temp.exo, GEN exodus, sweep.hdf5, design.hdf5, receipt).
- **Root cause:** a background cell is selected for decomposition only if its corner nodes differ in geometric region (`projects/GEN/src/cl_GEN_Geometry_Engine.cpp:329-352`, called per cell at `projects/XTK/src/cl_XTK_Integration_Mesh_Generator.cpp:476`; the non-intersected list uses the same test at `:540`). An analytic field is then evaluated exactly at derived/child vertices (`projects/GEN/src/ADV/fields/cl_GEN_Field_Analytic.hpp:103-108`), so in a selected cell the regular-subdivision centre vertex sees the true field. At the D-074 design, cell [0.975,1.0]x[0.5,0.525] has corner φ = 7.4e-4, 3.4e-3, 8.2e-3, 3.0e-9 (all positive) but centre φ ≈ -0.028 (back-computed from the diagonal cut points): it is solid while uncut, and when N goes negative the cell is cut and a void diamond covering most of the cell appears at once (figure, top row). Snapping is ruled out: `isocontour_tolerance` = 1e-12 absolute (`cl_GEN_Level_Set_Geometry.cpp:513`), and `intersection_tolerance` is read (`:34`) but never used. Ghost is ruled out: ghost γ 0.005 → 1e-9 leaves the jump (7.2e-4).
- **Code** (open; `projects/GEN/src/cl_GEN_Geometry_Engine.cpp:337-352`, develop-local):
  ```cpp
  // Test nodes for other geometric regions
  for ( uint iNodeInEntityIndex = 0; iNodeInEntityIndex < aNodeIndices.length(); iNodeInEntityIndex++ )
  {
      Geometric_Region tTestGeometricRegion = mGeometries( aGeometryIndex )->get_geometric_region( aNodeIndices( iNodeInEntityIndex ), *( *aNodeCoordinates )( aNodeIndices( iNodeInEntityIndex ) ) );
      if ( tTestGeometricRegion != tFirstNodeGeometricRegion )
      {
          return true;
      }
  }
  // If no differences were found, this entity is not intersected
  return false;
  ```
- **Fix:** 2ae1dfb3f. Cell selection now agrees with the vertices the decomposition evaluates, for exactly evaluated fields only. `Field::has_exact_derived_node_values` is true for `Field_Analytic` (propagated by Combined/Scaled, false for B-spline/stored fields). `Level_Set_Geometry::can_change_region_inside_cell` = exact field and not multilinear.
  - `determine_(non_)intersected_background_cells` also test the regular-subdivision vertices: QUAD4 centre; HEX8 cell centre + 6 face centres (shared, so neighbours agree), computed as N(xi)·X.
  - Before each node-hierarchy geometry pass, `add_background_cells_intersected_at_child_vertices` adds decomposed cells whose current child vertices change region (the MB-016 part).
  - Octree decompositions are not covered.
  - Rejected alternative (b), default multilinear for analytic fields: it changes the response by filling the sub-cell voids (−0.44 % objective here).
- **Verified (D-074 design, `runs/moris_bugs/MB-031/repro2d/`):** run with `run_probe.sh <MORISROOT> <case> <adv:delta>`. The build is `local/mb031-repro` = fix + a cherry-pick of the MB-029 ID fix (be553dd94, local only). Results are in `results_fix.txt`.
  - base obj 0.571092388475 (+0.039 % vs the unfixed 0.570867999591: 4 hidden interior-void cells now cut).
  - Mass jump gone: Δcon −5.1e-4 → −3.8e-9.
  - ADVs 240, 256, 304 no longer jump (Δobj −8.6e-8, −8.4e-9, −1.2e-7).
  - ADV 152 still jumps: +4.62e-4 (was +6.43e-4, −28 %). The jump is still at the flip of node N (−3.6e-5: 0.571092388535; −3.9e-5: 0.571554652463), with +8 IG cells, +1 IP cell and +2 DOFs. This residual is the corner-flip cut-topology change in the neighbouring cells, the same class as the 2.5e-4 left by the multilinear workaround; it is not a detection miss.
- **Code** (fix, 2ae1dfb3f, trimmed):
  ```diff
  # projects/GEN/src/cl_GEN_Geometry_Engine.cpp, new is_intersected( ..., aInteriorPoints ) overload
  +        if ( this->is_intersected( aGeometryIndex, aNodeIndices, aNodeCoordinates ) )
  +        {
  +            return true;
  +        }
  +        // Interpolated fields are bounded by the node values: nothing to find inside (MB-031)
  +        if ( aInteriorPoints.size() == 0 or not mGeometries( aGeometryIndex )->can_change_region_inside_cell() )
  +        {
  +            return false;
  +        }
  +        Geometric_Region tNodeGeometricRegion = mGeometries( aGeometryIndex )->get_geometric_region( aNodeIndices( 0 ), *( *aNodeCoordinates )( aNodeIndices( 0 ) ) );
  +        for ( const Matrix< DDRMat >& iPoint : aInteriorPoints )
  +        {
  +            if ( mGeometries( aGeometryIndex )->get_geometric_region_at_point( iPoint ) != tNodeGeometricRegion )
  +            {
  +                return true;
  +            }
  +        }
  +        return false;
  # projects/XTK/src/cl_XTK_Integration_Mesh_Generator.cpp, determine_intersected_background_cells()
  +                bool tIsIntersected = mXTKModel->get_geom_engine()->is_intersected( iGeom, tGeometricQuery.get_query_entity_to_vertex_connectivity(), tGeometricQuery.get_query_indexed_coordinates() );
  +                if ( not tIsIntersected and tTestInteriorPoints( iGeom ) )
  +                {
  +                    tIsIntersected = mXTKModel->get_geom_engine()->is_intersected( iGeom, tGeometricQuery.get_query_entity_to_vertex_connectivity(), tGeometricQuery.get_query_indexed_coordinates(), tInteriorPoints );
  +                    tNumInteriorOnly += tIsIntersected;
  +                }
  # projects/XTK/src/cl_XTK_Node_Hierarchy_Interface.cpp:78, Node_Hierarchy_Interface::perform()
  +            aMeshGenerator->add_background_cells_intersected_at_child_vertices( aMeshGenerationData, aCutIntegrationMesh, iGeom );
  ```
- **Regression test:** `projects/XTK/test/xtk/UT_XTK_Interior_Cut_Detection.cpp`, tag `[MB-031]` (`[MB-031_2D]`, `[MB-031_2D_JUMP]`, `[MB-016_3D]`, `[MB-016_2D]`). Results: `runs/moris_bugs/MB-031/gates/ut_prefix.log` (develop: 2D void area 0 vs 0.08; area 0 → ~1 at the corner flip; 3D face void missed; 2D two-geometry SIGSEGV) and `ut_postfix.log` (12/12 pass). Full XTK-test: 3285/29 → 3297/33 pass. GEN-test: 81801/28 pass on both.
- **Workaround:** `use_multilinear_interpolation = true` on the analytic design geometry. Measured at the D-074 design: mass jump removed (Δconstraint 5.1e-4 → 4.6e-9), objective jump 6.43e-4 → 2.52e-4 (-61 %); it changes the design's response (obj -0.44 %, mass constraint -5.1e-4 → +4.07e-3, because hidden sub-cell voids are filled). The remaining 2.5e-4 is a cut-topology change at the same node (2 newly cut cells with 1e-16-area void triangles, +4 enriched bases, +2 DOFs; ghost-independent) and is treated as an inherent XFEM discontinuity. Keeping the Gaussian width floor ≥ h should make sub-cell features rarer (this deck: cfg `s_min` = 0.02, h = 0.025; untested).

## MB-032 · SOL/NLA: absolute 1e20 residual guard aborts a well-posed linear adjoint solve
- **Status:** open
- **Symptom:** `*** Reason:  Assertion !( tResNorm > 1e20 ) failed.` / `Convergence::check_for_convergence(): Residual Norm has exceeded 1e20` (cl_NLA_Convergence.cpp:98), stack through `Monolithic_Time_Solver::solve_implicit_DqDs` <- `Time_Solver::solve_sensitivity` <- `Workflow_HMR_XTK::compute_dcriteria_dadv`, i.e. in the adjoint at the first Newton iteration (forward residual healthy, ~3e-3). Log line before the abort: `Newton - ResidualNorm: 6.42e+22` (h=0.01) / `7.25e+29` (h=0.005).
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/triage_1e20/run_arms.sh a0_asis` (1 eval, h=0.01, ~1 min, Lbracket_Stress_Xiga.cpp, moris-reinit-accept binary, restart `runs/moris_bugs/MB-032/restart_D077_Iter_23.hdf5`, XIGA_MODE=compliance XIGA_AGG=ks XIGA_KS_FORM=sum XIGA_KS_RHO=50). Control `a1_pnorm` (only XIGA_AGG=pnorm changed; same design, se 0.0145676, vol 0.185468): adjoint residual 0.22, rc 0.
- **Evidence:** `runs/moris_bugs/MB-032/` (SUMMARY.txt, per-arm env + receipts + log excerpts, D-107 probe excerpt). Toggles that do NOT remove it: XIGA_ANALYTIC_FWD=0 (5.3e22), XIGA_SDF_PRESERVE unset (6.4e22). Occurrences: D-107 (h=0.005, eval 1, sigma_agg 1.82), D-077 (h=0.01, eval 24, sigma_agg 1.57), same signature D-057 (KS 2.2) and D-040 (KS arm, severed design).
- **Root cause:** the trigger is the deck: the KS stress IQI (`IQI_Max_Stress` with `KSParameter`, integrand exp(rho(sigma/sigma* - shift)), shift defaults to 1.0, cl_FEM_IQI_Max_Stress.cpp:174/201) stays in GEN `IQI_types` (criterion 0) in compliance mode. With rho=50 and the aggregate at ~1.5-2.4 sigma*, the integrand is e^25-e^70, so its dQ/du column of the adjoint RHS is 1e22-1e30. MORIS forms and solves an adjoint for every criterion, even though compute_dobjective_dcriteria/compute_dconstraint_dcriteria weight criterion 0 by zero. The MORIS defect: the NLA convergence check applies an absolute magnitude cap (1e20) to the residual of the linear adjoint system, whose scale is set by the criterion, not by divergence; a finite, well-posed linear solve is aborted. (The forward value 8.6e25 passes because IQI values are not checked.)
- **Code:** `projects/SOL/NLA/src/cl_NLA_Convergence.cpp:97-99` (moris-develop-local; :98 in moris-reinit-accept)
  ```cpp
  // check that residual does not exceed upper limit (should defined in MORIS global definitions)
  MORIS_ERROR( !( tResNorm > 1e20 ),
          "Convergence::check_for_convergence(): Residual Norm has exceeded 1e20" );
  ```
- **Fix:** (to write) apply the absolute cap only to the forward nonlinear solve (or make it relative to the reference norm / a parameter), never to the linear adjoint; optionally skip adjoints of criteria whose dobjective/dconstraint weights are identically zero.
- **Regression test:** to add.
- **Workaround:** in stages where the stress aggregate is monitor-only (XIGA_MODE=compliance), run with `XIGA_AGG=pnorm` (and drop XIGA_KS_FORM/XIGA_KS_RHO); verified by arm a1_pnorm. For KS stages that can run far above sigma* (no inner loop, infeasible marches), pass a `Shift` property near the expected peak (e.g. the current aggregate) so the integrand stays O(1).

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
- **Regression test:** `GEN-test.exe "[linear intersection]"` (`ut_GEN_Intersection_Level_Set.cpp`) checks D(A) = A·D(I) for every derived node with the non-commuting A = [[1, 2], [-3, 0.5]]; on `a3b602862` it fails for nodes 21 and 22 (2 of 510 assertions), with the fix it passes. The three-level fixed-branch coordinate FD test is `GEN-test.exe "[MB-022]"` (`Nested Linear Intersections Finite Difference`, PR #32, branch `test/mb022-three-level-fd`): P1/P2 → Q → R on ADV-driven lines vs central FD of the coordinates; max rel. error ≤ 1.7e-10 with the fix, 0.30 at R (9 failed assertions) with the fix reverted. Logs: `runs/moris_bugs/MB-022/fd_test/`. Writing it exposed MB-030.
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
- **Fix:** after exhausted backtracking, `Algorithm_MMA::restore_accepted_state` re-evaluates the held design `xval` (one extra evaluation) so the shared Problem state (ADVs, geometry, analysis, criteria) matches it before `continue`. In reinit mode the restore evaluates R(xval), so the Problem holds that redistance with consistent criteria. If the restore evaluation is itself non-finite (the evaluation is not reproducible) it is a `MORIS_ERROR` (7a947d669; 307236d81 returned silently, which the Manager could not see), so the run exits nonzero for the sentinels.
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
  +                    return istop;    // restart requested (7a947d669; a non-finite restore is a MORIS_ERROR)
  +                }
  +                MORIS_LOG_INFO( "Native MMA iter %d: holding previous design and continuing.", (int)iter );
  +                continue;
  ```
- **Regression test:** `OPT-test.exe "[MB-024]"`: NaN at every design other than the initial x = 0.5 (a design-dependent failure, unlike the audit harness, which also failed the accepted design and so could not be restored). Asserts 18 evaluations (initial, trial, 15 backtracks, restore), no gradient on a failed evaluation, and the Problem returning x = 0.5 with objective 0.25. Fails on `a3b602862`, passes with the fix.
- **Workaround / scope:** Treat exhausted-backtracking output as a failed solve and recover from a known accepted state. max_inner_its=0 does not avoid this branch; the default recipe can reach it if all retries fail.
- **Follow-ups (code review, not changed):** (1) a held iteration restarts the next outer iteration at the full move limit from the same history, so a design-local failure costs 1 trial + 15 backtracks + 1 restore per outer iteration until `max_its`; carrying the shrunk move scale (or a consecutive-hold limit) across iterations would stop the burn. (2) With `max_inner_its > 0`, an inner trial that exhausts its backtracks discards the whole iteration even if an earlier inner trial was finite; committing the last finite trial instead would need a re-evaluation and subproblem multipliers that match it for `kktcheck`. Both are retry-policy changes in the failure path, deliberately left out of this correctness fix.


## MB-025 · OPT: SWEEP single finite-difference modes silently skip all perturbations
- **Status:** fix-written (`fix/mb022-025-audit @307236d81`; found by the implementation audit 2026-10-10 at `a3b602862`). Not yet in the main `build_opt`.
- **Symptom:** forward, backward, and central each perform only one baseline criteria evaluation and produce HDF5 files with no gradient datasets. The all control performs three evaluations and writes analytic plus all three FD gradients.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/implementation_audit_20261010/opt_invalid_audit.cpp`, Catch test `AUDIT SWEEP individual FD modes perform perturbations`. Build and run instructions in that directory's `README.md`.
- **Evidence:** `runs/moris_bugs/MB-025/`; native log `opt_sweep_test.log`; defect intake D-031.
- **Root cause:** `projects/OPT/src/cl_OPT_Algorithm_Sweep.cpp:291`. The parameter API documents all three individual modes, but sweep_solve executes finite differences only under mFiniteDifferenceType == all. No alternative dispatch handles individual modes. The manual already notes this defect; the audit independently reproduces it.
- **Fix:** `sweep_solve` dispatches `forward`/`backward`/`central` to `SA_Type::forward/backward/central` and saves `fd_<mode>` datasets. Single modes are FD-only (analytic gradients only for `none`/`all`, the original intent; 307236d81 always computed them, reverted in 7a947d669 after code review so FD-only checks work on decks whose analytic sensitivities are unavailable or expensive); unknown `finite_difference_type` strings abort in the constructor; `all` writes `fd_central` only for the gradient types requested. `none`/`all` output is otherwise unchanged (no in-repo deck uses an individual mode).
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
- **Regression test:** `OPT-test.exe "[MB-025]"` (`[MB-025][MB-026]`): for none/forward/backward/central/all at two points and two epsilons, asserts the criteria-evaluation count and that every requested dataset exists and equals the exact quadratic FD value (f' ± 100 eps, central f') to 1e-6, and that the analytic dataset exists only for none/all. Fails on `a3b602862`, passes with the fix.
- **Workaround / scope:** Use finite_difference_type=all and inspect only the central-difference datasets (on builds before 7a947d669 the `all` forward/backward datasets are also mislabelled, MB-026). The manual reported sweeps use this workaround, so their results are not invalidated by this finding.


## MB-026 · OPT: FD perturbation sign toggles on every backward call
- **Status:** fix-written (`fix/mb022-025-audit @7a947d669`; found by the code review of that branch, 2026-10-10). Present since before `a3b602862`. Not yet in the main `build_opt`.
- **Symptom:** In a SWEEP with more than one (evaluation point, epsilon) pair, every second forward/backward FD call computes the opposite one-sided scheme and saves it under the wrong label. In `all` mode `fd_forward` and `fd_backward` are swapped for the 2nd, 4th, … pair; in `backward` mode every second `fd_backward` is a forward difference. `fd_central` is unaffected (its formula is sign-symmetric), as is the analytic gradient. No error or warning.
- **Reproducer:** `OPT-test.exe "[MB-026]"` against `cl_OPT_Algorithm.cpp` from `a3b602862`: 6 of 63 assertions fail, e.g. `eval_1-2 epsilon_2-2 fd_forward` = 9.99 where the forward difference of 100 (x - 0.45)^2 at x = 0.5, eps = 1e-4 is 10.01; `eval_1-2 epsilon_2-2 fd_backward` = 10.01 (expected 9.99); the same at x = 0.6 (29.99/30.01 swapped).
- **Root cause:** `projects/OPT/src/cl_OPT_Algorithm.cpp:327` (`a3b602862`). `set_sensitivity_analysis_type( SA_Type::backward )` negates `mFiniteDifferenceEpsilons` in place and `SA_Type::forward` does not reset it, so the sign toggles with each backward call and leaks into the next forward call. SWEEP is the only caller that selects forward/backward.
- **Fix:** set the sign absolutely: backward `-|eps|`, forward `+|eps|` (central uses whatever sign; its quotient is invariant).
- **Code** (fix, `git show 7a947d669`, trimmed):
  ```diff
               case opt::SA_Type::backward:
  -            {
  -                mFiniteDifferenceEpsilons = -1.0 * mFiniteDifferenceEpsilons;
  -                [[fallthrough]];
  -            }
               case opt::SA_Type::forward:
               {
  +                const real tSign = aType == opt::SA_Type::backward ? -1.0 : 1.0;
  +                for ( uint iEntry = 0; iEntry < mFiniteDifferenceEpsilons.numel(); iEntry++ )
  +                {
  +                    mFiniteDifferenceEpsilons( iEntry ) = tSign * std::abs( mFiniteDifferenceEpsilons( iEntry ) );
  +                }
  ```
- **Regression test:** `OPT-test.exe "[MB-026]"` (shared with MB-025): exact forward/backward/central values at two points and two epsilons.
- **Affected outputs:** any pre-fix SWEEP output read for `fd_forward`/`fd_backward` at an even-numbered call. In this workspace: `runs/studies/lunch_and_learn_xiga/triage_bridge/scripts/fd_table.py:12` prints `fwd e2`/`bwd e2` from `eval_1-1 epsilon_2-4`, the swapped pair, so its two one-sided columns are exchanged (its central columns and the analytic gradient are correct). The 27 in-repo EXA SWEEP decks that use `all` may store swapped forward/backward gold values for their 2nd+ epsilon.
- **Workaround (builds before 7a947d669):** use `fd_central` only, or one epsilon and one evaluation point per sweep.


## MB-027 · WRK/OPT: non-idempotent in-place SDF redistance defeats the GCMMA inner loop
- **Status:** merged (#34, squash 4900ba1a8 on develop, 2026-10-13; PR #20 landed first as 10aaa6d2b, then the branch was refreshed onto it and merged). Previously: fix-written, PR #34 (https://github.com/kdoble48/moris/pull/34, branch `fix/mb027-preserve-interface`, base develop de7990fc2). The PR carries the three MB-027 commits cherry-picked from `fix/reinit-inner-consistency-xdb` (`90f189d2d`, `96e68a44b`, `67911d133`, the rebased copies of `9b35e3a32`, `1760b6e62`, `ce64ec9ab`) plus their prerequisites, which never reached develop: PR #20's slice-write commit `d50ede932` (its XTK `get_index()` hunk dropped) and MB-017/018/019 (`7a3115ee7`, `bdcd33884` without the ledger, `a0fe9502e`); it supersedes PR #20. The fast SDF distance `f0729c12b` is not in the PR (own branch `perf/sdf-bucket-grid`, pushed, no PR yet). Gate: new tests pass; GEN/XTK/INT/SDF/WRK/OPT suites at 1/2/4 procs have the same failures as develop (WRK_sensitivity_test missing .so; OPT [MB-023]/[MB-025] at 2/4 procs). The original dev line `fix/reinit-inner-consistency` (worktree `moris-reinit-accept`, on feature/ks-stress-iqi) is unchanged.
- **Symptom:** native GCMMA with `reinit_in_place` + `sdf_reinit_every_step` + `max_inner_its` 10 (X18 v8 KS mimic, h = 0.01): 145 of 558 iterations hit the inner-loop cap (11 evaluations each, 1598 of 3760 evaluations). In capped stretches (iterations 496-519) the trial step goes to ~0 but every trial reads the same constraint, above the approximation, e.g. iteration 504 trials 1.05306 → 1.05350 against an accepted 503 value of 1.04467; the accepted design's KS constraint climbs 1.01 → 1.14 monotonically while the volume rises 0.2291 → 0.2318. Signature: runs of `n_evals = max_inner_its + 1` in `studies/lunch_and_learn_xiga/scripts/inner_loop_stats.py <console.log> --cap 11`.
- **Reproducer:** `runs/studies/lunch_and_learn_xiga/perf/x18v8/rchain` (queue job 3104): restart from the iteration-500 design (`perf/x18v8/restarts/ks_v2_iter500_ckpt3270.hdf5`), MMA step 1e-8, no inner loop, 8 iterations. The design does not move, so evaluation k reads f(R^(k-1) x): KS aggregate pn 1.02362 → 1.02997 → 1.03636 → 1.04286 → 1.04930 → 1.05563 → 1.06177 → 1.06778 → 1.07368 (+0.0064 per application) and volume 0.147839 → 0.148721 (+0.00011 per application).
- **Evidence:** `runs/studies/lunch_and_learn_xiga/perf/x18v8/` (PROFILE.md "Inner-loop root cause"; arms `ab_inner/a_asis` and `ab_inner/b_sdf_off`). With the redistance off (b), the same 30 iterations from iteration 495 take 5.13 evaluations per iteration (6 capped) and lower the volume 0.2279 → 0.2049 at a constraint of 1.000; as-is (a) takes 6.13 (8 capped) and reaches 0.2226. The arm with no inner loop (c) diverges (KS 1.03 → 2.2, residual > 1e20 at evaluation 7, D-057), so the inner loop is needed; the redistance is what defeats it.
- **Root cause:** (1) `Reinitialize_Performer::compute_sdf_from_interface` (nodal distance to the linear-chord cut facets, then an L2 projection onto the quadratic B-splines) does not reproduce the zero contour of an already-signed-distance field, so R(R(x)) != R(x); the drift is systematic (thickening, higher peak stress). (2) `Workflow_HMR_XTK::perform` applies R on every criteria evaluation, including inner-loop trials, while `Algorithm_MMA::is_conservative` compares the criteria of R(xmma) with the approximation built for xmma. When the subproblem step is small the R drift dominates, no trial is conservative, and the trial accepted at the cap is R(x_k): the run then iterates x ← R(x).
- **Code:** `projects/WRK/src/cl_WRK_Workflow_HMR_XTK.cpp:444-445` (moris-xiga-ks, feature/ks-stress-iqi fbf62f3cb)
  ```cpp
  bool tEveryStep  = tReinitPerformer->get_reinit_every_step();
  bool tRedistance = tEveryStep ? ( tOptIter > 0 && tOptIter % tReinitFreq == 0 ) : tScheduled;   // every evaluation, trials included
  ```
  `projects/OPT/src/cl_OPT_Algorithm_MMA.cpp:222-227`
  ```cpp
  if ( ( f0app + epsimin ) < f0valnew ) return false;                 // f0valnew, fvalnew are criteria of R(xmma)
  for ( int j = 0; j < mNumCon; ++j )
      if ( ( fapp[ j ] + epsimin ) < fvalnew[ j ] ) return false;     // fapp is the approximation at xmma
  ```
- **Regression test:** `[sdf_preserve_interface][MB-027]` (WRK-test, `UT_WRK_Reinitialize_Preserve.cpp`): 8 redistance-like applications (sign-preserving nodal distance targets + L2 projection, as `compute_sdf_from_interface`) to an ellipse level set on 24x24 HMR at B-spline orders 1 and 2; with the freeze band the cut cells, every nodal value on them and every edge intersection point are bit-identical, without it the cut moves (negative control). `[reinit_on_accept_only]` (OPT-test): with the inner loop on, every trial is evaluated unredistanced and every redistanced evaluation re-evaluates the accepted trial. `[sdf_taper]` (WRK-test, OPT-test): wean monitor (tolerance + patience, reset, stop iteration, capped iterates not counted) and the MMA accepted-iterate notification.
- **Fix:** `sdf_preserve_interface` (OPT problem list, default false) resets to psi every coefficient nonzero on a cell whose psi coefficients are not all of one strict sign (a superset of the cells any cut test can intersect, by the B-spline convex-hull property), then grows that set until every other cell keeps a strict sign; the cut is then bit-identical and R is idempotent on everything the criteria see. `reinit_on_accept_only` (MMA list) evaluates the inner trials unredistanced and redistances only the accepted iterate (one extra evaluation per iteration); on its own it still drifts once per accepted iterate.
  ```diff
  +                                if ( tReinitPerformer->get_preserve_interface() )
  +                                {
  +                                    tReinitPerformer->preserve_interface( mPerformerManager->mGENPerformer, tOldCoeffs, tCoeffs );
  +                                }
  ...
  +        // 1. possibly cut cells (psi coefficients not all of one strict sign): freeze all their coefficients
  +            tCellSign[ iCell ] = tAllPos ? 1 : ( tAllNeg ? -1 : 0 );
  +            if ( tCellSign[ iCell ] == 0 )
  +                    tFrozen[ tCoeffs( iC ) ] = true;
  +        // 2. strictly signed cells must stay strictly signed with the same sign; grow the frozen set until they do
  +                    tOk         = tCellSign[ iCell ] > 0 ? tValue > aTolerance : tValue < -aTolerance;
  ```
- **Verification (MBB 60x20 testbed, `runs/studies/lunch_and_learn_xiga/reinit_rootcause/`, 120 its from X10_a60_d30_rcp Iter_30, max_inner_its 10):** inner-loop cap hit in 93 % of iterations as-is (10.4 evaluations per iteration, constraint ratchets to +0.018 with resets at iterations 44, 74, 119, final g = +0.011), 16 % with `sdf_preserve_interface` (4.2 evaluations per iteration, g = 0, objective 1.0961, merit gate CONVERGED at iteration 115), 19 % with the redistance off. Fixed design redistanced 8 times: criteria unchanged to all printed digits with the band (band max |dc| = 0), volume +0.78 % without. Per-trial J(R(psi)) - J(psi) = 0.0 exactly in all 60 trials of the preserve debug arm, against 4e-4 / 1.3e-3 (objective / constraint) as-is.
- **Verification (L-bracket KS mimic h = 0.01, `runs/studies/lunch_and_learn_xiga/reinit_rootcause/lbracket/`, 40 its from the iteration-495 checkpoint, fast switches):** as-is 45 % capped, 7.75 evaluations and 33.5 s per iteration, g climbs to +0.11; `sdf_preserve_interface` 15 % capped (redistance off: 12.5 %), 4.7 evaluations and 13.2 s per iteration, V 0.228 -> 0.210; accept-only + preserve gives the identical trajectory at one more evaluation per iteration. Fixed design (iteration 500) redistanced 8 times: volume and KS aggregate unchanged (0.147839, 1.02362) with the band (2690 of 10404 coefficients frozen, order 2), KS +0.0064 per application without (9 nodal sign flips, 6 cut cells changed). Full run from the S_stress seed (`X18/v8/mimic_h1/ks_v2_fixP`): 21/591 iterations capped (parent 145/558), 2399 evaluations (parent 3753), no constraint ratchet (max c·A over the last 300 iterations 1.029, parent 1.143); that run also used XIGA_STRESS_DESIGN_ONLY=1, so its volume (0.188 vs 0.229) is not a like-for-like comparison.
- **Residual:** with the redistance off the inner loop still caps 12-19 % of iterations (L-bracket and MBB): a separate, non-ratcheting approximation miss (first-iteration KS approximation error ~0.018 after a restart), not part of MB-027.
- **Workaround:** until merged, `XIGA_SDF_PRESERVE=1` (L-bracket) / `sdf_reinit.preserve_interface` (MBB) on the moris-reinit-accept binary; on other binaries the deck switch `XIGA_SDF_FREQ=0` (Lbracket_Stress_Xiga.cpp, 2026-10-11) turns the in-place redistance off and keeps the recipe's phases and native MMA; `XIGA_SDF_FREQ=k` redistances every k-th evaluation only.

## MB-038 · MTK/HMR: GCC 13 misaligned stores in base ctors of classes with a virtual `mtk::Mesh` base
- **Status:** fixed 86ce7f9b2 on `perf/vector-move-semantics` (PR #52 open); built, XTK/GEN/CNT/INT/OPT/HMR/MTK suites pass (2026-10-20).
- **Symptom:** SIGSEGV (si_code -6 in the OpenMPI handler) in `moris::mtk::Integration_Mesh::Integration_Mesh()` ← `Integration_Mesh_STK` ctor ← `create_integration_mesh_from_interpolation_mesh` ← `Geometry_Engine` ctor; XTK-test aborts in its first case "XTK Cell Clusters". Faulting instruction `movaps %xmm0,0xd0(%rdi)` with `rdi` = STK `this` + 0x28.
- **Reproducer:** `runs/moris_bugs/MB-038/gcc13_virtual_base_movaps_repro.cpp` (standalone, `g++-13 -O3`, crashes with or without move members); in MORIS, XTK-test on the move-semantics branch without the alignas change.
- **Evidence:** `runs/moris_bugs/MB-038/` (repro, alignas variant that runs, XTK segfault log, README with the gdb registers).
- **Root cause:** `mtk::Mesh` has `alignof == 16` (armadillo `Mat` local buffer). A class with it as a *virtual* base has a 16-aligned complete type but an 8-aligned non-virtual part, so the Itanium ABI may place it at 8 mod 16 inside a derived class (`Integration_Mesh` at 0x28 in `Integration_Mesh_STK`). GCC 13 vectorizes the zeroing in the base-object ctor (C2) with aligned stores assuming `alignof` of the complete type. Latent on develop: there the ctor happens to be inlined; a change that makes GCC emit and call the out-of-line C2 (here the moris::Vector move members) exposes it.
- **Fix:** `alignas( alignof( Mesh ) )` on every class with a direct virtual `Mesh` base (`mtk::Integration_Mesh`, `mtk::Interpolation_Mesh`, `mtk::Mesh_Core_STK`, `hmr::Mesh`, `hmr::Interpolation_Mesh_HMR`); the non-virtual part is then 16-aligned and placed accordingly (`Integration_Mesh` moves to 0x30). Descendants inherit it.
- **Code** (open, `projects/MTK/src/mesh/cl_MTK_Integration_Mesh.hpp:29` on 97c863cb7):
  ```cpp
  class Integration_Mesh : public virtual Mesh
  ```
- **Regression test:** XTK-test "XTK Cell Clusters" (crashes without the fix on the move-semantics branch).
- **Workaround / scope:** any class with an over-aligned virtual base; `Amesos_GenOp : virtual Epetra_Operator` is not over-aligned.

## MB-039 · GEN: `gen::ADV` holds a raw pointer into the ADV vector
- **Status:** fixed 4b6a87b31 on `perf/vector-move-semantics` (PR #52 open); built, XTK/GEN/CNT/INT/OPT/HMR/MTK suites pass (2026-10-20).
- **Symptom:** GEN-test Circle/Sphere/User-defined Geometry return values of the old ADVs or garbage (e.g. Circle φ = 1.236 vs 0) after `tADVManager.mADVs = { ... }`, then SIGSEGV; surfaced when moris::Vector gained move assignment (the braced-list temporary's buffer is moved in).
- **Reproducer:** `GEN-test.exe "Circle"` on the move-semantics branch without the ADV change; `runs/moris_bugs/MB-039/gen_field_circle_fail.log`.
- **Root cause:** `ADV::ADV( Vector< real >&, sint )` stores `&aADVs( aID )` (`projects/GEN/src/ADV/cl_GEN_ADV.cpp:20`), valid only while the vector's buffer never moves. Also `Geometry_Engine` ctor (`cl_GEN_Geometry_Engine.cpp:71-72`) copies then `clear()`s `mADVManager.mADVs` before `distribute_advs()` repoints the designs, so initial refinement reads cleared storage (UB that worked because `clear()` keeps the buffer).
- **Fix:** `ADV` stores the owning `Vector< real >*` and reads `( *mADVs )( mID )` on access (Dist_Vector and constant ADVs unchanged); `Geometry_Engine` no longer clears `mADVManager.mADVs`.
- **Code** (open, `cl_GEN_ADV.cpp:19-23` on 97c863cb7):
  ```cpp
  ADV::ADV( Vector< real >& aADVs, sint aID )
          : mValue( &aADVs( aID ) )
          , mID( aID )
  ```
- **Regression test:** GEN-test `ut_GEN_Field.cpp` cases Circle, Sphere, User-defined Geometry (reassign `mADVs` after creating the fields).

## MB-036 · GEN: nested-intersection d(xi)/d(parent) uses the edge-projected field gradient
- **Status:** open (found 2026-10-11 while writing the MB-022 FD regression). Renumbered from a duplicate MB-030 on 2026-10-19; the hidden Catch case on develop (PR #32) is still tagged `[MB-030]`, so rename the tag to `[MB-036]` in the fix PR.
- **Symptom:** Silent. Analytic intersection-coordinate sensitivities (`Node_Manager::append_dcoordinate_dadv_from_derived_node`) of a level ≥ 3 intersection node disagree with finite differences of its coordinates in the columns of ADVs that move a parent off its own parent edge: max relative error 0.38 (level 3) and 0.64 (level 4) in the reproducer. Levels 1 and 2 are exact. Example: node R, d R_x/d c0x analytic 0.486 vs FD 0.131.
- **Reproducer:** `GEN-test.exe "[MB-030]"` (hidden case `Nested Linear Intersections Finite Difference Off-Edge Parent`, `projects/GEN/test/ut_GEN_Intersection_Level_Set.cpp`, PR #32): four ADV-driven LINE geometries on one QUAD4, chain P1/P2 (line 0) → Q (line 1 on P1→P2, Q lies on line 0) → R (line 2 on Q→corner) → T (line 3 on P1→R).
- **Evidence:** `runs/moris_bugs/MB-036/`: `mb036_offedge_develop.log` (Catch output on develop de7990fc2 + this test, 26 failed assertions), `projected_gradient_repro.py/.txt` (independent numpy chain: with the edge-projected gradient it reproduces MORIS's analytic 0.4860, -0.0582 exactly; with the full gradient it reproduces the FD 0.1307, -0.4357).
- **Root cause:** `projects/GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:598-630`, `get_dfield_dcoordinates( aParentNode, ... )`, called from `Intersection_Node_Linear::get_dxi_dcoordinate_{first,second}_parent`. For an intersection-node parent it returns Δφ·v/|v|² over the parent's own locator edge v, i.e. the projection of ∇φ onto v. The chain term d(xi)/d(X_parent)·dX_parent/dADV is then exact only when dX_parent ∥ v. A level-2 node's parents sit on fixed background edges, so this holds there (existing tests pass); from level 3 a parent's edge itself moves (e.g. Q on a moving interface), and the component of ∇φ normal to v is dropped.
- **Fix:** not written. For analytic fields the full gradient is available (`Field_Analytic::get_dfield_dcoordinates`); for discrete/interpolated fields the consistent derivative of the derived-node value has to be defined (needs a theory decision on how φ at a derived node is represented).
- **Code** (open, `cl_GEN_Level_Set_Geometry.cpp:625-629` on de7990fc2):
  ```cpp
  // Sensitivities: dPhi/dx_i  = delta(Phi) / L_i where L_i = PaerentVectorLenth^2 / (ParentVector * e_i)
  for ( uint tCoordinateIndex = 0; tCoordinateIndex < tNumDim; tCoordinateIndex++ )
  {
      aSensitivities( tCoordinateIndex ) = tDeltaPhi * tParentVector( tCoordinateIndex ) / ( tParentLengthSquared + MORIS_REAL_EPS );
  }
  ```
- **Regression test:** `GEN-test.exe "[MB-030]"` (tag as merged in PR #32; rename to `[MB-036]` with the fix), hidden (`[.]`) until fixed; remove the `[.]` with the fix.
- **Workaround / scope:** Affects only intersections on intersections on intersections (≥ 3 levels of ADV-dependent nesting) whose intermediate parent moves transversally to its edge. Single-geometry XTK cuts and two-level nesting are unaffected.


## MB-033 · WRK: remeshing restart uses a stale MDL field after a failed XTK decomposition
- **Status:** merged (#38, landed on develop as 2c99051e1). Fixed 63c42e266 on `fix/remesh-stale-field-mb033`; built, WRK-test + reproducer pass.
- **Symptom:** in an adaptive-remeshing optimization that feeds an FEM field to the refinement (ZZ error indicator), a
  second consecutive `xtk::Model::decompose() - Decomposition Failed. Not all child meshes on the same level.` is followed,
  on the next remesh, by SIGSEGV in `Refinement_Mini_Performer::perform_refinement_2` ←
  `Remeshing_Mini_Performer::perform_refinement_mode_2` ← `Workflow_HMR_XTK::initialize` ← `opt::Problem::initialize` ←
  `opt::Manager::reinitialize` (gdb: `rax=0`, the parent climb ran past the level-0 root; the logged failing address is
  meaningless, si_code -6).
- **Reproducer:** `runs/docs_gallery_runs/zz_triage/r0_baseline_gdb/run.sh` (stock `Zienkiewicz_Zhu_Adaptive_Refinement.so`,
  `MORIS_EXA_MODE=showcase`, 1 rank, ~2 min): 27 evaluations, decomposition failures at evaluations 24/24/25, crash on the
  third remesh.
- **Evidence:** `runs/moris_bugs/MB-033/` (objective history with the remesh requests and failures, gdb excerpt,
  instrumented-run excerpt). The instrumented deck (`e1`) shows the GEN fields moving to the new mesh (15270 cells) while
  the ZZ field stays on the old one (`0x55555bc42720`, 5400 cells).
- **Root cause:** `Workflow_HMR_XTK::initialize()` appends `tModelPerformer->get_mtk_fields()` to the refinement fields
  unconditionally, but the MDL model's fields are only rebuilt by a successful forward solve; an evaluation that dies in
  XTK leaves them on the mesh from before the previous remesh. `perform_refinement_2` then maps that mesh's cells through
  a reference pattern that does not cover them: `tMap.find` can miss and the `get_parent()` climb runs off the root.
- **Fix:** `Workflow_HMR_XTK` tracks `mMDLFieldsOnCurrentMesh` (set after the MDL solve, cleared by a remesh) and passes
  MDL fields to the refinement only when it is set; otherwise it remeshes on the geometry fields with a warning. A
  scheduled remesh follows a solve, so normal runs are unchanged. `perform_refinement_2` skips an absent field and turns
  the unguarded map lookup and parent climb into `MORIS_ERROR`s naming the field. (`create_refinement_input_list_2` builds
  the refinement list from the passed fields by position; the MDL fields come after the GEN fields, so withholding them
  drops only trailing entries.)
- **Code** (fix, trimmed from `git show 63c42e266`):
  ```diff
  @@ cl_WRK_Workflow_HMR_XTK.cpp -207,9 +207,17 @@
  -            tFieldsIn.append( tModelPerformer->get_mtk_fields() );
  +            if ( mMDLFieldsOnCurrentMesh )
  +            {
  +                tFieldsIn.append( tModelPerformer->get_mtk_fields() );
  +            }
  +            else
  +            {
  +                MORIS_LOG_WARNING( "Workflow_HMR_XTK::initialize() - no forward solve on the current mesh since the last remesh; "
  +                                   "remeshing on the geometry fields only" );
  +            }
  @@ -223,6 +231,9 @@
  +            // the MDL fields now live on the replaced mesh
  +            mMDLFieldsOnCurrentMesh = false;
  @@ -688,6 +699,7 @@
           mPerformerManager->mMDLPerformer( 0 )->perform();
  +        mMDLFieldsOnCurrentMesh = true;
  @@ cl_WRK_perform_refinement.cpp -238,6 +248,11 @@ / -246,6 +261,12 @@
  +                    MORIS_ERROR( tIter != tMap.end(),
  +                            "perform_refinement_2(): cell of field '%s' has no background element in the current mesh "
  +                            "(field older than the last remesh?)",
  +                            tFieldName.c_str() );
                       ...
                           tBGElementNew = tBGElementNew->get_parent();
  +
  +                        MORIS_ERROR( tBGElementNew != nullptr,
  +                                "perform_refinement_2(): field '%s' lives on a mesh not covered by reference pattern %d "
  +                                "(field older than the last remesh?)",
  +                                tFieldName.c_str(),
  +                                (int)tThisRefernecePattern );
  ```
- **Regression test / verification:** no unit test reaches this path (it needs a full remeshing optimization with a
  failed decomposition); the reproducer is the check. Fixed build (`moris-remesh-stale-field`), stock deck, showcase mode:
  the third remesh logs the warning, refines on the geometry fields, and the run continues (objective 2.868 -> 2.835;
  `runs/docs_gallery_runs/zz_triage/r1_mb033_fixed/run.log`, stopped at 32 evaluations by a guard because the stock
  schedule asks level 10 for iterations 30-99). EXA `Levelset_Boxbeam_Adaptive_Refinement-{1,4}-procs` and
  `Zienkiewicz_Zhu_Adaptive_Refinement-1-procs` pass (`runs/docs_gallery_runs/GATE6_mb033_remesh_exa_ctest.log`).
- **Workaround:** keep every newly cut cell one remesh away from the interface level: floor `minimum_refinement_level` on
  the interface pattern at (max interface level - 1). The ZZ showcase (level-1 floor) ran 196 evaluations with 17
  decomposition failures, each recovered by a single remesh (`runs/docs_gallery_runs/zz_triage/fA_minlvl1`).
- **Related (suspect, not ledgered):** a restart requested on the first evaluation of a restarted GCMMA is not honoured
  (`aIter = -aIter` is a no-op at 0, `cl_OPT_Algorithm_GCMMA.cpp:246-249`); GCMMA steps once on cached gradients.

## MB-034 · WRK: `minimum_refinement_level` pairs read with stride 1
- **Status:** merged (#38, landed on develop as 2c99051e1). Fixed f0fb28d47 on `fix/remesh-stale-field-mb033`; built, WRK-test [MB-034] passes.
- **Symptom:** silent. Each `minimum_refinement_level` entry is `pattern, level_1, until_1, level_2, until_2, ...`. The
  loop reads pair `Ia` at `(Ia + 1, Ia + 2)` instead of `(2 Ia + 1, 2 Ia + 2)`, so only the first pair is read correctly;
  later lookups straddle pairs. `0, 1, 30, 2, 100` gives level 0 instead of 2 for iterations 30-99.
- **Reproducer / evidence:** `runs/moris_bugs/MB-034/stride_table.md` (three schedule strings, develop vs intended).
- **Root cause** (`cl_WRK_perform_remeshing.cpp:244-248` on develop 97fbd9450), fixed as below.
  ```cpp
  for ( uint Ia = 0; Ia < std::floor( ( (real)mParameters.mMinimumRefinementLevel( Ik ).size() ) / 2.0 ); Ia++ )
  {
      if ( mParameters.mMinimumRefinementLevel( Ik )( Ia + 2 ) > tOptIter )
      {
          tMinRefinement( Ik ) = mParameters.mMinimumRefinementLevel( Ik )( Ia + 1 );
  ```
- **Fix / Code** (trimmed from `git show f0fb28d47`):
  ```diff
  +    uint
  +    minimum_refinement_level_at(
  +            const Vector< uint >& aEntry,
  +            uint                  aOptIter )
  +    {
  +        // aEntry = pattern, level_1, until_1, level_2, until_2, ...: pair k sits at ( 2k + 1, 2k + 2 ) (MB-034)
  +        for ( uint iPair = 0; 2 * iPair + 2 < aEntry.size(); iPair++ )
  +        {
  +            if ( aEntry( 2 * iPair + 2 ) > aOptIter )
  +            {
  +                return aEntry( 2 * iPair + 1 );
  +            }
  +        }
  +        return 0;
  +    }
  @@ -241,15 +259,7 @@
  -                for ( uint Ia = 0; Ia < std::floor( ( (real)mParameters.mMinimumRefinementLevel( Ik ).size() ) / 2.0 ); Ia++ )
  -                {
  -                    if ( mParameters.mMinimumRefinementLevel( Ik )( Ia + 2 ) > tOptIter )
  -                    {
  -                        tMinRefinement( Ik ) = mParameters.mMinimumRefinementLevel( Ik )( Ia + 1 );
  -
  -                        break;
  -                    }
  -                }
  +                tMinRefinement( Ik ) = minimum_refinement_level_at( mParameters.mMinimumRefinementLevel( Ik ), tOptIter );
  ```
- **Regression test:** `WRK-test.exe "[MB-034]"` (`projects/WRK/test/UT_WRK_Remeshing_Mini_Performer.cpp`): two-pair,
  three-pair and bare entries; passes (9 assertions).
- **Scope:** decks with one (level, until) pair are unaffected. The ZZ deck's regression schedule coincides by chance.

## MB-035 · MTK/GEN: singular node-to-B-spline projection when the Lagrange order is below the B-spline order
- **Status:** open (2026-10-13)
- **Symptom:** a GEN level set discretized on B-splines of order 2 (`discretization_mesh_index` 0) with a Lagrange mesh of order 1: at h = 1/32 in 3D (39304 coefficients) the GEN nodal field is garbage, with no error or warning: max |φ_GEN − φ_exact| = 2373 for a field in [−0.58, 0.29], sign wrong at 15.6 % of the nodes, cut interface area +1049 %, sphere volume +94 %, 270992 interface facets instead of ~24700. Log signature in the `MTK - Mapper - Map Node-to-Bspline` block: `Belos - IterativeSolverConverged: 1` (Belos::Unconverged), `Belos - LinearSolverIterations: 1000`, `Belos - LinearResidualNorm_All_RHS: 1.37023002642362e+01` against an initial `Newton - ResidualNorm: 1.15800200738482e-03`. The same setup at h = 1/8 and 1/16 (≤ 25000 coefficients, direct solver) gives a field that matches φ at every Lagrange node to 5e-12. Independent of the intersection mode (linear, multilinear and cell order give the same broken field).
- **Reproducer:** `runs/moris_bugs/MB-035/reproducer.sh <outdir>` (deck `studies/lunch_and_learn_xiga/decks/Sphere3D_Precision.cpp`, `SPH_N=32 SPH_P=2 SPH_LAG=1`, forward only, moris-reinit-accept binary, ~4 min, ~4 GB, 1 rank). Original run: `runs/studies/lunch_and_learn_xiga/sphere3d_precision/n32_p2_l1_linear` (also `_multilinear`, `_cellorder`); defect D-158.
- **Evidence:** `runs/moris_bugs/MB-035/` (`MB-035_gen_field_error.png`: nodal error on the z = 0.5 slice at h = 1/16 vs 1/32 from the real `gen_field.exo`; mapper log excerpts; both `gen_field.exo`).
- **Root cause:** the projection solves an L2 system on the Lagrange (interpolation) mesh, where a B-spline field is represented through the T-matrix as its Lagrange-1 nodal interpolant. The mass matrix is Tᵀ M_L T with T of size (N+1)³ × (N+2)³, so its rank is at most (N+1)³ < (N+2)³: it is singular (null space dimension 3367 at N = 32) whenever the Lagrange order is below the B-spline order. The direct solver happens to return a solution that interpolates the nodal values (that is what the small meshes show); GMRES + ILU(1) on the singular system diverges. `Mapper::map_node_to_bspline` switches solver on size only and never checks the solve.
- **Code** (open; `projects/MTK/MAP/src/cl_MTK_Mapper.cpp:672-695`, develop-local):
  ```cpp
  if ( sum_all( tNumberOfCoefficients ) < 25000 && par_size() < 25 )
  {
      tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL ) );
      ...
  }
  else
  {
      tParameterlist( 0 ).add_parameter_list( moris::prm::create_linear_algorithm_parameter_list( sol::SolverType::BELOS_IMPL ) );
      tParameterlist( 0 )( 0 ).set( "preconditioners", "0" );
      tParameterlist( 7 ).add_parameter_list( moris::prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK ) );
      tParameterlist( 7 )( 0 ).set( "ifpack_prec_type", "ILU" );
      tParameterlist( 7 )( 0 ).set( "fact: level-of-fill", 1 );
  }
  ```
- **Fix:** (to write) (a) GEN: reject (MORIS_ERROR) a B-spline design discretization whose Lagrange mesh order is below the B-spline order, or project on a Lagrange mesh of at least the B-spline order; (b) Mapper: check the linear-solver return status (and the final residual) and abort instead of setting coefficients from an unconverged solve.
- **Regression test:** to add: MTK/GEN, 3D HMR 32³, B-spline order 2 on Lagrange order 1, project a sphere SDF: must raise an error (a) or reproduce the nodal values (b); also a mapper test that an unconverged Belos solve raises.
- **Workaround:** Lagrange order = B-spline order for any B-spline design field (the sphere study's p = 2 / Lagrange 2 arms are clean to h = 1/32). At p = 2 on Lagrange 1 the cut sees only the nodal (corner) values, i.e. the same field as p = 1, even when the solve works.

## MB-037 · XTK/GEN: near-coincident interfaces give inverted IG cells and a half-snapped layer
- **Status:** open (found 2026-10-20 in the XTK walk-through pathology sweep).
- **Evidence:** `runs/moris_bugs/MB-037/parallel_planes_sweep.log` (sweep table). Builder for
  `studies/lunch_and_learn_xiga/scripts/build_bug_evidence.py`: to add (render of the δ = 1e-11 cut mesh).
- **Symptom:** one HEX8 [0,1]³, two parallel planes A = n·(x − p) and B = A − δ (n ∝ (1, 0.3, 0.2), stored regime,
  `isocontour_tolerance` 1e-12). For δ = 0.1 … 1e-9 the layer volume is exact and the smallest cell volume falls as δ³.
  At **δ = 1e-11, 4 IG tets have volume ≤ 0** (smallest −4.8e-20): inverted cells. At **δ = 1e-12** the B-values at
  A's cut vertices fall inside the tolerance, so those vertices are labelled INTERFACE for B. The layer is then
  partly collapsed: it keeps 7.6e-14 of its 1.06e-12 volume (93 % lost) and still forms its own subphase. At δ ≤ 1e-13
  the layer snaps away cleanly. No error or warning in any case.
  **It needs only one geometry:** corner values in {−1, −δ, +δ, 1} with δ = 1e-11 give cells with volume ≤ 0 in
  29 962 of 65 536 corner patterns (46 %; smallest −1.9e-18, first at corners [1e-11, −1, …, −1]), and slivers
  (< 1e-12) in 63 906. One interface plane 1e-11 off the face/cell-centre subdivision vertices (`x=0.5+1e-11`) gives
  inverted cells in 13 072 of 13 122 two-geometry cases, whatever the other geometry
  (`near_zero_sweep.log`; `runs/studies/xtk_pathology/2026-10-20/XTK_PATHOLOGY_HEX_TWO_GEOMETRIES.csv`).
- **Impact (measured 2026-10-21):**
  - Quadrature: the inverted cells weigh ~1e-18; phase moments (1, x, y, z, x²) are unchanged with signed or absolute
    cell volumes (single-hex sweep, `[XTK_PATHOLOGY_HEX_NEAR_ZERO]`).
  - FEM patch test (2×2×2 HEX8, same material in every phase, exact linear displacement imposed by Nitsche, Nitsche
    interfaces; `runs/studies/xtk_pathology/patch_test/results.md`): the displacement stays exact (L2 error ≤ 3.4e-13).
    But the near-coincident layer makes the system near singular: condition number 1.7e13 (δ = 1e-9), 1.7e15 (1e-11),
    1.2e17 (1e-12) vs 3.3e5 for the control, ghost on or off; the Newton residual stalls at 2e-5..1.7e-4 of its start
    and the run still reports success (D-183). Element-averaged stress output is 0 on 40–50 % of the cut cells whenever
    degenerate cells are present (D-182). At δ ≤ 1e-13 (layer snapped away) all of this disappears.
- **Reproducer:** `XTK-test.exe "[XTK_PATHOLOGY_HEX_PARALLEL]"` (`projects/XTK/test/xtk/UT_XTK_Cut_Pathologies.cpp`,
  branch `docs/xtk-walkthrough`, commit ca29cad68; a hidden `[.]` test).
- **Root cause (from code reading):** the INTERFACE label is an absolute test on the field value,
  `|φ − c| ≤ isocontour_tolerance` (`GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:669-686`). A cut point is admitted at
  any distance from an existing vertex (`XTK/src/cl_XTK_Node_Hierarchy_Interface.cpp:235-283`). The distance tolerance
  `intersection_tolerance` is parsed into `mIntersectionTolerance` (`GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:34`) and
  never used. A second geometry's cut vertex can therefore land within roundoff of a first geometry's vertex: the tets
  between them invert, or the layer is labelled away vertex by vertex.
- **Code:**
  ```cpp
  // GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:672
  if ( std::abs( aLevelSetValue - mParameters.mIsocontourThreshold ) <= mParameters.mIsocontourTolerance )
  {
      return Geometric_Region::INTERFACE;
  }
  ```
- **Fix:** (to write) distance-based snapping in the node hierarchy: label a vertex INTERFACE when the cut point would
  fall within `intersection_tolerance`·h of it along the edge (|ξ ± 1| ≤ tol). Decide per parent entity so it stays
  parallel-consistent. See `studies/xtk_holistic_review/findings/14_cut_meshing_alternatives.md` §4A.
- **Regression test:** `[XTK_PATHOLOGY_HEX_PARALLEL]` must report 0 cells with volume ≤ 0 and a layer that is either
  exact or fully snapped, for every δ.
- **Workaround:** none at the deck level beyond keeping static geometries away from each other by more than ~1e-9 h.

## MB-040 · XTK/GEN: node-hierarchy labels disagree with the cut test for exact/multilinear fields → out-of-bounds cell write
- **Status:** merged (#53), landed on develop as a6c1905d7 (2026-10-22). Before merge: fix-written, branch `fix/mb040-nh-template-hardening` (commit 132e55764), draft PR #53 to develop
  (2026-10-20). XTK/GEN/INT/OPT suites pass and match a develop-level baseline; EXA sensitivity ctests not run.
- **Evidence:** `runs/moris_bugs/MB-040/one_throw_case.log` (template warnings + assert), `one_throw_case_gdb.log`
  (backtrace at the throw), `two_geometries_gdb.log` (follow-on SIGSEGV in tcmalloc from
  `create_edges_from_element_to_node`: heap corrupted by the case before). Before/after sweep table:
  `runs/studies/xtk_pathology/mb040_gate/VERIFICATION_RECORD.md`.
- **Symptom:** single HEX8 [0,1]³, analytic regime (User_Defined_Field, discretization −2), geometries in the order
  B = x + y + z − 1, then A = trilinear interpolant of corners [−1, 1, 1, −1, 1, −1, 0, 0] (corner i at
  (i&1, (i>>1)&1, (i>>2)&1)). Develop prints `WARNING UNDEFINED TEMPLATE ID: 532 / 403 / 430 / 523 / 2147483647`, then
  `Assertion aCellIndex < (moris_index)mIntegrationCells.size() failed` (`XTK/src/cl_XTK_Cut_Integration_Mesh.cpp:1032`).
  The order A, B works. Over the exhaustive two-geometry sweep: analytic 828 and multilinear 808 failing cases of
  104 976 (5 or 6 cut edges on a TET4, or undefined template IDs 34 … 541), plus 404 / 472 cases with an IG cell
  holding a positive and a negative vertex. The stored (discretized) regime is unaffected.
- **Reproducer:** test `[MB-040]` in `projects/XTK/test/xtk/UT_XTK_Interior_Cut_Detection.cpp` (on the fix branch); or
  `MORIS_XTK_PATHOLOGY_ONE="analytic|x+y+z=1|BA|-1;1;1;-1;1;-1;0;0" XTK-test.exe "[XTK_PATHOLOGY_ONE]"`
  (`UT_XTK_Cut_Pathologies.cpp`, branch `docs/xtk-walkthrough`).
- **Root cause:** `Level_Set_Geometry::get_geometric_region` gave a derived node with 2 locators, one of them on the
  interface, the region of the other locator without looking at its value (`GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:99-110`
  on develop). That is exact when the derived value is the linear interpolation of the locators (discretized fields).
  For an exact (analytic) or multilinear value the node's own value can have the other sign, while the node hierarchy
  decides cut edges from the values (`:461-471`): the undefined IDs decode (TET4 edges 0:(0,1) 1:(1,2) 2:(0,2) 3:(0,3)
  4:(1,3) 5:(2,3)) to three cut edges forming one face (532/523 = {2,3,5}, 403/430 = {0,3,4}), which no two-sign
  labelling allows. The damage then came from three latent defects: `load_2d/3d_template` only warned on an unknown
  ID (`XTK/src/cl_XTK_Node_Hierarchy_Template_Library.cpp:100-106`, `:824-830`), `Node_Hierarchy_Template::mNumCells`
  had no initializer (`XTK/src/cl_XTK_Node_Hierarchy_Interface.hpp:40`), and `sort_nodes_2d/3d` had no `else` for an
  unexpected cut-edge count. `commit_new_ig_cells_to_cut_mesh` then sized the cell list from the undefined count and
  called `set_integration_cell` out of range (`XTK/src/cl_XTK_Integration_Mesh_Generator.cpp:794`) before the assert
  (`:797`).
  A first fix (always label from the value for exact/multilinear fields) removed the failures but produced ~25 000
  mixed-sign cases: a vertex created by a later geometry was then labelled against the exact value of an *earlier*
  geometry whose cut is already final. The right rule depends on the order.
- **Fix:** a two-locator node takes its locators' region only if the value interpolates the locators, or if this
  geometry was decomposed before the geometry whose intersection created the node; otherwise the value decides
  (geometries get `set/get_decomposition_index` from the geometry engine; `Intersection_Node::get_interface_decomposition_index`).
  Hardening: `MORIS_ERROR` on undefined template IDs (2D/3D) and on unexpected cut-edge counts; `mNumCells = 0`.
- **Code** (commit 132e55764, `GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp`, trimmed):
  ```diff
  -                // If we only have 2 locators, can use special logic if at least one node is on the interface
  -                if ( tLocators.size() == 2 )
  +                bool tValueInterpolatesLocators = not this->use_multilinear_interpolation() and not mField->has_exact_derived_node_values();
  +
  +                bool                     tDecomposedBeforeNode = false;
  +                const Intersection_Node* tIntersectionNode     = dynamic_cast< const Intersection_Node* >( &tDerivedNode );
  +                if ( tIntersectionNode != nullptr )
                   {
  -                    if ( this->get_geometric_region( tLocators( 0 ).get_index(), tLocators( 0 ).get_global_coordinates() ) == Geometric_Region::INTERFACE )
  +                    tDecomposedBeforeNode = this->get_decomposition_index() < tIntersectionNode->get_interface_decomposition_index();
  +                }
  +
  +                if ( tLocators.size() == 2 and ( tValueInterpolatesLocators or tDecomposedBeforeNode ) )
  +                {
  ...
  +                    else if ( tDecomposedBeforeNode and tFirstRegion == tSecondRegion )
  +                    {
  +                        // Both locators on the same side of a final cut: the node is on that side as well
  +                        return tFirstRegion;
                       }
  ```
- **Regression test:** `[MB-040]` (fails on the hardening alone with `undefined template ID 532`, passes with the fix).
  The exhaustive two-geometry sweep goes to 0 failing and 0 mixed-sign cases in all three regimes.
- **Workaround (until merged):** `discretization_mesh_index = -1` (stored) for analytic geometries combined with other
  geometries, and avoid `use_multilinear_interpolation` with several geometries.

## MB-041 · FEM/VIS: dropped IG elements visualized as zero
- **Status:** resolved by #55 (landed on develop as 0f191448f, 2026-10-22): the volume-drop threshold is now a hard-coded 0,
  so no IG cell is dropped and VIS has no dropped cells to write. PR #54 closed as superseded (2026-10-22); its VIS fix stays on
  `fix/vis-dropped-cells-output` (ce59f590f) and is needed again only if dropping is re-enabled.
  Before: fix-written, branch `fix/vis-dropped-cells-output` (commit ce59f590f), draft PR #54 to develop
  (2026-10-21). XTK/GEN/INT/OPT suites identical to a develop-level baseline; new EXA ctest passes; other EXA ctests
  not run.
- **Evidence:** `runs/moris_bugs/MB-041/ctest_without_fix.log` (240 of 1032 cells off), `ctest_with_fix_final.log`;
  study: `runs/studies/xtk_pathology/patch_test/` (cases `*_visfix`), defect row D-182.
- **Symptom:** immersed linear patch test (same material in every phase, exact linear displacement): the displacement is
  exact, but `ELEMENTAL_AVG` stress output is 0 in many cut cells (21, 240, 329 cells at plane separations 1e-3, 1e-6,
  1e-11; 144 of 324 in a single-plane case). Every zeroed cell has a volume below 1e-6 of its cluster. The VIS mesh
  type does not matter.
- **Reproducer:** EXA ctest `Immersed_Patch_Test-1-procs` (`projects/EXA/structure/linear/Immersed_Patch_Test`, on the
  fix branch).
- **Root cause:** `Cluster::determine_elements_for_residual_and_iqi_computation` (`FEM/INT/src/ELEM/cl_FEM_Cluster.cpp`)
  sets `mComputeResidualAndIQI( iElem ) = false` for the smallest IG elements whose summed relative volume stays below
  `mVolumeError = 1e-6` (`cl_FEM_Cluster.hpp:73`). `Cluster::compute_quantity_of_interest( aFemMeshIndex, aFieldType )`
  applied the same flag to VIS output, and `Set::compute_quantity_of_interest_elemental` initializes the elemental
  field to 0.0 (`FEM/INT/src/CORE/cl_FEM_Set.cpp:3583`; the NaN initialization is commented out at :3582, and the
  "set unused IQI values to not be NAN" test in `Element_Bulk::compute_quantity_of_interest_elemental` compares with
  `== quiet_NaN()`, which is never true, `FEM/INT/src/ELEM/cl_FEM_Element_Bulk.cpp:907`).
- **Code** (commit ce59f590f, `FEM/INT/src/ELEM/cl_FEM_Cluster.cpp`):
  ```diff
  +        bool tEvaluateDroppedElements = aFieldType != vis::Field_Type::GLOBAL;
  +
           // loop over the IG elements
           for ( uint iElem = 0; iElem < mElements.size(); iElem++ )
           {
               // check whether to compute QI
  -            if ( mComputeResidualAndIQI( iElem ) )
  +            if ( mComputeResidualAndIQI( iElem ) or tEvaluateDroppedElements )
  ```
- **Fix:** VIS nodal/elemental fields evaluate every IG element (the solution is defined there); global values keep the
  drop rule.
- **Regression test:** `Immersed_Patch_Test-1-procs`: every output cell's σxx exact to 1e-4 relative; L2 displacement
  error < 1e-8. Fails on develop (240 cells at 0), passes with the fix (max deviation 6.7e-7).
- **Note:** the drop rule itself limits a patch test to ~1e-6 relative accuracy once cells are dropped (by design).

## MB-042 · WRK/HMR: refinement-pattern file written one iteration earlier than the restart test expects
- **Status:** open (2026-10-21). Seen on develop 97c863cb7 and on the 2026-10-12 showcase line (D-151); fails identically
  with and without PR #55.
- **Evidence:** `runs/moris_bugs/MB-042/` — `ctest_restart_full.log` (full ctest output of the test), `signature.txt`,
  `bin_listing.txt` (files the create run left: `HMR_Background_Refinement_Iter_0/10.hdf5`, `ADV_Alg_0_Iter_0..12.hdf5`),
  `HMR_Background_Refinement_Iter_10.hdf5`, `Parameter_Receipt.xml`; defect row D-151.
- **Symptom:** the Create_File case passes (12 GCMMA iterations, `reinitialize_interface_iter = 10`,
  `write_refinement_pattern_file = true`). At iteration 10 the first evaluation returns NaN criteria after
  `xtk::Model::decompose() - Decomposition Failed. Not all child meshes on the same level.` (the remesh signal, cf.
  D-152), the problem reinitializes, iteration 10 is redone, and the pattern is written as `Iter_10`. The Restart case then
  aborts:
  ```
  *** Reason:  Assertion tFile failed.
  ***          open_hdf5_file - could not open HDF5 file HMR_Background_Refinement_Iter_11.hdf5
  *** Where:   ... projects/MRS/IOS/src/HDF5_Tools.hpp:119
  example_test_case.cpp:105: FAILED
  ```
- **Reproducer:** `ctest -R Levelset_Boxbeam_Adaptive_Refinement_Restart-1-procs` in any develop `build_opt` with EXA built
  (~18 s).
- **Code** (develop 97c863cb7):
  - writer, `HMR/src/cl_HMR.cpp:248-254` (unchanged vs the older moris/ line):
    ```cpp
    uint tOptIter = gLogger.get_iteration( "OPT", "Manager", "Perform" );
    ...
    tHDF5.create( "HMR_Background_Refinement_Iter_" + std::to_string( tOptIter ) + ".hdf5" );
    ```
  - reader, `EXA/optimization/LevelSet_Boxbeam_Adaptive_Refinement_Restart/Levelset_Boxbeam_Restart.cpp:531` (regression
    mode; the hardcoded 11 predates #36, which only factored it into `restart_iteration()`):
    ```cpp
    return 11;
    ```
- **Root cause:** not yet located. The writer is unchanged, so the logger's OPT iteration at the reinit (or the iteration
  the reinit fires at) changed between the older moris/ line and develop. Files that differ between the two lines and
  touch this path: `WRK/src/cl_WRK_Reinitialize_Performer.cpp`, `cl_WRK_Workflow_HMR_XTK.cpp`,
  `cl_WRK_perform_remeshing.cpp`, `OPT/src/cl_OPT_Manager.cpp`, `cl_OPT_Algorithm_MMA.cpp`,
  `MRS/IOS/src/cl_GlobalClock.cpp`. Next step: bisect those between the two lines.
- **Fix options:** (a) if develop's numbering is right (the pattern belongs to the iteration that triggered the reinit),
  update the test to restart from 10 (ADV file, `restart_index`, pattern) and re-gold `example_test_case.cpp`; (b) if the
  counter regressed, restore it in the core and keep the test.
- **Workaround:** none needed for production runs. Showcase mode already picks the newest pattern present.

## MB-043 · HMR: short-circuited collective in Database::perform_refinement deadlocks on asymmetric manual flags
- **Status:** open (2026-10-21). Reachable only through the manual flag API (`HMR::flag_element` →
  `Database::flag_element`, `Database::flag_parent`); production refinement uses `flag_elements_on_working_pattern`,
  which does not set the flag, and the flag stored after a refinement is the global result, so production runs stay
  consistent. Unit tests use the manual API (`HMR/test/ut_HMR_L2.cpp:112`, `ut_HMR_BSpline_Mesh_Private.cpp:85`).
- **Evidence:** `runs/moris_bugs/MB-043/np2_stack.txt` (gdb backtraces of both ranks), `np2_probe.log`; defect D-186;
  reproducer run dir `runs/studies/xtk_pathology/2026-10-21_mpi/levels/`.
- **Symptom:** 2 ranks, 4×2×2 box, `HMR::flag_element` on the 4 elements with centre x > 1.5 (all on rank 1; rank 0 flags
  none), then `perform_refinement_based_on_working_pattern( 0 )`: no progress, killed by a 600 s timeout. Stacks:
  ```
  rank 0: PMPI_Allreduce <- Logger::memory_usage <- Logger::sign_in <- hmr::Background_Mesh_Base::perform_refinement
          <- hmr::Database::perform_refinement
  rank 1: PMPI_Wait <- communicate_mats <- Background_Mesh_Base::synchronize_refinement_queue
          <- Background_Mesh_Base::collect_refinement_queue <- hmr::Database::perform_refinement
  ```
- **Reproducer:** `mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI_LEVELS]"` with `Box_Mesh( 1.5 )` flagging through
  `mHMR->flag_element( iElement )` (walk-through branch `docs/xtk-walkthrough`, test file history before the workaround).
- **Code** (develop, `HMR/src/cl_HMR_Database.cpp`):
  ```cpp
  // :894
  bool tFlag = mHaveRefinedAtLeastOneElement;
  // :963, inside the loop over levels
  tFlag = tFlag || mBackgroundMesh->perform_refinement( aActivePattern );
  // :1471, Database::flag_element (rank-local)
  mHaveRefinedAtLeastOneElement = true;
  ```
- **Root cause:** `||` short-circuits, so `Background_Mesh_Base::perform_refinement` (collective: Tracer allreduce,
  `collect_refinement_queue` allgathers) runs only on ranks where `tFlag` is false. Even in serial, once the flag is
  set the per-level refinement is skipped and the later `#MINREF` loop (`:967-970`) does all queued refinement at once.
- **Fix direction:** evaluate the refinement unconditionally, `bool tRefined = mBackgroundMesh->perform_refinement( aActivePattern );
  tFlag = tFlag or tRefined;` (and reduce `mHaveRefinedAtLeastOneElement` globally where it is read). Regression test: 2 ranks,
  `flag_element` on one rank only, refinement completes and both ranks see the same element count.
- **Workaround:** flag elements without the manual API: `put_on_refinement_queue()` + `set_refined_flag( working pattern )`
  on the background element (what `Database::flag_element` does, minus the flag), or `flag_elements_on_working_pattern`.

## MB-044 · XTK: refinement-level check of the intersected cells is per rank
- **Status:** open (2026-10-21).
- **Evidence:** `runs/moris_bugs/MB-044/` (`summary_np1_2_4.txt`, `np{1,2,4}.log`); defect D-187.
- **Symptom:** box 4×2×2 on [0,2]×[0,1]×[0,1], refined one level for x > 1.5; planes x = 0.3 (level-0 cells, rank 0)
  and x = 1.8 (level-1 cells, rank 1). `Model::decompose()` returns false on 1 rank and true on 2 and 4 ranks:
  ```
  np=1  x=0.3 and x=1.8  cut levels { 0 1 }  decomposed 0
  np=2  rank 0 cut levels { 0 }, rank 1 { 1 }  decomposed 1
  ```
  Controls behave on every rank count: a plane across the level transition (both levels on one rank) is refused,
  single-level cuts decompose.
- **Reproducer:** `mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI_LEVELS]"` (walk-through branch `docs/xtk-walkthrough`,
  `projects/XTK/test/xtk/UT_XTK_Cut_Pathologies.cpp`); CHECK fails on 2 and 4 ranks, passes on 1.
- **Code** (develop, `XTK/src/cl_XTK_Integration_Mesh_Generator.cpp:530-547`):
  ```cpp
  bool tFlag = true;
  moris_index tReferenceLevel = 0;
  if ( aMeshGenerationData.mAllIntersectedBgCellInds.size() > 0 )
  {
      tReferenceLevel = aBackgroundMesh->get_mtk_cell( aMeshGenerationData.mAllIntersectedBgCellInds( 0 ) ).get_level();
  }
  for ( size_t iBgCellIndex = 1; iBgCellIndex < aMeshGenerationData.mAllIntersectedBgCellInds.size(); iBgCellIndex++ )
  {
      moris_index tLevel = aBackgroundMesh->get_mtk_cell( aMeshGenerationData.mAllIntersectedBgCellInds( iBgCellIndex ) ).get_level();
      if ( tReferenceLevel != tLevel )
      {
          tFlag = false;
      }
  }
  aCutIntegrationMesh->mSameLevelChildMeshes = all_land( tFlag );
  ```
- **Root cause:** the reference level is rank-local. The aura makes adjacent level changes visible to both ranks
  (the slanted control), so the gap opens only when the coarse and fine cut cells are more than one aura layer apart.
- **Impact:** a parallel run decomposes a configuration the serial check rejects, so serial and parallel runs of the same
  deck take different paths: the refusal is the remesh signal in adaptive optimization runs (D-152). Whether the
  mixed-level decomposition is itself wrong (hanging nodes between child meshes of different levels) is not yet tested.
- **Fix direction:** reduce the levels globally: local min/max over the intersected cells (empty rank: neutral values),
  then `min_all` / `max_all`, and `mSameLevelChildMeshes = ( tGlobalMin == tGlobalMax )`. Regression test: the
  reproducer above must refuse on every rank count.

## MB-045 · XTK/FEM: a near-zero-width neck at a grid node keeps a load-carrying connection
- **Status:** open (found 2026-10-21; limitation of the cut/enrichment discretization, mechanism not isolated).
- **Evidence:** `runs/moris_bugs/MB-045/` (both designs: restart ADVs, `xtk_temp.exo`, GEN and VIS exodus, criteria;
  `ligament_zoom.png`; `mma_trials_excerpt.log`). Study context: `runs/studies/ls_smoothing_km_helmholtz_coarse/diag/helm_L2_s1_stall/`
  (D-188). Builder for `build_bug_evidence.py`: to add.
- **Symptom:** 2D cantilever (80 x 40 quadratic cells, h = 0.025), B-spline level set with the Helmholtz ADV filter
  (r = 0.0345, branch `study/ls-smoothing`). Two designs whose design fields differ by at most 3.7e-4 (no visible change)
  give strain energy **54.7761** and **76.3424**. The difference is at two mirror nodes (0.775, 0.375) and (0.775, 0.625),
  where the ligament between three holes has its neck (a saddle of phi): phi = **+1.31e-6** in the stiff design and
  **-6.98e-6** in the flexible one. With +1.3e-6 the material on the four sides of the node is connected through cut
  cells of area 4.8e-12, 8.5e-12, 1.8e-11 (and 1.8e-8, 8.0e-8): a neck a few 1e-6 wide that carries enough load to
  lower the strain energy by 28 %; with -7.0e-6 no material cell touches the node. The cut is correct geometry (values far
  above `isocontour_tolerance` 1e-12, so not MB-037's snapping); the analysis gives a sub-cell neck finite stiffness.
- **Consequence in optimization:** the volume gradient drives the neck to zero width while the connection keeps the
  stiffness; every MMA trial that crosses phi = 0 at the node is non-conservative by +39 %, `raaupdate` raises raa
  up to 10x per rejected trial and `asymp` relaxes it 0.1x per iteration, so the move shrank from 1e-4 to 1e-8
  (`MMA-TRIAL`: first trial app 0.384 vs true 0.536 in every iteration). The run stalled at strain energy 54.79 for ~25
  iterations, then one of the two mirror necks severed first and the design broke symmetry (final 30.22 vs 29.5-29.8
  for the symmetric runs).
- **Reproducer:** the study deck `runs/moris_bugs/MB-045/LS_Smoothing.cpp` with `XS_COARSEN=0 XS_FILTER=helmholtz
  XS_FILTER_R=0.0345 XS_SEED=1 XS_SEED_R=0.18 XS_MAX_ITS=1 XS_RESTART=<eval9|eval10>/restart.hdf5` (first evaluation =
  the restart design), binary from `moris-ls-adv-filter` (branch `study/ls-smoothing`, 76ff4d738).
- **Hypotheses (not tested):** (1) the material on both sides stays one subphase group through the tiny cells, so it
  shares the enriched B-spline basis, which spans the neck and transmits load regardless of the neck width; (2) the
  ghost penalty on facets next to the node couples the sides. Test: the same two designs with ghost off, and with the
  neck node moved off the grid (phi offset), comparing strain energy.
- **Code:** not yet isolated (candidates: subphase-group connectivity in XTK enrichment; ghost facet selection).
- **Fix:** (to decide) e.g. treat material connected only through cells below a relative-volume threshold as
  disconnected for enrichment, or report sub-cell necks. Any fix changes optimization paths near pinch-off.
- **Workaround:** none at deck level; a merit-gate step-collapse check flags the resulting stall (study runner).

## MB-046 · FEM/INT: 3D principal stress third invariant uses σyz for σzz

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, verified (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: open (found 2026-10-21).
- **Evidence:** `runs/moris_bugs/MB-046/` (stress-path deck, generator, analyzer, bulk exodus of the control case); study
  `runs/studies/xtk_pathology/patch_test/` (`stress_paths.md`, column SP1).
- **Symptom:** immersed linear patch test (uniform exact stress, every other stress output exact to 1e-12): the
  `PRINCIPAL_STRESS` output reads σ1 = 1.31612e-3 everywhere instead of 1.30070e-3 (eigenvalue), +1.2 %, at every node
  and cell, NODAL and ELEMENTAL alike. The code's formula in numpy gives 1.31611756006e-3, matching MORIS to 12
  digits; with σzz in the term it gives the eigenvalue.
- **Reproducer:** `bash runs/studies/xtk_pathology/patch_test/run_stress_paths.sh keepsmall control`, then `analyze_stress_paths.py` (column SP1_N).
- **Root cause:** I3 = det σ = σxx σyy σzz + 2 σyz σxz σxy − σxx σyz² − σyy σxz² − **σzz** σxy². The code has
  `tStressVector( 3 )` (σyz) in the last product. The error vanishes only when σxy = 0 or σzz = σyz.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Stress.cpp:193`; the same term at `cl_FEM_IQI_Max_Stress.cpp:241` and
  `cl_FEM_IQI_Zienkiewicz_Zhu.cpp:216`):
  ```cpp
  real tI3 = tStressVector( 0 ) * tStressVector( 1 ) * tStressVector( 2 ) - tStressVector( 0 ) * std::pow( tStressVector( 3 ), 2.0 ) - tStressVector( 1 ) * std::pow( tStressVector( 4 ), 2.0 ) - tStressVector( 3 ) * std::pow( tStressVector( 5 ), 2.0 ) + 2.0 * tStressVector( 3 ) * tStressVector( 4 ) * tStressVector( 5 );
  //                                                                       the fourth product should use tStressVector( 2 ) ^
  ```
- **Fix:** `tStressVector( 2 )` in all three copies; unit test against `eigvalsh` for a full 3D tensor. Also guard
  `acos` against |R/√Q³| > 1 from round-off and Q = 0 (hydrostatic stress gives 0/0 = NaN).
- **Workaround:** compute principal stresses in post-processing from the six components (`STRESS_VECTOR`, or the
  NORMAL_STRESS and SHEAR_STRESS IQIs).

## MB-047 · FEM/INT: one CM shared by the leader and follower side

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, verified (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: open (found 2026-10-21). The cache and the FI-manager overwrite are confirmed in code; which side's field
  the IWG residual ends up using was not traced.
- **Evidence:** `runs/moris_bugs/MB-047/` (decks: shared, split, shared with 10 Newton iterations;
  `newton_shared_vs_split.txt`).
- **Symptom:** the same CM name in `leader_constitutive_models` and `follower_constitutive_models` (here
  `STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE`; the workspace L-bracket stress decks do it on the
  `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` design|pad interface). With non-polynomial Dirichlet data on the
  2x2x2 patch-test mesh (planes 0.1 apart):
  - 2 Newton iterations: the same displacement as a run with a second, identical CM on the follower side (max diff
    1.5e-17), but the residual at that solution is 4.1e-4 (drop 2.8e-4), vs 3.1e-15 with split CMs.
  - 10 iterations: linear convergence (factor ~0.04, 9 iterations) to a different solution, 5.2e-4 relative in
    displacement and 1.6 % in nodal σxx.
  - The MB-037 thin-layer runs (D-183): drop 1.7e-5 (δ = 1e-11) and 9.4e-5 (δ = 1e-12) shared, vs 7.4e-15 and 1.7e-15
    split, at the same condition number. D-183's Newton stall is this defect, not conditioning.
  - The linear patch test passes either way (σ_L = σ_F in the exact solution).
- **Root cause:** `Model_Initializer_Legacy` hands the same `mConstitutiveModels( tCMIndex )` shared_ptr to both sides
  (`FEM/INT/src/CORE/cl_FEM_Model_Initializer_Legacy.cpp:472-483` for IWGs, `:753-770` for IQIs).
  `IQI::set_field_interpolator_manager` sets each CM's FI manager for the leader and then for the follower
  (`FEM/INT/src/IQI/cl_FEM_IQI.cpp:286-294`; called from `cl_FEM_Set.cpp:1319-1342` and `:1450-1473`), so the last call
  wins. `Constitutive_Model::traction( aNormal )` caches on `mTractionEval` and returns the cached value for any later
  normal (`FEM/INT/src/CM/cl_FEM_Constitutive_Model.cpp:2855-2874`).
  `IWG_Isotropic_Struc_Linear_Interface::compute_residual` sums `w_L * tCMLeader->traction( n ) + w_F * tCMFollower->traction( n )`
  (`cl_FEM_IWG_Isotropic_Struc_Linear_Interface.cpp:122-124`) from one cached value.
- **Code** (`projects/FEM/INT/src/CM/cl_FEM_Constitutive_Model.cpp:2863-2873`):
  ```cpp
  if ( mTractionEval )
  {
      this->eval_traction( aNormal );
      mTractionEval = false;
  }
  return mTraction;
  ```
- **Fix:** (to decide) clone the CM per side in the initializer when the leader and follower names match, or reject a
  shared CM; independently, make `traction()` / `testTraction()` cache per normal (or not cache).
- **Phase-based initializer** (IQI audit C-4): the same sharing of properties and of a phase's MMs/CMs; fixed in PR #57
  (per-phase follower `Phase_User_Info`), A/B bit-identical on four phase-based EXA decks; no deck with a side-dependent
  shared property found, so that case is untested.
- **Workaround:** define a second, identical CM and name it on the follower side. Verified: split-CM runs reach a
  residual drop of 1e-15 in one Newton step in every case, including δ = 1e-12.
- **Open:** effect on the L-bracket stress decks (1 Newton step, so the forward solve takes the Jacobian's solution, but
  dR/dp in the sensitivities uses the inconsistent residual); D-185 (basis agglomeration, "residual and Jacobian
  inconsistent?") used the shared CM too.

## MB-048 · FEM/INT: JUMP_TRACTION sums the two tractions

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, verified (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: open (found 2026-10-21).
- **Evidence:** `runs/moris_bugs/MB-048/faceted_splitcm_control.exo` (JT_F on `dbl_iside_*` = 2.5507392e-3 = 2|σn|
  on every facet, ratio 1 ± 1.4e-13), `faceted_sharedcm_control.exo` (JT_F = 0 identically).
- **Symptom:** in exact equilibrium (uniform stress, no jump) `JUMP_TRACTION` reads 2|σn| with separate CMs and exactly
  0 with one shared CM (MB-047 cache), so it never measures the traction jump.
- **Root cause:** t_L(n) − t_F(−n) = σ_L n + σ_F n. The jump is t_L(n) + t_F(−n) = (σ_L − σ_F) n.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Jump_Traction.cpp:49-51`; the same at `:87-89`):
  ```cpp
  const Matrix< DDRMat >
          tTractionJump =
                  tCMLeader->traction( mNormal ) - tCMFollower->traction( -1.0 * mNormal );
  ```
- **Fix:** `+` instead of `-` (or `traction( mNormal )` on both sides); check the second copy and any EXA gold values
  (`Heated_Sphere` uses it for a thermal flux jump).
- **Workaround:** output `TRACTION` on both single-sided interface sets and add them in post-processing.

## MB-049 · FEM/INT: IQI_Max_Stress rejects every valid component index

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, verified (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: open (found 2026-10-21).
- **Evidence:** `runs/moris_bugs/MB-049/log_signature.txt`, deck `MB037_Patch_Test_MaxStress.cpp`.
- **Symptom:** `MAX_NORMAL_STRESS` with `vectorial_field_index = 0` aborts with "IQI_Max_Stress::compute_QI -
  mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." (`cl_FEM_IQI_Max_Stress.cpp:57`).
  The same check guards MAX_SHEAR and MAX_PRINCIPAL; only MAX_VON_MISES_STRESS skips it.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Max_Stress.cpp:57-58`, again in `compute_QI( real )`):
  ```cpp
  MORIS_ERROR( mIQITypeIndex > 2,
          "IQI_Max_Stress::compute_QI - mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." );
  ```
- **Fix:** `mIQITypeIndex >= 0 && mIQITypeIndex <= 2`, as in `IQI_Stress` (`cl_FEM_IQI_Stress.cpp:46`).
- **Also:** the same inverted check is in `IQI_Zienkiewicz_Zhu::compute_QI( real )` (`cl_FEM_IQI_Zienkiewicz_Zhu.cpp:101`), so a GLOBAL
  ZZ principal/normal/shear indicator aborts too (the elemental `compute_QI( Matrix )` path has the correct check, `:46`).
- **Workaround:** none for those stress types; von Mises works.

## MB-050 · FEM/INT: SP follower constitutive models land in the leader slots

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, NOT exercised by any test (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: fix-written (2026-10-21), branch `fix/stress-iqi-mb046-049` (worktree `moris-stress-iqi-fixes`), not yet built.
- **Evidence:** code reading only; no deck in the workspace uses a follower CM on an SP.
- **Symptom:** an SP with `follower_constitutive_models` gets those CMs in its *leader* CM slots (overwriting the leader
  CMs of the same role) and a null follower slot. Affects `SP_Spalart_Allmaras_Nitsche_Interface`, which reads
  `mFollowerCM( SPALART_ALLMARAS_TURBULENCE )` (`FEM/INT/src/SP/cl_FEM_SP_Spalart_Allmaras_Nitsche_Interface.cpp:97, :158, :230`).
- **Root cause:** the leader/follower loop of `create_stabilization_parameters` calls `set_constitutive_model` without
  the side, whose default is `LEADER` (`FEM/INT/src/SP/cl_FEM_Stabilization_Parameter.hpp:437-440`).
- **Code** (`projects/FEM/INT/src/CORE/cl_FEM_Model_Initializer_Legacy.cpp:356-358`):
  ```cpp
  mStabilizationParameters( iSP )->set_constitutive_model(
          mConstitutiveModels( tCMIndex ),
          tCMNamesPair( iCM )( 1 ) );
  ```
- **Fix:** pass `tIsLeader` (and the follower-side CM instance, MB-047).

## MB-051 · FEM/INT: 2D principal stresses ignore the in-plane shear

- **Status:** merged (#57), landed on develop as 0117bd3ad (2026-10-22). Before merge: fix-written, verified (branch `fix/stress-iqi-mb046-049` on 41f63bcc9, uncommitted; runs/studies/xtk_pathology/patch_test/VERIFICATION.md). Previously: fix-written (2026-10-21), branch `fix/stress-iqi-mb046-049`, not yet built. Unit test `[MB-051]` in
  `projects/FEM/INT/test/UT_FEM_IQI_Stress_3D.cpp`.
- **Symptom:** in 2D, `PRINCIPAL_STRESS` (and the max/ZZ variants) return ( σxx + σyy )/2 ± |σxx − σyy|/2 in plane stress
  (σzz = 0), and use σzz in place of σxy in plane strain / axisymmetry. The in-plane shear never enters.
- **Root cause:** `get_stress_vector` standardizes the stress to ( xx, yy, zz, yz, xz, xy ), so the 2D shear is entry 5
  (`cl_FEM_IQI_Stress.cpp:272-283`), but the 2D radius reads entry 2.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Stress.cpp:165`; the same at `cl_FEM_IQI_Max_Stress.cpp:204-205` and
  `cl_FEM_IQI_Zienkiewicz_Zhu.cpp:179-180`):
  ```cpp
  real tMaxShearStress = std::sqrt( std::pow( ( tStressVector( 0 ) - tStressVector( 1 ) ) / 2.0, 2.0 ) + std::pow( tStressVector( 2 ), 2.0 ) );
  ```
- **Fix:** `tStressVector( 5 )`.
- **Workaround:** compute 2D principal stresses in post-processing from σxx, σyy, σxy.

## MB-052 · SOL/DLA: Belos "Output Frequency" read with the wrong type
- **Status:** open (2026-10-21).
- **Evidence:** `runs/moris_bugs/MB-052/` (`repro/` deck, `run_repro.sh`, `repro/logs/03_running_moris.log`).
- **Symptom:** any deck that sets the Belos parameter `"Output Frequency"` to a non-default value aborts in
  `LinearSolver - Belos - Solve`:
  ```
  terminate called after throwing an instance of 'std::bad_variant_access'
    what():  std::get: wrong index for variant
  ```
- **Reproducer:** `bash runs/moris_bugs/MB-052/run_repro.sh` (patch-test deck, `tLinearSolver = 1`, `"Output Frequency", 10`).
  Without the key the same run converges.
- **Code** (develop, `SOL/DLA/src/cl_DLA_Linear_Solver_Belos.cpp:177-179`; the default is the integer −1,
  `PRM/src/fn_PRM_SOL_Parameters.hpp:654`):
  ```cpp
  if ( mParameterList.get< moris::sint >( "Output Frequency" ) != -1 )
  {
      mMyPl->set( "Output Frequency", mParameterList.get< moris::real >( "Output Frequency" ) );
  ```
- **Fix:** read it as `moris::sint` in both places (Belos expects an int).
- **Workaround:** leave `"Output Frequency"` at its default.

## MB-053 · FEM/INT: `IQI_J_Integral` is wrong in every use
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, verified (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, worktree `moris-iqi-formulas`, 2026-10-22; unit test `[IQI_Formulas]` fails on develop and passes after; patch decks match the closed forms, `runs/studies/iqi_audit/verify_f2/SUMMARY_fixed.txt`; INT/GEN/XTK/OPT = develop + new cases; EXA: same 10 failures as develop, uBend reference updated). Previously: open (2026-10-21). Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1). No workspace or EXA deck uses J_INTEGRAL today.
- **Evidence:** `runs/moris_bugs/MB-053/` (patch deck `Patch_IQI_Audit.cpp`, exodus, `check.txt`).
- **Symptom:** uniform plane-stress patch (u = Gx, E = 1, ν = 0.3), J on the side x = 1 (n = (1, 0)): MORIS (VIS GLOBAL) **+1.978e-7**,
  exact W − t·∂u/∂x = **−4.923e-7**; MORIS equals the code's formula 2W − σxx εxx − σxy ∂ux/∂y to 10 digits.
- **Root cause:** (1) W = σᵀε with engineering shear is 2W (missing ½); (2) `gradx( 1 )( 1, 0 )` is ∂ux/∂y, the term needs ∂uy/∂x
  `( 0, 1 )`; (3) in `compute_QI( real )`, `aWStar *` multiplies only the W term, so both traction terms are added unweighted at every
  integration point; (4) normal hard-coded (1, 0), 2D only.
- **Reproducer:** `XP_DIM=2 moris-run Patch_IQI_Audit.cpp --run-moris --skip-animation` (from `runs/moris_bugs/MB-053/`), read GLOBAL `J`.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_J_Integral.cpp:61-74`):
  ```cpp
  mSet->get_QI()( tQIIndex ) += aWStar * ( trans( mLeaderCM( tElastLinIsoIndex )->flux() ) * mLeaderCM( tElastLinIsoIndex )->strain() ) * tN( 0 )
                              - mLeaderCM( tElastLinIsoIndex )->traction( tN )( 0 ) * mLeaderCM( tElastLinIsoIndex )->strain()( 0 )
                              - mLeaderCM( tElastLinIsoIndex )->traction( tN )( 1 ) * mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX )->gradx( 1 )( 1, 0 );
  ```
- **Fix:** J = ∫ ( ½σ:ε n₁ − tᵢ ∂uᵢ/∂x₁ ) aWStar with the facet normal (`mNormal`) and t = σn, gradx(1)( 0, i ); 2D/3D.
- **Regression test:** the patch above: J = −4.923e-7 (and the same on a rotated side with the facet normal).

## MB-054 · FEM/INT: EIGEN_VALUE criterion always returns the first eigenvalue
- **Status:** merged (#59, develop `1b07bc236`, 2026-10-22). Previously: fix-written, verified (PR #59, branch `fix/iqi-index-handling` @ `5a2d76170`, worktree `moris-iqi-index`, 2026-10-22; unit test `[IQI_Component_Index]` fails before and passes after, INT/GEN/XTK/OPT suites = develop baseline + the new cases; `runs/studies/iqi_audit/verify_f1/`). Previously: open (2026-10-21), confirmed by code reading. Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1).
- **Symptom:** an EIGEN_VALUE IQI declared as in EXA `Column_Buckling` (`vectorial_field_index = i`, no `dof_quantity`) returns
  λ₁ for every i when evaluated as a criterion (`compute_QI( real )`); after that call the reset index also makes the VIS value λ₁.
  The EXA test passes because the deck is not an optimization (`is_optimization_problem = false`) and reads only the VIS path.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Eigen_Value.cpp:55-63`):
  ```cpp
  if ( mQuantityDofType.size() > 1 )
  {
      MORIS_ERROR( mIQITypeIndex != -1, "IQI_Eigen_Value::compute_QI - mIQITypeIndex not set." );
  }
  else
  {
      mIQITypeIndex = 0;
  }
  ```
  `mEigenVectorIndex` (parsed in `set_parameters`) is never used; with the default index −1 the VIS path reads out of bounds.
- **Fix:** select the eigenvalue by `vectorial_field_index` (or `mEigenVectorIndex`) independent of `dof_quantity`; error if unset.
- **Regression test:** Column_Buckling with IQIEIGENVALUE2 as a GEN IQI / OPT criterion (SWEEP, no ADVs): criterion = λ₂ (1.76392 × 10 per the review).

## MB-055 · FEM/INT: `IQI_Stabilization` ignores `vectorial_field_index`
- **Status:** merged (#59, develop `1b07bc236`, 2026-10-22). Previously: fix-written, verified (PR #59, branch `fix/iqi-index-handling` @ `5a2d76170`, worktree `moris-iqi-index`, 2026-10-22; unit test `[IQI_Component_Index]` fails before and passes after, INT/GEN/XTK/OPT suites = develop baseline + the new cases; `runs/studies/iqi_audit/verify_f1/`). Previously: open (2026-10-21), confirmed by code reading. Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1).
- **Symptom:** with `dof_quantity` of 0 or 1 entries (the usual declaration), both `compute_QI` overloads force `mIQITypeIndex = 0`,
  so any value but the first of a multi-valued SP (tauC of `SP_Incompressible_Flow`) silently returns the first (tauM); the
  gradient (`:106`) always takes column 0 too.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Stabilization.cpp:37-44`, same at `:62-69`):
  ```cpp
  if ( mQuantityDofType.size() > 1 )
  {
      MORIS_ERROR( mIQITypeIndex != -1, "IQI_Stabilization::compute_QI - mIQITypeIndex not set." );
  }
  else
  {
      mIQITypeIndex = 0;
  }
  ```
- **Fix:** use `vectorial_field_index` when set (default 0), independent of `dof_quantity`; take the same column in `compute_dQIdu`.

## MB-056 · FEM/INT: bulk power dissipation drops the factor 2 on shear
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, verified (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, worktree `moris-iqi-formulas`, 2026-10-22; unit test `[IQI_Formulas]` fails on develop and passes after; patch decks match the closed forms, `runs/studies/iqi_audit/verify_f2/SUMMARY_fixed.txt`; INT/GEN/XTK/OPT = develop + new cases; EXA: same 10 failures as develop, uBend reference updated). Previously: open (2026-10-21), confirmed by code reading. Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1).
- **Symptom:** `IQI_Power_Dissipation_Bulk` returns 2μ(εxx² + εyy² + εxy²) instead of 2μ ε:ε, so it is not rotation invariant (simple
  shear: μγ²/2; the same flow rotated 45°: μγ²). Poiseuille (h = 1, μ = 1, dp/dx = −12): exact ΔpQ = 12, the code gives 6 + ρ-term.
  It is the objective of EXA `UBend_Pseudo_Time_Continuation_Sensitivity_Test`; its FD unit test passes because dQIdu is the correct
  derivative of the wrong value.
- **Root cause:** `CM_Fluid_Incompressible` stores tensor shear strain εxy = ½(∂u/∂y + ∂v/∂x) (`cl_FEM_CM_Fluid_Incompressible.cpp:328`,
  `:345-347`); σᵀε then counts σxy εxy once. The NS bulk IWG applies `tPre = diag(1, 1, 2)` (`IWG_Incompressible_NS_Velocity_Bulk.cpp:128-140`); the IQI does not.
  Separately (suspect), it adds ρ|v|², which is not a power (Dilgen 2018 uses the Brinkman α|u|²).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Power_Dissipation_Bulk.cpp:80-81`):
  ```cpp
  aQI = tSelect * ( trans( tCMFluid->flux() + tP ) * tCMFluid->strain() +    //
                    tPropDensity->val()( 0 ) * tFIVelocity->val_trans() * tFIVelocity->val() );
  ```
- **Numerical confirmation (2026-10-22):** `runs/moris_bugs/MB-056/` (`Patch_Fluid_IQI.cpp`, the exact Poiseuille fields as L2
  projections, so no solver error enters): POWER_DISSIPATION_BULK = **7.2** at ρ = 1 and **8.4** at ρ = 2, i.e. 6 + 1.2ρ, where the
  viscous dissipation is **12**; the surface POWER_DISSIPATION on inlet + outlet gives **12.0** (correct). Both parts are present: the
  halved shear (6 instead of 12) and the ρ∫|v|² = 1.2ρ term.
- **Fix:** weight the shear rows by 2 (as the IWG) in value and dQIdu. Second term: kept and documented as not part of the dissipation (user decision 2026-10-22).

## MB-057 · FEM/INT: L2_ERROR_ANALYTIC with a field reference compares the solution with itself
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, NOT exercised by any test (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, 2026-10-22; no deck uses the field branch). Previously: open (2026-10-21), confirmed by code reading. Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1). Same as `studies/iqi_audit/FINDINGS.md` B-7.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_L2_Error_Analytic.cpp:50-55` and `:79-92`):
  ```cpp
  tFIField = mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );   // should be mLeaderFieldTypes( 0 )( 0 )
  ...
  real tJumpNorm = norm( tFI->val() - tPropL2Check->val() );   // compute_QI( real ): no field branch, null property
  ```
- **Fix:** fetch the field with `mLeaderFieldTypes( 0 )( 0 )` and add the field branch to `compute_QI( real )`.

## MB-058 · FEM/INT/VIS: Exodus strain energy ignores Thickness
- **Status:** open (2026-10-21). Found by the IQI audit (`runs/studies/iqi_audit/SUMMARY.md`, reviews A/B/C; complements `studies/iqi_audit/` phase 1).
- **Evidence:** `runs/moris_bugs/MB-058/` (same patch as MB-053): GLOBAL STRAIN_ENERGY with `Thickness` = 2 reads **6.2857e-7**, the
  thickness-free value; expected **1.2571e-6**. The criterion overload `compute_QI( real )` multiplies aWStar by the thickness (`:100-104`).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Strain_Energy.cpp:59`, `compute_QI( Matrix )`):
  ```cpp
  aQI = 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain();   // no Thickness factor
  ```
- **Fix:** apply the Thickness property (and bedding as in the criterion path) in `compute_QI( Matrix )`; check the other IQIs that read
  Thickness only in the criterion overload. Related: `studies/iqi_audit` A-9 (2-arg dQIdu ignores Thickness).

## MB-059 · FEM/INT: drag and lift coefficients have the wrong sign
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, verified (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, worktree `moris-iqi-formulas`, 2026-10-22; unit test `[IQI_Formulas]` fails on develop and passes after; patch decks match the closed forms, `runs/studies/iqi_audit/verify_f2/SUMMARY_fixed.txt`; INT/GEN/XTK/OPT = develop + new cases; EXA: same 10 failures as develop, uBend reference updated). Previously: open (2026-10-22). IQI audit suspect (review C-3), confirmed numerically. No deck uses DRAG_COEFF or LIFT_COEFF.
- **Evidence:** `runs/moris_bugs/MB-059/` (`Patch_Fluid_IQI.cpp`, exact Poiseuille u = 6y(1−y), p = 12(1−x), μ = 1, RefPressure 1,
  evaluated on the bottom wall). Force of the fluid on the wall: x **+6** (∫μ ∂u/∂y), y **−6** (−∫p). The IQIs return DRAG **−6.0** and
  LIFT **+6.0**; with p = 0, LIFT 0 and DRAG −6.
- **Root cause:** the leader is the fluid, so `mNormal` is the fluid's outward normal and σn is the traction the body exerts on the fluid.
  The force on the body is −∫σn.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_2D_Drag_Lift_Coefficient.cpp:90, :94`, same at `:151, :155`):
  ```cpp
  aQI = { { tFluidCM->traction( mNormal )( 0 ) / tRefPressure } };
  ```
- **Fix:** negate (force on the body), or document that the value is the reaction on the fluid. Also correct the RefVelocity comment
  ("maximum velocity"): the standard benchmarks normalize by the mean velocity.

## MB-060 · FEM/INT: NS strong residual and crosswind IQIs ignore the index for a one-entry dof_quantity
- **Status:** merged (#59, develop `1b07bc236`, 2026-10-22). Previously: fix-written, verified (PR #59, branch `fix/iqi-index-handling` @ `5a2d76170`, worktree `moris-iqi-index`, 2026-10-22; unit test `[IQI_Component_Index]` fails before and passes after, INT/GEN/XTK/OPT suites = develop baseline + the new cases; `runs/studies/iqi_audit/verify_f1/`). Previously: open (2026-10-22). IQI audit suspect (review C-6), confirmed numerically for STRONG_RESIDUAL_INCOMPRESSIBLE_NS;
  SP_CROSSWIND_INCOMPRESSIBLE_NS has the same code (not run). The only deck blocks (uBend, `:2276-2300`) are commented out; they set no
  `dof_quantity`, so all three would return the continuity residual.
- **Evidence:** `runs/moris_bugs/MB-060/` (`Patch_Fluid_IQI.cpp`, velocity exact, p = 0, so the x-momentum strong residual is 12 per
  unit area): `dof_quantity "VX"`, index 1 → **0** (div v, identical to index 0); `"VX,VY"`, index 1 → **12.0**.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Strong_Residual_Incompressible_NS.cpp:51-59`; `cl_FEM_IQI_SP_Crosswind_Incompressible_NS.cpp:57-65`):
  ```cpp
  if ( mQuantityDofType.size() > 1 ) { MORIS_ERROR( mIQITypeIndex != -1, ... ); }
  else { mIQITypeIndex = 0; }
  ```
- **Fix:** reset only an unset index (−1), and error on an index out of range. Same class as MB-054 and MB-055.

## MB-061 · FEM/INT: contact pressure mixes stress measures and configurations
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, verified (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, worktree `moris-iqi-formulas`, 2026-10-22; unit test `[IQI_Formulas]` fails on develop and passes after; patch decks match the closed forms, `runs/studies/iqi_audit/verify_f2/SUMMARY_fixed.txt`; INT/GEN/XTK/OPT = develop + new cases; EXA: same 10 failures as develop, uBend reference updated). Previously: open (2026-10-22). IQI audit suspect (review A MED-1), confirmed numerically. Active users: the three linear contact EXA decks
  (Parabolic_Indenter_Linear, Struc_Contact_Patch_Test, Struc_Hertz_Contact) use CONTACT_PRESSURE_REFERENCE with a linear CM; the only
  CONTACT_PRESSURE_CURRENT block (Parabolic_Indenter_NonLinear) is commented out.
- **Evidence:** `runs/moris_bugs/MB-061/` (`Patch_Contact_Pressure.cpp`, u = (F − I)X as an L2 projection, Saint Venant–Kirchhoff,
  plane strain, E = 1, ν = 0.3, side X = 1, closed forms in `check_contact.py`):

  | F | CURRENT (SVK) | n·σn (Cauchy) | N·P N (nominal) | REFERENCE (linear CM) | N·σ_lin N |
  |---|---|---|---|---|---|
  | diag(0.9, 1.1) | −6.0577e-2 = n·Pn | −5.5070e-2 | −6.0577e-2 | −7.6923e-2 | −7.6923e-2 |
  | [[0.9, 0.1], [0.05, 1.1]] | −6.0223e-2 = n·Pn | −5.0655e-2 | −5.0889e-2 | −8.6066e-2 = n·σ_lin n | −7.6923e-2 |

  CURRENT is 1.10× and 1.19× the Cauchy pressure. REFERENCE uses the deformed normal, so even with the linear CM it is 12% off at this
  (large) gradient; the error is O(|∇u|) in small-strain contact. CONTACT_PRESSURE_REFERENCE with the SVK CM aborts:
  `CM_Struc_Nonlinear_Isotropic::traction - Unknown stress type` (`cl_FEM_CM_Struc_Nonlinear_Isotropic.cpp:1763`, DEFAULT not handled).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Contact_Pressure.cpp:57-69`):
  ```cpp
  Matrix< DDRMat > tNormal = tLeaderGeometry->get_normal_current( tLeaderDofs );   // comment says "undeformed"
  ...
  aQI = { dot( tCMLeader->traction( tNormal ), tNormal ) };                          // REFERENCE
  ...
  aQI = { dot( tCMLeader->traction( tNormal, CM_Function_Type::PK1 ), tNormal ) };  // CURRENT
  ```
- **Fix:** CURRENT = n·σn with `traction( n, CAUCHY )`; REFERENCE = N·(P N) with the IG normal `mNormal` (for a linear CM, N·σN).

## MB-062 · FEM/INT: analytic geometric sensitivities ignore spatially varying material properties
- **Status:** merged (#61, develop `5eb5acf63`, 2026-10-22). Previously: fix-written, verified (PR #61, branch `fix/mb062-analytic-property-variation` @ `564fc0bb9`, worktree `moris-geo-sens-props`, 2026-10-22; `[Analytic_Geo_Sens_Varying]` fails before / passes after; SWEEP se_vary_an 2.71e-2 → 6.73e-9 = FD, `runs/studies/iqi_audit/verify_f3/sweep_final/check.txt`; rule: analytic dP/dx when the property has space derivatives, constant when translation-invariant, FD fallback otherwise). Previously: open (2026-10-22). IQI audit suspect (review A MED-2), confirmed numerically. All decks that enable
  `is_analytical_sensitivity` (L-bracket stress 2D/3D, Sphere3D_Precision) use a constant E, so no result is affected today. A SIMP- or
  density-coupled level-set deck, or any E(x), on the analytic path gets a silently wrong gradient.
- **Evidence:** `runs/moris_bugs/MB-062/` (`Sweep_SE_Varying_E.cpp`: 8×8 linear, circular hole radius ADV, clamped left side, body load,
  ghost; SWEEP central FD ε = 1e-5):

  | E | FD geometry path | analytic path (PR #35) |
  |---|---|---|
  | 1 | 1.2e-8 | 1.2e-8 |
  | 1 + 2x | 6.7e-9 | **2.71e-2** (adjoint −7.1299e-2 vs FD −6.9421e-2) |

- **Root cause:** the analytic dR/dX and dQ/dX hold the properties fixed in space (stated in the PR #35 commit message), but
  `supports_analytical_*_geometry` checks only the CM type, so a property with a spatial (or field/dv) dependence is never routed to the
  FD fallback. The missing term is ½ εᵀ (∂C/∂E) ε · ∂E/∂x · δx at each moved integration point (and its residual counterpart).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Strain_Energy.cpp:279-286`; `IWG/cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp:228-238`):
  ```cpp
  const CM_Struc_Linear *tCMStrucLinear = dynamic_cast< const CM_Struc_Linear * >( ... );
  return tCMStrucLinear != nullptr && tCMStrucLinear->has_spatial_strain_derivatives();   // no check of property dependence
  ```
- **Fix:** return false (FD fallback) when any CM, bedding or thickness property is space-, field- or dv-dependent; or add the ∂C/∂x term
  through `Property::dPropdx`. The same gap applies to IQI_Volume with a varying density (studies/iqi_audit B-5).

## MB-063 · FEM/INT: heat-method penalty returns NaN far from the interface
- **Status:** merged (#60, develop `e8ab19bd0`, 2026-10-22). Previously: fix-written, verified (PR #60, branch `fix/iqi-formulas` @ `43b41ed9e`, worktree `moris-iqi-formulas`, 2026-10-22; unit test `[IQI_Formulas]` fails on develop and passes after; patch decks match the closed forms, `runs/studies/iqi_audit/verify_f2/SUMMARY_fixed.txt`; INT/GEN/XTK/OPT = develop + new cases; EXA: same 10 failures as develop, uBend reference updated). Previously: open (2026-10-22). IQI audit suspect (review B MED-4), confirmed numerically. Not reachable in the current decks: with
  band = 3h it needs a heat-method distance > ~400 h on the penalized side; the largest meshes in use reach ~134 h (MBB 120×60) and
  ~283 h (L-bracket diagonal).
- **Evidence:** `runs/moris_bugs/MB-063/` (`Patch_Heat_Penalty.cpp`: PHID = constant as an L2 projection, B = g = γ = 1, weights 1,
  zero references, expected QI = 1): PHID = −300 / −350 (exponent 600 / 700) → **1.0**; PHID = −360 / −400 (exponent 720 / 800) → **NaN**,
  no error in release.
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Heat_Method_Penalty.cpp:96-104`):
  ```cpp
  real tVal = std::exp( -2.0 * mPhiGradient * mLevelSetSign * tFI->val()( 0 ) / mPhiBound );   // inf beyond 709.8
  MORIS_ASSERT( std::isfinite( tVal ), ... );                                                 // debug only
  Matrix< DDRMat > tPhiTildeDx = ( 4.0 * mPhiGradient * tVal ) / std::pow( 1.0 + tVal, 2 ) * ...;   // inf / inf = NaN
  ```
- **Fix:** the algebraically identical tanh form: φ̃ = B tanh(y), ∇φ̃ = g s (1 − tanh² y) ∇PHID with y = g s PHID / B.

## MB-064 · FEM/INT: IQI_Dof and IQI_Jump_Dof mishandle the component of a scalar field
- **Status:** merged (#59, develop `1b07bc236`, 2026-10-22). Previously: fix-written, verified (PR #59, branch `fix/iqi-index-handling`, worktree `moris-iqi-index`, 2026-10-22; `[IQI_Component_Index]`).
  Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` B-1, B-3; `runs/studies/iqi_audit/SUMMARY.md` A-MED-4). Origin: the index
  block of `2316396cc` (lisenoel, 2020), a slip (`studies/iqi_audit/ORIGINS.md`). No workspace deck meets the trigger.
- **Symptom:** (B-1) DOF with a scalar `dof_quantity` (TEMP, PHID) as an OPT criterion with analytic sensitivities aborts:
  `IQI_Dof::compute_dQIdu - derivative cannot be computed as mIQITypeIndex not set.` (B-3) A spatial or time derivative of a scalar
  on the NODAL/ELEMENTAL output path reads column −1 (`Mat::operator(): index out of bounds` in the test, garbage otherwise). JUMP_DOF
  has the same −1 read for a scalar derivative.
- **Evidence:** `runs/moris_bugs/MB-064/` (unit test and its output without the fix).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Dof.cpp:208, 223` and `:89-96, 116, 124` on develop):
  ```cpp
  else if ( mQuantityDofType.size() > 1 && mIQITypeIndex != -1 ) { /* selection */ }
  else { MORIS_ERROR( false, "... mIQITypeIndex not set." ); }      // a scalar always lands here
  aMat = { tSpatialGradient( mSpatialDerivativeDirection, mIQITypeIndex ) };   // -1 on the VIS path
  ```
- **Fix:** `IQI_Dof::get_dof_component_index` (0 for a scalar field) with `IQI::get_component_index`; JUMP_DOF likewise.

## MB-065 · FEM/INT: DOF IQI with dof_quantity "UY" returns UX
- **Status:** merged (#59, develop `1b07bc236`, 2026-10-22). Previously: fix-written, verified (PR #59, branch `fix/iqi-index-handling`, worktree `moris-iqi-index`, 2026-10-22; `[IQI_Component_Index]`).
  Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` B-2). Origin: `2316396cc` (lisenoel, 2020), the scalar reset applied to a
  member of a vector group (intended for scalars; `ORIGINS.md`). Decks writing `dof_quantity "UY"` without a group: not yet counted.
- **Symptom:** `dof_quantity "UY"` resolves to the {UX,UY} field interpolator; `compute_QI( real )` resets the index to 0, so the
  criterion is ∫UX, and the analytic dQIdu is zero (the dof-type comparison is against UY, never the group key UX). On the output path
  the whole vector is returned and NODAL stores component 0 (UX).
- **Evidence:** `runs/moris_bugs/MB-065/` (same unit test; without the fix the "single member" section fails).
- **Code** (`projects/FEM/INT/src/IQI/cl_FEM_IQI_Dof.cpp:70-77`; `IP/cl_FEM_Field_Interpolator_Manager.cpp:350`):
  ```cpp
  if ( mQuantityDofType.size() > 1 ) { MORIS_ERROR( mIQITypeIndex != -1, ... ); }
  else { mIQITypeIndex = 0; }                                  // "UY" -> component 0 = UX
  if ( tDofType( 0 ) == mQuantityDofType( 0 ) )                 // UX group vs UY: never true -> dQIdu = 0
  ```
- **Fix:** use the position of the dof type in its group, and compare against the group key `tFI->get_dof_type()( 0 )`.

## MB-066 · FEM/INT: 2D geometric stiffness of the linear CM uses σzz as the shear stress in plane strain
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` A-12). Origin:
  `6ebde3640`, kkmaute, 2024 ("added functionality to compute geometric stiffness"), a slip (`ORIGINS.md`). The body is adapted from
  the SVK version, which converts its 3-entry vector with `voigt_to_full_sym_stress`; the assert message still names the SVK function.
  Exposure: 0 decks meet the trigger. The only users of `STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS` with the linear CM are EXA
  `Column_Buckling.cpp:209-210` and `Column_Buckling_Immersed.cpp:313-314`, and both are PLANE_STRESS, where the code is correct.
- **Symptom:** in 2D PLANE_STRAIN or AXISYMMETRIC, the linear CM's flux is (σxx, σyy, σzz|σθθ, σxy) (`CM/cl_FEM_CM_Struc_Linear.cpp:166-210`,
  the strain layout at `:382-386`), but the geometric stiffness builds the 2×2 stress from entry 2. Pure shear (σxx = σyy = 0, σxy = τ)
  gives σzz = 0 in place of τ. For σxx = σyy = s with ν = 0.3, it injects σzz = 2νs as a fake shear. The geometric stiffness, and with it
  the linear buckling eigenvalues, are wrong. Used by `IWG/cl_FEM_IWG_Isotropic_Struc_Nonlinear_Geometric_Stiffness.cpp:127`
  (linear-CM branch).
- **Code** (`CM/cl_FEM_CM_Struc_Linear.cpp:1110-1117`, origin/develop):
  ```cpp
        if ( mSpaceDim == 2 )
        {
            tStressMatrix.set_size( mSpaceDim, mSpaceDim );
            tStressMatrix( 0, 0 ) = mFlux( 0 );
            tStressMatrix( 1, 1 ) = mFlux( 1 );
            tStressMatrix( 0, 1 ) = mFlux( 2 );
            tStressMatrix( 1, 0 ) = mFlux( 2 );
        }
  ```
  Compare `:627` in the same file, which handles the layout: `const uint tShearRow = ( mPlaneType == Model_Type::PLANE_STRAIN ) ? 3 : 2;`
- **Fix:** take the shear from entry 3 when `mFlux.numel() == 4` (plane strain, axisymmetric), entry 2 when it is 3. Or reuse a
  shared Voigt→tensor helper (FINDINGS fix family 2, with MB-051).
- **Regression test:** a plane-strain element with u = A x (pure shear and equal-biaxial cases). Compare `GeometricStiffness` with
  ∫ ∇Nᵀ σ ∇N built from the full 2×2 σ. Also run it in plane stress as a control (passes today).

## MB-067 · FEM/INT: 2D STRESS_VECTOR adds a 6-component vector into a 3- or 4-component QI
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` A-1). Origin:
  `60ff03463`, momo, 2021 ("Modified stress IQI to be able to return a vector"), a slip: in the same commit, `get_matrix_dim`
  returns the CM layout and `get_stress_vector` returns the standardized 6-vector (`ORIGINS.md`). Exposure: 0 decks use
  `STRESS_VECTOR` or `STRESS_VECTOR_CAUCHY`.
- **Symptom:** a 2D GLOBAL or criterion `STRESS_VECTOR`:
  - `compute_QI( real )` adds the standardized (xx, yy, zz, yz, xz, xy) vector into `mSet->get_QI()( tQIIndex )`. That QI is sized
    by `get_matrix_dim`: (3, 1) for plane stress, (4, 1) for plane strain (`CORE/cl_FEM_Model.cpp:1117-1120`, copied per set in
    `CORE/cl_FEM_Set.cpp:2603-2607`).
  - A debug build throws an Armadillo size mismatch.
  - An opt build defines `ARMA_NO_DEBUG` (`CMakeLists.txt:783-785`), so the `+=` is unchecked and the result is undefined.
  - AXISYMMETRIC falls to the `default` of `get_matrix_dim` and aborts. 3D is correct.
- **Code** (`IQI/cl_FEM_IQI_Stress.cpp:117-126` and `:256-259, :326-340`, origin/develop):
  ```cpp
            case Stress_Type::STRESS_VECTOR:
                this->get_stress_vector( tStressVector );
                break;
  ...
        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tStressVector );
  ...
    void IQI_Stress::get_stress_vector( Matrix< DDRMat >& aStressVector )
    {
        // create stress vector
        aStressVector.set_size( 6, 1, 0.0 );
  ...
                    case Model_Type::PLANE_STRESS:
                    {
                        return std::make_pair( 3, 1 );
                        break;
                    }

                    // 2D plane strain and axisymmetric
                    case Model_Type::PLANE_STRAIN:
                    {
                        return std::make_pair( 4, 1 );
  ```
- **Fix:** make the two agree. Either `get_matrix_dim` returns (6, 1) in every dimension, or STRESS_VECTOR assembles the CM-layout flux.
  Add the AXISYMMETRIC case (4 entries). Fix after PR #57 lands (same files).
  Note (2026-10-24): PR #63 keeps this behaviour; after it lands the vector comes from `standardize_stress_vector`
  (`CORE/fn_FEM_Stress_Measures`) via `IQI_Stress::get_stress_vector`. Confirmed again by the PR #63 code review.
- **Regression test:** plane stress, plane strain and 3D, each a single element with u = A x and a GLOBAL `STRESS_VECTOR`. Expect the
  number of `get_matrix_dim()` rows to equal the length of the assembled vector, and the integrated QI to equal σ · area in the
  chosen layout. Fails today in 2D (debug abort).

## MB-068 · FEM/INT: 2D SHEAR_STRESS index is not ignored, as documented: 0 and 1 return 0
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` A-2). Origin:
  `012616842` and `756844958`, wunsch, 2020 (the stress-IQI overhaul), a slip. The header written in the same commit says the input
  is ignored in 2D, and the code maps it (`ORIGINS.md`). Exposure: 0 decks meet the trigger. 30 decks use SHEAR_STRESS, and all pass
  `vectorial_field_index` 2, several with a comment on this trap:
  - `runs/studies/xtk_pathology/2026-10-21_connectivity/crosstalk_fem/deck/Crosstalk_Gap.cpp:504` and its 28 `cases/*/src` copies;
  - EXA `Struc_Hertz_Contact.cpp:656-665`.
- **Symptom:** in 2D the standardized vector has σxy at entry 5, and entries 3 and 4 are always 0 (`IQI/cl_FEM_IQI_Stress.cpp:270-275`).
  - SHEAR_STRESS at index 0 (the header's reading) returns σyz = 0, index 1 returns σxz = 0, and only index 2 returns σxy.
  - The index guard (`:42-49`, `:87-94`) accepts 0..2.
  - `MAX_SHEAR_STRESS` has the same mapping (`IQI/cl_FEM_IQI_Max_Stress.cpp:300`). On develop it is masked, because MB-049's
    inverted check aborts every index 0..2. It becomes reachable once PR #57 lands.
- **Code** (`IQI/cl_FEM_IQI_Stress.hpp:117-120`, `IQI/cl_FEM_IQI_Stress.cpp:250-251`, `IQI/cl_FEM_IQI_Max_Stress.cpp:299-300`;
  the same doc at `IQI/cl_FEM_IQI_Max_Stress.hpp:181-184`):
  ```cpp
         * evaluate the i-th shear stress using the stress vector provided by the constitutive model
         * @param[ in ] aStressIndex number of the principal stress to get,
         *               i.e. 1,2,3 for the yz-, xz-, and xy-shear stresses, respectively
         *               in 2D the input is ignored, as there's only one shear stress
  ...
        // pick stress value
        aStressVector.set_size( 1, 1, tStressVector( aStressIndex + 2 ) );
  ...
        // pick stress value
        return tStressVector( aStressIndex + 2 );
  ```
- **Fix:** in 2D, return entry 5 for any index (the documented contract), or reject indices 0 and 1 with a `MORIS_ERROR`. Either way,
  make the header match. The same applies to the copy in `IQI_Zienkiewicz_Zhu`. Fix after PR #57 lands (same files).
  Note (2026-10-24): after PR #63 the only copy is `eval_shear_stress` in `CORE/fn_FEM_Stress_Measures.cpp` (all stress IQIs and
  `IWG_Struc_Stress` call it); fix and document it there.
- **Regression test:** a 2D plane-stress element with γxy ≠ 0. Assert that SHEAR_STRESS at indices 0, 1 and 2 all equal σxy, or that 0
  and 1 raise an error. Repeat for MAX_SHEAR_STRESS after the MB-049 fix.

## MB-069 · FEM/INT: Cauchy stress IQIs drop σzz for a 2D nonlinear CM in plane strain
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` A-3; the trigger
  is corrected here). Origin: `014645d7b`, noel, 2022 ("Added Cauchy stress to stress IQI"), unclear: the 3-entry "plane stress"
  mapping of `012616842` was written for linear CMs (`ORIGINS.md`). Exposure: 0 decks meet the trigger. The only `*_CAUCHY` user is
  EXA `Parabolic_Indenter_NonLinear.cpp:742-777`. It uses VON_MISES_STRESS_CAUCHY and NORMAL_STRESS_CAUCHY at indices 0 and 1, as
  ELEMENTAL_AVG output, with SVK and no `model_type`, so the deck default PLANE_STRESS applies.
- **Correction to FINDINGS:** `mPlaneType = PLANE_STRAIN` (`CM/cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp:112`) is only the class default.
  A deck CM always gets `set_model_type`:
  - the initializers pass `model_type` unless it is `END_MODEL_TYPE` (`CORE/cl_FEM_Model_Initializer_Phasebased.cpp:193, :202-205`);
  - an unset enum parameter is index 0 (`MRS/CNT/src/cl_Parameter.cpp:21-22`), which is `PLANE_STRESS` (`PRM/ENM/src/cl_FEM_Enums.hpp:251-252`).

  The trigger is therefore an explicit `model_type = PLANE_STRAIN` (or a CM built in C++ without `set_model_type`).
- **Symptom:**
  - The 2D nonlinear CMs return a 3-entry Voigt vector in both plane modes: `mVoigtSymMap = { { 0, 2 }, { 2, 1 } }`, with the size
    taken from `mConst`, which is 3×3 for both (`CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Bonet.cpp:79, :102`).
  - The Cauchy stress is F S Fᵀ / J from the 2×2 F, so σzz is never formed.
  - `IQI_Stress::get_stress_vector` maps 3 entries as "2D plane stress", with σzz = 0.
  - In plane strain (F33 = 1), σzz ≠ 0: for Bonet NH σzz = λ ln J / J, and for SVK σzz = λ tr(E) / J. So VON_MISES_STRESS_CAUCHY and
    PRINCIPAL_STRESS_CAUCHY are wrong, and NORMAL_STRESS_CAUCHY at index 2 returns 0.
  - Example: Bonet, λ = μ = 1, F = diag(1.1, 1.0) gives a missing σzz = λ ln 1.1 / 1.1 ≈ 0.0866.
  - Linear CMs are not affected, because their plane-strain flux has 4 entries (`:278-284`).
  - Related CM-level note (not part of this entry): the 2D Bonet PK2 (`…_Bonet.cpp:199-202`) does not depend on `mPlaneType`. Its
    "plane stress" is plane-strain kinematics with σzz discarded.
- **Code** (`IQI/cl_FEM_IQI_Stress.cpp:261-276`; `CM/cl_FEM_CM_Struc_Nonlinear_Isotropic.cpp:219-221`;
  `CM/cl_FEM_CM_Struc_Nonlinear_Isotropic_Compressible_Neo_Hookean_Bonet.cpp:230-231`, origin/develop):
  ```cpp
        // get stress vector from Constitutive model
        const Matrix< DDRMat >& tCMStress =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux( mFluxType );

        // pull apart stress vector into components
        uint tNumStressComponents = tCMStress.length();
        switch ( tNumStressComponents )
        {
            // 2D plane stress
            case 3:
            {
                aStressVector( 0 ) = tCMStress( 0 );
                aStressVector( 1 ) = tCMStress( 1 );
                aStressVector( 5 ) = tCMStress( 2 );
                break;
            }
  ...
                // set Voigt notation map
                mVoigtNonSymMap = { { 0, 1 }, { 2, 3 } };
                mVoigtSymMap    = { { 0, 2 }, { 2, 1 } };
  ...
        // cast full Cauchy stress to voigt notation
        this->full_to_voigt_sym_stress( this->deformation_gradient() * t2PKStressFull * trans( this->deformation_gradient() ) / this->volume_change_jacobian(), mCauchyStress );
  ```
- **Fix:** branch the IQI mapping on `get_plane_type()`, not on the vector length, and have the nonlinear CMs supply σzz in plane strain.
  If they do not supply it, raise a `MORIS_ERROR` for a Cauchy stress IQI with a nonlinear CM in PLANE_STRAIN. Fix after PR #57 lands
  (same files).
- **Regression test:** a 2D Bonet NH element with `set_model_type( PLANE_STRAIN )`, homogeneous F = diag(1.1, 1.0), λ = μ = 1. Compare
  VON_MISES_STRESS_CAUCHY with the closed form from the full 3×3 Cauchy tensor, including σzz ≈ 0.0866. It fails today. Plane stress
  serves as the control.

## MB-070 · FEM/INT: max-stress and MAX_DOF integrands return NaN for a fractional exponent below the shift
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` A-4, with B-14).
  Origin: `fe7ef6598` and `10d68929a`, wunsch, 2020 (first max-von-Mises IQI). The default was kept by `bb5bf35f4`, wunsch, 2020
  ("added ability to change constant shift"), intended (implicit). The MAX_DOF `mShift = 1.0` (`:35`) is in the same row of
  `ORIGINS.md` and was not traced separately.
- **Exposure:** 0 decks reach the NaN.
  - `MAX_*_STRESS`: 348 workspace files, and every one binds a `Shift` property. No EXA deck uses `MAX_*_STRESS`.
  - `MAX_DOF`: 70 files. 20 rely on the default shift 1:
    - 18 EXA decks: Heat_Method, SA_Cut_PCBar_Transient, the Shape_Sensitivity_Circle/Two_Material sweeps, Thermal_*_Flow_About_Sphere,
      Solver_Examples, Two_Bar_Truss, Single_Phase_Hollow_Cylinder_Static and Stefans_Problem_conform;
    - 2 copies in `runs/studies/xtk_pathology/stress_iqi_ab/exa/`.

    All 20 use integer exponents (1, 2, 8, 30), so they get the (u/ref − 1)^p form without a NaN.
- **Symptom:**
  - The integrand is pow( q/ref − shift, p ), with q the stress measure (MAX_*_STRESS) or the dof value (MAX_DOF). Shift defaults to
    1.0, and the MAX_*_STRESS header does not mention it.
  - `std::pow` of a negative base with a non-integer exponent is NaN. So with p = 8.5 and the shift omitted, every point with q < ref
    gives NaN, and the criterion is NaN.
  - The geometric derivative (`:422`) and the MAX_DOF derivatives (`:205, :269`) use the same base with p − 1.
  - MAX_DOF with `sign = -1` clamps the base to ≤ 0, so a fractional p gives NaN wherever the integrand is active. With `sign = 1` and
    p < 1, `pow( 0, p − 1 )` is inf.
  - With an even integer p and the default shift, the integrand (q/ref − 1)^p is minimal at q = ref and penalizes low stress, which
    is not a p-norm.
  - MAX_NORMAL, MAX_SHEAR and MAX_PRINCIPAL_STRESS with shift 0 and a negative component give NaN for a fractional p.
  - The only guards are `MORIS_ASSERT( isfinite … )` in MAX_DOF (debug only). `ReferenceValue = 0` is not checked.
- **Code** (`IQI/cl_FEM_IQI_Max_Stress.cpp:98-104`, `:166-167`, `:422`; `IQI/cl_FEM_IQI_Max_Dof.cpp:35`, `:82`, `:104`, `:107-108`,
  origin/develop):
  ```cpp
        real tShift = 1.0;

        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr )
            tShift = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );

        // evaluate the QI
        aQI = { { std::pow( 1.0 / tRefValue * tStressValue - tShift, tExponent ) } };
  ...
        mSet->get_QI()( tQIIndex ) += aWStar * ( std::pow( 1.0 / tRefValue * tStressValue - tShift, tExponent ) );
  ...
        const real tdQIdVM = tExponent * std::pow( tVonMises / tRefValue - tShift, tExponent - 1.0 ) / tRefValue;
  ...
            mShift = 1.0;
  ...
        real tArgument = tFIMaxDof->val()( mIQITypeIndex ) / mRefValue - mShift;
  ...
            aQI = { { std::pow( tArgument, mExponent ) } };
  ...
        MORIS_ASSERT( isfinite( aQI ),
                "IQI_Max_Dof::compute_QI - QI is nan, exiting!" );
  ```
- **Fix:** document `Shift` and its default in both headers, without changing the default: that is a behaviour change for every
  MAX_DOF deck above, which is a decision (`ORIGINS.md` item 5). Raise a `MORIS_ERROR` at initialization for a non-integer exponent
  when the base can be negative (`sign` ≠ 1, or MAX_* without a shift ≥ the possible minimum), or use |base|^p · sign(base). Make the
  finiteness checks `MORIS_ERROR`, and check ref ≠ 0. Fix after PR #57 lands (same files, `IQI_Max_Stress`).
  Note (2026-10-24): confirmed again by the PR #63 code review (now also reaches the analytic dQIdu and dQIdx). After #63 the
  MAX_*_STRESS integrand and its derivative are `IQI_Max_Stress::eval_integrand` / `eval_dintegrand_dstress`; the KS form added there
  is finite for any shift.
- **Regression test:** on a harness element with σ_vm < ref, assert that MAX_VON_MISES_STRESS is finite and equals (σ/ref)^p when the
  shift is 0, and that a fractional p with the default shift is an error, not NaN. Repeat with a compressive σxx for MAX_NORMAL_STRESS.
  For MAX_DOF: p = 2.5 and u = 0.5 with the default parameters, expect an error at init or a finite value. Add `sign = -1` with p = 2.5.

## MB-071 · FEM/INT: IQI normalization does nothing, and "design" returns nproc in parallel
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` C-1). Origin: the
  framework is `c769494bd`, Adam, 2020 ("IQI normalization framework"); the divisions by `mReferenceValue` were deleted with the
  functions that held them in `0cc8dfbb1` and `afe975b6b`, lisenoel, 2021. A slip (`ORIGINS.md`).
- **Exposure:** 0 active decks. The only `"design"` is commented out in EXA `Mach_Leading_Edge.cpp:1311`
  ("<-- what does this do?"). The parameter is advertised in `PRM/src/fn_PRM_FEM_Parameters.hpp:249` and exposed by pymoris
  `FEMBlock.normalization` (`pymoris/src/pymoris/blocks/fem.py:147`, default `"none"`). `test/UT_FEM_IQI_Strain_Energy.cpp:278` asserts
  the serial "design" value 1.0.
- **Symptom:**
  - `normalization = "2.5"`: silently ignored, and the criterion stays raw.
  - `"design"`, serial: the criterion is 1.0 at iteration 0 and raw from iteration 1 on, and gradients are never scaled. The objective
    is discontinuous and inconsistent with its gradient.
  - `"design"`, N ranks: `normalize_IQIs` runs on the rank-local `mGlobalIQIVal` (`Equation_Model::compute_IQIs`), including ranks with
    no elements, and the workflow then sums over ranks, so the iteration-0 criterion is N.
  - The reference is set on `mIQIs( tRequestedIQIIndex )`, but the name lookup above it computed `tIQIIndex`.
- **Code** (`IQI/cl_FEM_IQI.cpp:237-245`; `CORE/cl_FEM_Model.cpp:974-978`; `projects/FEM/MSI/src/cl_MSI_Equation_Model.cpp:110-114`;
  `projects/WRK/src/cl_WRK_Workflow_HMR_XTK.cpp:902-905`, origin/develop; `mReferenceValue` has no reader in `projects/`):
  ```cpp
    void
    IQI::set_reference_value( real aReferenceValue )
    {
        if ( not mNormalized )
        {
            mReferenceValue = aReferenceValue;
            mNormalized     = true;
        }
    }
  ...
            else if ( tNormalization == "design" )
            {
                mIQIs( tRequestedIQIIndex )->set_reference_value( mGlobalIQIVal( tRequestedIQIIndex )( 0 ) );
                mGlobalIQIVal( tRequestedIQIIndex )( 0 ) = 1.0;
            }
  ...
        // Normalization
        if ( gLogger.mIteration == 0 )
        {
            this->normalize_IQIs();
        }
  ...
        for ( uint iIQIIndex = 0; iIQIIndex < mNumCriteria; iIQIIndex++ )
        {
            tVal( iIQIIndex )( 0 ) = sum_all( tVal( iIQIIndex )( 0 ) );
        }
  ```
- **Fix:** apply 1/`mReferenceValue` to the value and to dQI/du, dQI/dp at assembly, and compute the "design" reference from the
  globally summed iteration-0 value. Use `tIQIIndex`. Alternatively remove the option and raise an error when it is set. Update
  `UT_FEM_IQI_Strain_Energy.cpp:278`.
- **Regression test:** a 2-rank ctest with `normalization = "design"`:
  - assert `sum_all( value ) == 1.0` at iteration 0;
  - assert that iteration 1 equals Q1/Q0 and that dQ/dp is scaled by 1/Q0;
  - a numeric reference of 2.0 halves both the value and the gradient;
  - cover requested IQIs in a different order from their definition.

## MB-072 · FEM/INT: NODAL output of normal-dependent IQIs uses a stale or empty normal
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` C-2; the same root
  as A-8 and B-13). Origin: `b033c7ca3`, wunsch, 2023 ("Support for output of IQIs on (dbl) side-sets"), intended (documented): the
  asserts added in `517296282` say such IQIs "must be evaluated elementally" (`ORIGINS.md`). Enforced only in debug builds.
- **Exposure:** 0 decks meet the trigger. 128 files use normal-dependent IQIs (TRACTION, THERMAL_ENERGY_DIFFUSIVE_FLUX, GAP,
  CONTACT_PRESSURE_*), as criteria, GLOBAL or FACETED_*, with one exception. EXA `Parabolic_Indenter_NonLinear.cpp:940-942` requests
  GAP as NODAL on a NONCONFORMAL_SIDESET (`:793-794`). It is unaffected: GAP takes its normal from the IG geometry interpolator
  (`IQI/cl_FEM_IQI_Gap.cpp:37`), and the nonconformal NODAL path positions that interpolator. The NODAL ContactPressure block next to
  it (`:934-937`) is commented out.
- **Symptom:**
  - `IQI::mNormal` is a persistent member, set only by `Element_Sideset` and `Element_Double_Sideset` (the integration and VIS-integral
    paths).
  - `compute_nodal_QIs_standard` (bulk and side sets) and `compute_nodal_QIs_double_sided` position only the IP field interpolators.
    They set neither the normal nor the IG point.
  - `compute_nodal_QIs_nonconformal_side` positions the IG point (`set_space_time_from_local_IG_point`, `:1457`) but sets no normal.
  - So TRACTION, JUMP_TRACTION, THERMAL_ENERGY_DIFFUSIVE_FLUX and other IQIs that read `mNormal`, requested as NODAL fields, give
    σ(x_v)·n_last at every vertex, where n_last is the normal of the last facet evaluated, possibly on another set.
  - If nothing set a normal before, `traction( mNormal )` gets a 0×0 normal: a debug assert, and an out-of-bounds read in
    `flatten_normal_*` in release.
  - IQIs that read the IG geometry interpolator (NORMAL_VECTOR, and on the standard path CONTACT_PRESSURE and GAP) get the value at
    the last quadrature point.
  - VIS computes NODAL (enum 1) before GLOBAL (enum 7) on each set (`projects/FEM/VIS/src/cl_VIS_Output_Manager.cpp:848-895`).
- **Code** (`ELEM/cl_FEM_Interpolation_Element.cpp:1366-1367, 1392-1397`; `IQI/cl_FEM_IQI_Traction.cpp:40-47`;
  `IQI/cl_FEM_IQI_Thermal_Energy_Diffusive_Flux.cpp:45`, origin/develop):
  ```cpp
            // set vertex coordinates for field interpolator
            mSet->get_field_interpolator_manager()->set_space_time( tNodalPointLocalCoords );
  ...
                // reset the requested IQI
                tReqIQI->reset_eval_flags();

                // compute quantity of interest at evaluation point
                Matrix< DDRMat > tQINodal( 1, 1, 0.0 );
                tReqIQI->compute_QI( tQINodal );
  ...
        MORIS_ASSERT( mNormal.numel() > 0,
                "IQI_Traction::compute_QI() - "
                "Normal is not set. IQIs requiring a normal must be evaluated elementally "
                "and averaged such that there is a well-defined normal." );

        // evaluate traction
        const Matrix< DDRMat >
                tTraction = tCMLeader->traction( mNormal );
  ...
        aQI = -1.0 * dot( tCMDiffusion->flux(), mNormal );
  ```
- **Fix:** reject normal-dependent IQIs as NODAL with a `MORIS_ERROR` at VIS setup, and clear `mNormal` before the NODAL loop. Or
  compute a vertex normal (averaged over the vertex's facets), set it, and position the IG interpolator at the vertex. Make the
  existing asserts `MORIS_ERROR`.
- **Regression test:** a quarter-circle Neumann side set with σ = constant uniaxial. Request TRACTION as NODAL and as FACETED_AVG, and
  assert that the nodal t_x at every vertex equals σ n_x(vertex) to 1e-10. Today they all equal σ n_x of one facet. The alternative
  is an error at setup. B-13: `compute_QI` without `set_normal` raises an error.

## MB-073 · FEM/INT: EIGEN_VECTOR read at a stale point outside the criterion and NODAL paths
- **Status:** open (2026-10-22), confirmed by code reading. Found by the IQI audit (`studies/iqi_audit/FINDINGS.md` C-3). Origin:
  `dda7b2636`, maute, 2023 ("added support for eigen vectors in IQI computation"), a slip (`ORIGINS.md`). Exposure: 0 decks meet the
  trigger. 80 files use EIGEN_VECTOR, all as NODAL VIS fields, which position the point correctly: the ghost_study decks in
  `studies/xtk_holistic_review/spg_explainer/` and EXA `Column_Buckling.cpp:541` and `Column_Buckling_Immersed.cpp:683`.
- **Symptom:**
  - `IQI_Eigen_Vector` reads the eigen field interpolator manager.
  - That manager is positioned only in `Element_Bulk::compute_QI` (`:486-490`) and the NODAL routines (`Interpolation_Element.cpp:1372`,
    `:1467-1470`).
  - `compute_quantity_of_interest_global` (`:839`), both `compute_quantity_of_interest_elemental` overloads (`:933`, `:1013`) and
    `compute_dQIdu` (`:782`) position only the main manager.
  - So an EIGEN_VECTOR as GLOBAL, ELEMENTAL_INT/AVG or field output is the eigenvector at one leftover point (the last NODAL vertex or
    criterion quadrature point), constant over the quadrature loop. It is not the integral or the average.
- **Code** (`ELEM/cl_FEM_Element_Bulk.cpp:485-490` vs `:838-839`; `IQI/cl_FEM_IQI_Eigen_Vector.cpp:86-87`, origin/develop):
  ```cpp
            // if eigen vectors
            if ( mSet->mNumEigenVectors > 0 )
            {
                // set evaluation point for interpolators (FIs and GIs)
                mSet->get_field_interpolator_manager_eigen_vectors()->set_space_time_from_local_IG_point( tLocalIntegPoint );
            }
  ...
            // set evaluation point for interpolators (FIs and GIs)
            mSet->get_field_interpolator_manager()->set_space_time_from_local_IG_point( tLocalIntegPoint );
  ...
        Field_Interpolator* tFI =
                mLeaderEigenFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ), mEigenVectorIndex );
  ```
- **Fix:** position the eigen manager next to the main manager in every quadrature loop of `Element_Bulk`, and in the side-set elements
  if they gain eigenvector support. Better, do it in one helper that sets both.
- **Regression test:** a 2×1 QUAD4 patch with eigenvector coefficients set to φ = x. Assert that ELEMENTAL_AVG equals each element's
  centroid x (exact for linear φ) and that GLOBAL equals ∫x dA.

---

## MB-074 · FEM/INT: analytic sensitivity mode ignores the deck's FD scheme and step
- **Status:** merged (#61, develop `5eb5acf63`, 2026-10-22). Previously: fix-written, verified (PR #61, branch `fix/mb062-analytic-property-variation` @ `564fc0bb9`, worktree `moris-geo-sens-props`,
  2026-10-22). Found while fixing MB-062: the SWEEP fallback arm (`se_vary_an`) gave 1.99e-7 against FD 6.73e-9 until this was fixed;
  after, the two arms agree bitwise. Origin: `00023a85fb` (lisenoel, 2021), upstream.
- **Symptom:** with `is_analytical_sensitivity = true`, every part of dR/dp, dQI/dp that is still computed by finite differences
  (material pdvs, the FD geometry fallback of IWGs/IQIs without an analytic integrand, the cluster measures of Nitsche SPs) uses the
  default scheme POINT_1_FORWARD and default step instead of `finite_difference_scheme` / `finite_difference_perturbation_size`.
  Exposure: decks in analytic mode that set a non-default scheme, e.g. `Sphere3D_Precision.cpp` (`SPH_SA_FD_*`), and any analytic-mode
  run whose gradients mix analytic and FD parts. EXA sensitivity ctests not run against the fix.
- **Code** (`projects/FEM/INT/src/ELEM/cl_FEM_Element.hpp:211-228` on develop):
  ```cpp
  if ( tIsAnalyticalSA )
  {
      m_compute_dRdp  = &Element::select_dRdp;
      m_compute_dQIdp = &Element::select_dQIdp;
  }
  else
  {
      ...
      mSAFDScheme       = mSet->get_finite_difference_scheme_for_sensitivity_analysis();   // only read on the FD path
      mSAFDPerturbation = mSet->get_finite_difference_perturbation_size();
  }
  ```
- **Fix:** read the scheme and perturbation in both branches.

## MB-076 · FEM/VIS: nonconformal output IQIs evaluated at the frozen reference pairing
- **Status:** fix-written, PR #64 (branch `feature/nonconformal-output-gap-remap-develop` @ `f3d00b51b`, worktree
  `moris-pr-nc-gap-remap`, 2026-10-23). Ported from the never-PR'd `feature/nonconformal-output-gap-remap` (WIP snapshot
  `177880a8c`, 2026-07-18). `UT_FEM_Element_Nonconformal_Sideset_Gap` passes; INT/GEN/OPT/XTK suites pass. Not yet rerun on
  an EXA contact deck with VIS output.
- **Symptom:** in a nonconformal contact run, the Exodus `GAP` / `CONTACT_PRESSURE` / traction fields are ~0 everywhere,
  while the IWG residual sees the correct gap. Fails silently.
- **Cause:** the VIS clusters copy the leader/follower point pairs from the FEM clusters when the VIS mesh is built
  (near the reference configuration). The solver's ray remap runs only inside the IWG residual/Jacobian evaluation, so the
  output path positions the follower field interpolator at the stale pairing.
- **Code** (`projects/FEM/INT/src/ELEM/cl_FEM_Interpolation_Element.cpp:1460-1462` on develop `c45d336ee`):
  ```cpp
  // the coordinate of the follower side is known from the mapper
  // it is not a IG point, but the method can still be used (a parametric point on the follower side geometry)
  mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
          ->set_space_time_from_local_IG_point( tFollowerCoordinates.get_column( iNode ) );
  ```
- **Fix:** output-only virtual `Element_Double_Sideset::remap_follower_for_output()` (no-op), overridden by
  `Element_Nonconformal_Sideset` to re-fire the ray with `remap_nonconformal_ray()`; called after positioning the
  interpolators in the double-sideset IQI loops and the VIS nodal loop. 2D line facets only.
- **Workaround (until merged):** post-recover the gap from the displacement field (see the
  nonconformal-contact-iqis-zero-at-vis memory).

## MB-077 · OPT: restart on a different rank count loads the ADVs in file order (reorder overwritten)
- **Status:** merged (#70): landed on develop as `2b35a72b1` (2026-10-24). Was: fix-written (2026-10-24): PR #70 (`fix/opt-restart-layout`, be4c05d84), together with MB-081. A serial-written restart file now gives identical criteria on 1, 2 and 4 ranks; develop on 2 ranks starts from a scrambled design (SE 21783.9 vs 2998.18). Evidence: `runs/studies/xtk_pathology/2026-10-24_mb077_mb081/`.
- **Symptom:** a run restarted from an `ADV_Alg_0_Iter_N.hdf5` written on a different number of ranks starts from a
  scrambled design: the level-set field differs from the written design everywhere (max |dLVLSET| = 1.0 on a [0, 1]
  field), with no error. Restarts on the same rank count, and the uniform iteration-0 design, look fine, which hides it.
  Evidence: two partitions (np2, np4) restarted from the same serial-written file produce byte-identical (wrong)
  criteria, because both take the file order; serial restarted from it reproduces the design exactly (dLVLSET = 0).
- **Code** (`projects/OPT/src/cl_OPT_Problem.cpp:521-558` on develop `4a7f493a4`; present before the 2024 clang-tidy
  refactor `ffa74e15d`):
  ```cpp
  // reorder adv for new proc layout. FIXME see above
  if ( par_size() != tParSizeFile )
  {
      ...
              mADVs( tIndxed )        = tRestartADVs( Ik );      // reordered by IjklId
      ...
              mADVs( Ik )        = tRestartADVs( Ik );           // id not found: written at the file position
      ...
  }

  // Copy restart vectors on member variables
  mADVs        = tRestartADVs;                                    // overwrites the reorder with the file order
  mUpperBounds = tRestartUpperBounds;
  mLowerBounds = tRestartLowerBounds;
  ```
- **Fix direction:** copy the file vectors only when `par_size() == tParSizeFile` (and per the existing FIXME, reorder by
  IjklIds whenever the layout can differ, i.e. always); make a missing id an error instead of writing at the file
  position. Regression test: write ADVs on np1, restart on np2/np4, compare the loaded ADVs by IjklId.
- **Reproducer:** `runs/studies/heatmethod_serial_parallel/same_design/` (`same_design.sh`: the EXA SIMP hole-seeding
  deck restarted from serial ADVs on np 1/2/4/8; `compare_same_design.py`). Evidence: `runs/moris_bugs/MB-077/`.
- **Workaround (until fixed):** restart only on the rank count that wrote the file, or rewrite the restart file in the
  target run's own layout (its IjklIds order and NumProcs) before restarting.
- **Impact:** any warm restart that changes the rank count (stage chaining probe → smoke → full, recovery restarts,
  `RUNE_RESTART_FILE`) silently starts from a scrambled design.

## MB-078 · FEM/INT: a repeated set name in `mesh_set_names` assembles the IWG twice, silently
- **Status:** merged (#66), warning only: landed on develop as `8bcb7dd21` (2026-10-24, squash of
  `fix/fem-duplicate-iwg-warning` `3482c405a` + `b07e86572`; warn once per run per IWG/IQI and set). Repeated
  definitions are still assembled; the EXA deck strings are unchanged. `Set_User_Info::add_IWG` / `add_IQI` log
  `IWG <name> is defined more than once on mesh set <set> and will be assembled once per definition; check its mesh_set_names.`
  Assembly is unchanged. Checked:
  - the hole-seeding deck prints both IWG warnings with SE 3044.21, unchanged;
  - a deck without repeats prints none;
  - INT-test.exe 60976/305 passes, the same as the baseline.

  Evidence: `runs/studies/heatmethod_serial_parallel/dup_warning_e2e/`.
  Severity low: the input is a deck mistake, but MORIS accepted it
  silently and the physics changes.
- **Symptom:** listing a side set twice in an IWG's `mesh_set_names` multiplies that IWG's contribution. The EXA deck
  `Level_Set_Beam_SIMP_Hole_Seeding` sets its traction on `"SideSet_4_n_p1,SideSet_4_n_p1"`. On a fixed design that
  gives SE 3044.21; listed once it gives 761.053 (ratio 4.0000: load doubled). The Dirichlet set is duplicated the same
  way (`"SideSet_2_n_p1,SideSet_2_n_p1"`), which doubles the Nitsche terms.
- **Code** (`projects/FEM/INT/src/CORE/cl_FEM_Model_Initializer_Legacy.cpp:909-948`, `cl_FEM_Set_User_Info.hpp:423-426`
  on develop `4a7f493a4`):
  ```cpp
  auto tMeshSetNames = string_to_vector< std::string >( tIWGParameterList.get< std::string >( "mesh_set_names" ) );
  for ( auto &tMeshSetName : tMeshSetNames )
  {
      ...
      else
      {
          // if the fem set already exister, only add the IWG
          this->mSetInfo( aMeshToFemSetIndex[ tMeshTuple ] ).add_IWG( tIWG );   // same IWG again
      }
  }
  // Set_User_Info::add_IWG: mIWGs.push_back( aIWG );
  ```
  The IQI path (`:982` onward) has the same pattern.
- **Fix direction:** deduplicate `tMeshSetNames` (or error with the IWG name and the repeated set).
  Regression test: one IWG with a repeated set name assembles the same residual as with the set listed once.
- **Reproducer:** `runs/studies/heatmethod_serial_parallel/d193_single_set/` vs `same_design_native/d234np2/np1/`
  (first criteria line).
- **Related:** the EXA deck string looks like a typo for `"SideSet_4_n_p1,SideSet_4_c_p1"`, which would also have
  prevented D-193 (traction lost when an interior pinhole cuts the loaded element).
- **Effect on this deck** (`runs/studies/heatmethod_serial_parallel/dirichlet_dup/`, serial, fixed designs, SE at the
  first criteria line):

  | design | traction ×2, Dirichlet ×2 (original) | traction ×1, Dirichlet ×2 | traction ×1, Dirichlet ×1 |
  |---|---|---|---|
  | uniform solid (np1 ADV iteration 0, fixed-design deck) | 1524.97 | 381.242 | 381.068 |
  | np2 iteration 234 | 3044.21 | 761.053 | 760.736 |

  - Duplicated traction: SE × 4.0000 on both designs (linear elasticity, load × 2), a pure rescaling. The deck
    normalizes SE by a constant, `tInitialStrainEnergy = pick( 1437.01, 7043.0 )`, measured with the duplicate. So the
    objective's penalty shares are as intended, and the design path is that of the deduplicated deck with the
    normalizer / 4.
  - Duplicated Dirichlet (Nitsche terms × 2): SE +0.042-0.046 % (more compliant), SE gradient norm +0.02 %
    (66.048 vs 66.034 at iteration 234).
  - It matters only when the strings are fixed: then `tInitialStrainEnergy` must be divided by 4, or the perimeter and
    regularization shares grow 4×. The EXA test's reference SE values change by the same factor.
- **Origin of the deck string** (git history, 2026-10-24):
  - Introduced with the deck in `13a7e97bb` (schmidt, 2021-08-18, "added hole seeding example ( no checks yet )").
  - Template: `One_Element_Mat_Geo_Pdv_sweep.cpp`, the only earlier deck with the misspelled `tDirchletSets`. It reads
    `"SideSet_4_n_p0,SideSet_4_c_p0"` / `"SideSet_2_n_p1,SideSet_2_c_p1"` and wires Dirichlet → `tDirchletSets`,
    traction → `tLoadSets`.
  - The copy turned `_c_` into `_n_` in both strings and swapped the wiring (Dirichlet IWG on `tLoadSets`, traction on
    `tDirchletSets`). Both are in the initial commit.
  - Later touches only renamed or realigned the strings: `73724dee9` (wunsch, 2023-05-30, typo `tDirchletSets` →
    `tDirichletSets`) and `dcc6faa1b` (momo, 2023-10-23).
  - **Already fixed upstream on a side branch:** `f2d2d661c` (kkmaute, 2024-11-22, "some changes") sets both to
    `_n_p1,_c_p1`. In the same commit the deck restarts from `ADV_Alg_0_Iter_15.hdf5`, and the EXA test gains
    `test_pause()` ("for debugging in parallel"). It is only on `upstream/directSA` (201 commits ahead of `upstream/main`,
    tip 2025-03-09). It is not in `upstream/main` and not in `origin/develop`.

## MB-079 · FEM/INT: stress IQIs skip the index check when the first dof group has one member
- **Status:** open (2026-10-24), confirmed by code reading. Found by the medium code review of PR #63 (the PR keeps the behaviour; the
  check now lives in `IQI_Stress::eval_stress_value`). Exposure: 0 EXA decks. All 7 EXA NORMAL / SHEAR / PRINCIPAL stress IQIs set
  `vectorial_field_index` and list the displacement group first. Workspace decks not swept.
- **Symptom:** the guard on `mIQITypeIndex` is conditioned on `mLeaderDofTypes( 0 ).size() > 1`, the size of the first leader dof
  group, not on the stress type. With a single-member first group (a thermo-elastic IQI declared `"TEMP;UX,UY"`) and no
  `vectorial_field_index` (default −1), the measure is called with index −1 + 1 = 0:
  - NORMAL_STRESS reads standardized entry `aStressIndex − 1` = (uint)−1: out of bounds (unchecked in release, `ARMA_NO_DEBUG`);
  - SHEAR_STRESS reads entry 2 (σzz), not a shear stress, silently;
  - PRINCIPAL_STRESS falls to the `default` of the index switch and aborts with "Only 1st and 2nd principal stresses known for 2D"
    (or "1st, 2nd, and 3rd" in 3D), which does not name the missing index.
  Same guard in `IQI_Max_Stress` and `IQI_Zienkiewicz_Zhu` on develop.
- **Code** (`IQI/cl_FEM_IQI_Stress.cpp:41-49`, `:248-255`, origin/develop; same shape at `:86-94`):
  ```cpp
        // check if dof index was set
        if ( ( mLeaderDofTypes( 0 ).size() > 1 ) && ( mStressType != Stress_Type::VON_MISES_STRESS ) && ( mStressType != Stress_Type::STRESS_VECTOR ) )
        {
            MORIS_ERROR( mIQITypeIndex != -1,
                    "IQI_Max_Stress::compute_QI - mIQITypeIndex not set, but is needed for principal, normal and shear stresses." );
  ...
    void IQI_Stress::eval_normal_stress( uint aStressIndex, Matrix< DDRMat >& aStressVector )
  ...
        aStressVector.set_size( 1, 1, tStressVector( aStressIndex - 1 ) );
  ```
- **Fix:** condition the check on the stress type only, using `IQI::get_component_index( 3, -1 )` (an unset index is an error) for
  NORMAL / SHEAR / PRINCIPAL. After PR #63 there is one place to change: `IQI_Stress::eval_stress_value`, which `IQI_Max_Stress` and
  `IQI_Zienkiewicz_Zhu` inherit.
- **Regression test:** an `IQI_Stress` NORMAL_STRESS with dof groups `{ { TEMP }, { UX, UY } }` and no index: expect a `MORIS_ERROR`
  (today: out-of-bounds read). Control: with index 1, σyy.

## MB-080 · GEN: debug validation calls every user-defined field with a 3D point
- **Status:** merged (#69): landed on develop as `c455c797b` (2026-10-24). Was: fix-written (2026-10-24): PR #69 (`fix/gen-udf-probe-dimension`, ee2abbaab). The functions are checked on the coordinates they are called with. Also: in release the sensitivity probe was a silent out-of-bounds read (Vector::operator() checks bounds only in debug). Evidence: `runs/studies/xtk_pathology/2026-10-24_mb080/`.
- **Symptom:** on develop `d355343dd`, debug `XTK-test.exe "MB-031*"` fails both 2D cases ("MB-031 2D void inside one
  background cell, all corners positive" and its sibling) with
  `vector::_M_range_check: __n (which is 2) >= this->size() (which is 2)`. Release passes. Backtrace:
  `User_Defined_Field::User_Defined_Field` → `validate_user_defined_functions` (`cl_GEN_User_Defined_Field.cpp:102`) →
  test field `mb031_disc` → `distance_minus_radius` (`UT_XTK_Interior_Cut_Detection.cpp:53`) indexes the 2-entry centre
  with `iDim = 2`, because the probe point has 3 components.
- **Code** (`projects/GEN/src/ADV/fields/cl_GEN_User_Defined_Field.cpp:100-102`, develop `d355343dd`):
  ```cpp
  MORIS_ERROR( get_field_value_user_defined, "No field evaluation function was provided to a user-defined field." );
  MORIS_ASSERT( std::isfinite( this->get_field_value_user_defined( { { 0.0, 0.0, 0.0 } }, mFieldVariables ) ),
          "There is an error in a user-defined field (field evaluates to nan/infinity)." );
  ```
- **Fix direction:** probe with a point of the field's spatial dimension (or skip the probe when the dimension is
  unknown). The test field could also loop over `aCentre.size()`, but the probe is what is wrong: deck user functions
  written for 2D can fail the same way in debug.
- **Reproducer:** `moris-cache-warm-develop/build_dbg/projects/XTK/test/bin/XTK-test.exe "MB-031*"` (cwd
  `build_dbg/projects/XTK/test`). Log: `runs/studies/xtk_pathology/2026-10-24_ipv_ids/review_fix_dbg/mb031_cache-warm-develop.log`.

## MB-081 · OPT: the restart file is opened with write permission; parallel restarts abort in debug
- **Status:** merged (#70): landed on develop as `2b35a72b1` (2026-10-24). Was: fix-written (2026-10-24): PR #70 (with MB-077): the restart file is opened read-only. Debug 2-rank restart runs; develop's debug build stops on the assert.
- **Symptom:** a 2-rank debug run with `restart_file` set stops right after
  `Problem - Reading ADVs from restart file: ADV_restart.hdf5` with
  `Assertion par_size() > 1 && !aReadOnly ? aAddParExt : true failed` in `open_hdf5_file()`
  (`projects/MRS/IOS/src/HDF5_Tools.hpp:107`).
- **Code** (`projects/OPT/src/cl_OPT_Problem.cpp:481-483`, develop `d355343dd`):
  ```cpp
  // Note: only processor 0 reads this file; therefore no parallel file name extension is used
  hid_t tFileID = open_hdf5_file( mRestartFile, false );    // aReadOnly defaults to false
  ```
- **Fix direction:** `open_hdf5_file( mRestartFile, false, true )`: read-only, which is what the comment describes.
- **Reproducer:** `runs/studies/xtk_pathology/2026-10-24_ipv_ids/review_fix_dbg/deck_np2/` (hole-seeding deck from
  `ADV_restart.hdf5`, `mpirun -n 2` with a debug `moris`). Related: MB-077 (same function, restart across rank counts).

## MB-082 · XTK: #67 numbers a merged enriched IP vertex once per attached cell (debug assert)
- **Status:** merged (#68): landed on develop as `78e70a57d` (2026-10-24, squash of deffa8872 + 56cb704ea).
  Introduced by #67 (`68efd8f97`, squash-merged 2026-10-24 before the review fix was pushed). Also on
  `work/develop-local` since a7059b01e and in production `build_opt_20261024b`, which is release-only and so harmless there.
- **Symptom:** in debug, any run that builds an enriched IP mesh stops on
  `Assertion tVertexTracker( tVertexIndex ) failed` (`cl_XTK_Enriched_Interpolation_Mesh.cpp:2510`). develop's
  `[XTK_PATHOLOGY_MPI_IPV_IDS]` fails 12/12 on 2 ranks in debug. In release, a vertex on k cells is listed and numbered
  k times (the last id wins, k - 1 ids are wasted: e.g. 225 owned vertices, ids up to 233), and a not-owned vertex is
  requested k times with the last answer kept unchecked. Ids still agree across ranks and solutions are unchanged.
- **Code** (develop `68efd8f97`, `sort_unzipped_vertices_into_owned_and_not_owned`):
  ```cpp
  for ( uint iCell = 0; iCell < tNumUIPCs; iCell++ )          // a merged UIPV is on several UIPCs
      for ( uint iVert = 0; iVert < tNumUIPVsOnCell; iVert++ )
      {
          MORIS_ASSERT( tVertexTracker( tVertexIndex ), "... already been assigned one." );
          ...
          mOwnedUnzippedVertices.push_back( tVertexIndex );   // once per cell
  ```
- **Why the #67 tests missed it:** they ran only in release and checked cross-rank id agreement, which still held.
  #68 adds a release-visible numbering check (owned ids exactly 1..N, no id repeated on a rank). That check fails
  #67's code in release in 12/12 pathology cases and 40/40 fuzz cases.
- **Evidence:** `runs/studies/xtk_pathology/2026-10-24_ipv_ids/{review_fix_dbg,review_fix_opt,numbering_check,followup_pr}/`.

## MB-075 · XTK/HMR: SIGSEGV when a cut coarse cell borders an uncut finer cell
- **Status:** merged (#71): landed on develop as `30c2cc00a` (2026-10-24). Was: PR #71 (`fix/mb075-transition-neighbor-facet`, 27f66d1cc, from develop
  `78e70a57d`). Design: `studies/xtk_holistic_review/findings/18_mb075_transition_segfault_fix.md`.
- **Symptom:** decomposition crashes when the interface leaves the refined band into the coarse cells next to it.
  Release: SIGSEGV (rc 139). Debug: `Null ptr on facet that should have children facets`. 1, 2 and 4 ranks, stored and
  analytic fields. No working run passed through this configuration before the fix.
- **Cause:** HMR fills the transition location only when looking from the coarse side
  (`cl_HMR_Mesh.cpp:1161`, `:1459-1470`). From the fine side the facet index is the fine cell's own facet, under which
  none of the coarse cell's child facets are stored.
- **Code** (develop, `cl_XTK_Integration_Mesh_Generator.cpp`, `construct_subphase_neighborhood`, same-level branch):
  ```cpp
  this->collect_subphases_attached_to_facet_on_cell(
          aCutIntegrationMesh, tOtherCell, tNeighborOrdinal,
          tFacetIndex,                               // my (fine) facet ...
          aFacetConnectivity,
          ( *aBgFacetToChildFacet )( tFacetIndex ),  // ... whose child-facet list is null
          ... );
  // collect_subphases_attached_to_facet_on_cell: MORIS_ASSERT( aBgFacetToChildrenFacets != nullptr, ... );  // debug only
  //                                              aBgFacetToChildrenFacets->size();                          // SIGSEGV
  ```
- **Fix (#71):** use the neighbour's own facet across levels; the null list is a `MORIS_ERROR`; a guard errors on an
  ambiguous cross-level link (two same-phase pieces of the cut coarser cell on the face; finding 21). With this
  subdivision scheme geometry does not seem to reach it, so it is a safety net.
- **Reproducer:** hidden `[XTK_CONNECTIVITY_TWO_LEVEL]` with `MORIS_XTK_CONNECTIVITY_CRASH_CASES=1` on develop (case A2);
  regression test `[MB-075]`. Evidence: `runs/studies/xtk_pathology/2026-10-24_mb075/`, earlier draft
  `runs/studies/xtk_pathology/2026-10-21_connectivity/defects/transition_coarse_cut_segfault/`.
- **Workaround (until merged):** keep the interface inside the refined band (widen the refinement buffer).
