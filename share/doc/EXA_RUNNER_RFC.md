# RFC: EXA-test.exe — one runner, test case = data

- **Status:** ACCEPTED 2026-07-02 — owner decisions recorded below; implementation started (Step 0 + Laplace_2D pilot)
- **Date:** 2026-07-01

## Owner decisions (2026-07-02)

1. **Catch2 test-name renames (§2.3c): approved.**
2. **The 8 disconnected examples: reconnect** during their physics-directory batch. The §5.3
   ctest parity check therefore becomes: old set ⊆ new set, with the only additions being
   exactly the reconnected examples' registrations (list them in the batch PR).
   **Amended 2026-07-02 (thermal batch finding):** the two thermal disconnects
   (`Channel_with_Four_Cylinders_Transient`, `Two_Channels_with_Separation_Wall_Transient`)
   do not compile — their DECKS use the pre-modern `ParameterList` API. That is evidently why
   they were disconnected. Their leaves are migrated and ready, but reconnection is deferred
   until the decks are ported to `Module_Parameter_Lists` (separate deck-modernization work,
   no build-time/size payoff → needs its own sign-off under the quantification policy).
   Protocol for remaining batches: test-compile a disconnected example's deck FIRST; reconnect
   only if it compiles, otherwise migrate the leaf, leave it disconnected, and report.
3. **The 2 stale `src/` copies: delete in Step 7.**
4. **Examples build whenever `MORIS_USE_EXAMPLES` is ON (its default).** No preset or
   configuration may flip examples off by default — the size/time win must come from this
   RFC's single-runner design, not from turning examples off.
- **Author:** cmake-expert (advisory; read-only — all patches below are drafts for an implementer)
- **Scope:** `moris/projects/EXA/**` build restructuring only. The `dynamic_link_input()` input-`.so` machinery, `MORIS_USE_EXAMPLES` gating, and all ctest registration names are explicitly preserved.
- **Related:** `moris/CMAKE_AUDIT.md` §"EXA link cost" (2025-09-28), `moris/share/doc/CMAKE_REVIEW.md` (2026-07-01), baseline `$MORIS_RUNS_DIR/benchmarks/build_metrics/baseline_build_opt_2026-07-01.json`

---

## Verdict

The 57 statically linked EXA example executables are 57 copies of the same ~34 MB LTO link whose only differences are one ~1.5 s test-case TU each (median TU compile 1.5 s, `.ninja_log`). Merging them into one `EXA-test.exe` is mechanically feasible today: all 65 live leaf `CMakeLists.txt` carry a byte-identical dependency list (md5-verified), the test-case TUs share one `main()` already, and Catch2 v2.13.9's tag filtering gives per-example selection for free. The make-or-break constraint is not CMake — it is the **dlopen symbol side channel**: every example's input `.so` resolves file-scope globals (`gInterpolationOrder`, `gTestCaseIndex`, …) against the executable's `-rdynamic` symbol table, so those globals must stay at global scope with their exact names, defined exactly once. The proposed fix — one canonical `EXA_Globals.{hpp,cpp}` for the ~23 deck-visible globals, per-TU namespace wrapping for everything else, and three TEST_CASE-name renames — handles every collision found in the survey of all 66 TUs. Predicted: **−1.75 GB executables and ≈ −134 CPU-min of link per clean `build_opt` pass, ≈ −20 GB in `build_dbg`**, at the cost of **+5 to +55 s on the single-example edit-relink loop** and a migration that must touch every test-case TU once (mechanically, per a fixed recipe).

---

## 1. Design summary, predicted deltas, and new costs

### 1.1 Design in one paragraph

Replace the per-example `add_executable(Main_<X>)` with a single `EXA-test.exe` that links `runner_main.cpp` + `EXA_Globals.cpp` + every migrated `example_test_case.cpp`. Each example's ctest registration keeps its exact current name and PROCS list but its COMMAND becomes `mpirun -n <P> $<TARGET_FILE:EXA-test> "[EXA_<Name>]"` with the same per-example `WORKING_DIRECTORY` (`<leaf-binary-dir>/bin`), where its input `.so`(s) — built by the untouched `dynamic_link_input()` targets — continue to land. A new `moris_add_example()` function shrinks each leaf `CMakeLists.txt` to one call; leaves accumulate their TU into a global property that a final `runner/` directory turns into the executable, so old-style and new-style leaves coexist during migration.

### 1.2 Measured baseline (all numbers re-measured 2026-07-01 in this workspace)

| Quantity | Value | Evidence |
|---|---|---|
| EXA executables built (opt) | 57 (`67 add_executable` in tree: 65 live + 2 stale `src/` copies; 8 live leaves never `add_subdirectory`'d) | `find build_opt/projects/EXA -name "*.exe"`; grep of leaf CMakeLists |
| EXA exe bytes (opt) | **1,793 MB** of the baseline's 2,414.9 MB total executables (83 files). (The motivating figure "2,021 MB" was a `du`-based measure; both methods are re-run in §1.5's protocol.) | `find … -printf "%s"`; `baseline_build_opt_2026-07-01.json /sizes/executables` |
| EXA exe bytes (dbg) | **20,643 MB** across 37 exes (build_dbg EXA dir: 22 GB) | `find build_dbg/projects/EXA -name "*.exe" -printf "%s"` |
| Per-exe link time (opt, `-flto=4`) | 145–151 s (top: `Main_HeatConduction.exe` 150.8 s) | `build_opt/.ninja_log` |
| EXA link CPU per clean opt build | 57 × ≈145 s ≈ **8,265 s ≈ 138 CPU-min** (baseline: link = 9,329.6 s over 177 edges = 74.4 % of 209.0 total edge CPU-min) | `.ninja_log`; baseline JSON |
| Test-case TU compile time | median **1.5 s**, max 2.5 s | `.ninja_log` |
| ctest registrations (opt) | **73** (not 67 — the in-tree 67 `add_test` *calls* expand by per-example PROCS lists; 73 is what `build_opt` actually registers) | `grep -c add_test build_opt/projects/EXA/**/CTestTestfile.cmake` → 146/2 |
| Input `.so` files (opt) | 68, from 64 `dynamic_link_input()` calls (some leaves build several: `Field_Example` 3, `…_Restart` 2) | `find build_opt/projects/EXA -name "*.so"` |

### 1.3-VERIFIED — Final A/B measurement (2026-07-02, migration complete)

Two full clean builds, identical config (BUILD_ALL=ON, EXAMPLES=ON, Ninja -j4, CI-mirror flags),
pure committed states via git worktrees (`7ea5ef129^` vs `7017b0a77`), same machine:

| Metric | Before | After | Δ (measured) | Predicted |
|---|---|---|---|---|
| Clean-build wall time | **51 m 35 s** | **22 m 03 s** | **−57%** | 25–32 min ✓ (beat) |
| Total edge CPU | 201.9 min | 83.5 min | **−118.4 CPU-min (−59%)** | ~−134 ✓ |
| Link share of build | 71.6% | 33.1% | −38.5 pts | ✓ |
| Executables | 2,411.9 MB (83) | 429.0 MB (27) | **−1,983 MB (−82%)** | −1,755 ✓ (beat) |
| Build tree | 3.89 GB | 1.62 GB | **−2.27 GB (−58%)** | ✓ |
| Static libs / shared libs | identical | identical | 0 — change touched only examples | ✓ |

Snapshots: `$MORIS_RUNS_DIR/benchmarks/build_metrics/exa-final-{before,after}_*.json`.
`build_dbg`'s ~−20 GB accrues on its next clean rebuild (examples were OFF there anyway;
now they can be ON for free).

### 1.3 Predicted deltas (superseded by 1.3-VERIFIED above; kept for the record)

| Metric | Before | After | Δ |
|---|---|---|---|
| `build_opt` executable bytes | 2,414.9 MB (83 files) | ≈ 660 MB (27 files: −57 EXA, +1 runner ≈ 40 MB) | **≈ −1.75 GB (−73 %)** |
| `build_opt` clean-build link CPU | 9,329.6 s | ≈ 1,250 s (−57 EXA links ≈ −8,265 s, +1 runner link ≈ +180 s) | **≈ −134 CPU-min (−87 % of link CPU)** |
| `build_opt` total edge CPU | 209.0 min | ≈ 75 min | **≈ −64 %** |
| Clean-build wall time (`ninja -j4`) | ≈ 58 min | *estimate* 25–32 min (link tail dominated; `-flto=4` links were never 4-way job-parallel anyway) | measure, do not quote until verified |
| `build_dbg` EXA exe bytes | 20,643 MB (37 exes) | ≈ 620 MB (1 dbg runner) | **≈ −20.0 GB** (build_dbg tree 36 GB → ~16 GB) |
| Core-lib change → EXA relink cost | 57 links ≈ 138 CPU-min | 1 link ≈ 3 CPU-min | **≈ −135 CPU-min per incremental core change** |
| ctest suite wall time | — | unchanged (same processes, same workloads, same working dirs) | ≈ 0 |

Runner size/link estimates: a single example exe is ≈ 34 MB / ≈ 148 s because it already contains the full static closure; the runner adds only ~66 more small TUs (median 1.5 s compile, ~100–300 KB of object each), so ≈ 40 MB / 150–200 s is the honest band.

### 1.4 New costs (explicit, so nobody discovers them later)

1. **Single-example iteration loop gets slower.** Today: touch `Laplace_2D/example_test_case.cpp` → 1.5 s compile + ≈ 148 s link of `Main_Laplace.exe` ≈ 150 s. After: 1.5 s + one runner link at 150–200 s ≈ **+5…+55 s**. This is the only scenario that regresses.
2. **Any example's test-case edit relinks the shared runner**, momentarily invalidating *all* EXA ctest entries (they now point at one file). Same-cost, wider blast radius.
3. **One serial link on the critical path.** 57 independent link edges become 1; ninja can no longer overlap EXA links with each other (it can still overlap the one link with compiles; and each old link already monopolized 4 LTO threads, so the practical parallelism loss is small — verify with the wall-time measurement).
4. **The 68 input `.so` targets remain** — no size or time change there (by design; that is the pymoris seam).
5. **Union linking:** the runner links the union of all example dependencies. Today that union equals every leaf's list (all 65 live leaves' `EXAMPLE_DEPENDENCIES` blocks are md5-identical). A future example needing an extra library imposes it on the runner as a whole.
6. **One-time migration labor:** every TU is edited once by a fixed mechanical recipe (§2.4); three Catch2-level test names are renamed (§2.3); ctest names are untouched.

### 1.5 Measurement protocol (owner policy compliance)

All claims above are re-verified before/after with `share/scripts/build_metrics.py` (`build_metrics.py <build_dir> --label <name>`; `--diff A.json B.json`). **Never run in `build_opt`/`build_dbg`** (pymoris fingerprints those caches); use throwaway build dirs. Handoff to **moris-build**:

```bash
spack env activate /home/doble/codes && source ~/.bashrc_moris
SCRATCH=$(mktemp -d)

# BEFORE (current main), and AFTER (migration branch), identically:
cmake -S /home/doble/moris_workspace/moris -B $SCRATCH/b -G Ninja \
      -DBUILD_ALL=ON -DMORIS_USE_EXAMPLES=ON            # mirror .github/workflows/ctest.yml:36
( cd $SCRATCH/b && /usr/bin/time -v ninja -j4 )          # scenario A: clean build (wall + CPU)
python3 share/scripts/build_metrics.py $SCRATCH/b --label exa-runner-{before|after} \
        --out $MORIS_RUNS_DIR/benchmarks/build_metrics/
# scenario B: single-example loop
touch projects/EXA/thermal/diffusion/Laplace_2D/example_test_case.cpp
( cd $SCRATCH/b && time ninja -j4 )
# scenario C: core-lib change → relink fan-out
touch projects/WRK/src/cl_WRK_Workflow_Factory.cpp
( cd $SCRATCH/b && time ninja -j4 )
# ctest parity (see §5.3)
python3 share/scripts/build_metrics.py --diff before.json after.json
rm -rf $SCRATCH
```

Report: Δ executables MB, Δ link seconds, Δ total edge CPU-min, and wall times for A/B/C. The `build_dbg` variant repeats A with the dbg flags from `ctest.yml:20`.

---

## 2. Mergeability findings (survey of all 66 test-case TUs)

Method: every `example_test_case.cpp` under `projects/EXA` was scanned for file-scope definitions, `extern "C"` functions, and `TEST_CASE` names; every input-deck `.cpp` was scanned for `extern` declarations.

### 2.1 The dlopen side channel is the binding constraint (do not namespace these)

Each test case builds a fake argv pointing at its input `.so` (`Laplace_2D/example_test_case.cpp:190`: `char tString2[] = "./Laplace.so";`) and calls `fn_WRK_Workflow_Main_Interface`, which dlopens it. The deck then resolves globals **back out of the executable**:

- `Laplace_2D/Laplace.cpp:25-27`: `extern uint gInterpolationOrder;  extern uint gTestIndex;`
- This works because example executables link with `-rdynamic` (verified in `build_opt/build.ninja`, `Main_Laplace` FLAGS: `… -std=c++17 … -rdynamic -fPIC …`), exporting the exe's global symbols to the dynamic loader.

Union of deck-side `extern` declarations across all EXA deck sources (regenerate with `grep -rhP '^\s*extern\s+(?!"C")' projects/EXA --include='*.cpp'` excluding test TUs):

- `uint`: `gInterpolationOrder` (20 decks), `gTestCaseIndex` (10), `gCaseIndex` (4), `gDim` (5), `gTestIndex`, `gOrder`, `tGeoModel`, `tDim`, `gLevelSetInterpolationOrder`, `gLagrMeshInterpolationOrder`, `gFEMInterpolationOrder`
- `bool`: `gPrintReferenceValues` (3), `gUseBspline`, `gUseMixedTimeElements`, `gUseBelosWithILUT`, `gInletVelocityBCFlag`, `gInletPressureBCFlag`, `gHaveStaggeredFA`, `gHaveStaggeredSA`
- `std::string`: `tOrder`, `tStressType`, `tOutputFileName`, `gPrecSolver`
- plus `moris::Logger gLogger` (defined in the runner main, as today).

**Consequence:** these names must exist at global scope, with these exact types, defined **exactly once** in the runner. They cannot be namespace-wrapped, `-D`-renamed, or made `static` — any of those changes the dynamic symbol and silently breaks every deck at dlopen. Types were cross-checked (e.g. `Column_Buckling/example_test_case.cpp:25-34` defines `std::string tOrder; uint tDim; std::string tStressType; std::string tOutputFileName;` matching its deck's externs); no same-name/different-type conflict exists in the tree today.

### 2.2 Collision inventory in the test-case TUs

1. **Global variable multiple-definitions** — ~40 TUs each *define* (not declare) overlapping globals: `uint gInterpolationOrder;` appears as a definition in ≥ 25 TUs, `bool gPrintReferenceValues = false;` in ≥ 35, `gTestCaseIndex`/`gTestIndex`/`gCaseIndex`, `gDim`, etc. Merged as-is → hard multiple-definition link errors (C++ has no common symbols; GCC ≥ 10 is `-fno-common` anyway).
2. **`extern "C"` helper functions** — 80 `extern "C"` helpers across TUs with massive name reuse: `check_results` ×27, `check_linear_results` ×11, `check_linear_results_serial` ×11, `check_results_serial` ×6, `check_results_parallel` ×5, `check_quadratic_results{,_serial}` ×5 each, plus stragglers. `extern "C"` ignores namespaces for linkage, so **namespace wrapping alone cannot fix these** — the `extern "C"` must be removed. Verified safe: no deck source references any `check_*` symbol; the `extern "C"` is template cargo-cult from `example_test_case.cpp.template`.
3. **Duplicate `TEST_CASE` names** — Catch2 v2.13.9 (vendored `moris/include/catch.hpp`) hard-fails on duplicates: `enforceNoDuplicateTestCases` at `catch.hpp:14278`. Live collisions (3 pairs):
   - `"Channel_with_Four_Cylinders_Static_Linear"` and `"…_Quadratic"`: `thermal/advection/Channel_with_Four_Cylinders_Static/example_test_case.cpp:176,217` **vs** `thermal/diffusion/Channel_with_Four_Cylinders_Static_Temp_Only/example_test_case.cpp:369,462`
   - `"Field_example_write"`: `structure/linear/Beam_Temperature_Field/example_test_case.cpp:112` **vs** `thermal/diffusion/Field_Example/example_test_case.cpp:175` (copy-paste artifact)
   - Two more pairs come only from the **stale, never-built copies** `optimization/LevelSet_Boxbeam_Adaptive_Refinement/src/` and `optimization/Level_Set_Beam_SIMP_Hole_Seeding/src/` (no `add_subdirectory` references them) — exclude, don't migrate.
4. **File-scope initializers that become load-bearing** — `Stefans_Problem_cut/example_test_case.cpp:25` (`uint gInterpolationOrder = 1;`), `Stefans_Problem_conform:25,28`, `Single_Phase_Hollow_Cylinder_Static:18`, `Channel_2D_Static:18-19` (`gInletVelocityBCFlag = true; gInletPressureBCFlag = false;`). Today the fresh process guarantees these values; with shared zero-initialized globals, a deck reading `gInterpolationOrder` would see `0` unless the TEST_CASE sets it. **Migration recipe rule R4 makes every TEST_CASE set every global it or its deck consumes.**
5. **What is *not* a problem:** `int fn_WRK_Workflow_Main_Interface(int, char*[]);` appears in every TU but is a declaration (defined once in WRK-lib). Catch2 registration statics are per-TU and anonymous. `#ifdef MORIS_HAVE_SLEPC`-gated cases compile identically in the merged TU. Multiple workflow invocations per process are already exercised today (`Solver_Examples` runs 4 `fn_WRK_Workflow_Main_Interface` calls in one process; the Restart example runs 2 with different `.so`s), so the "many solves in one process" pattern is proven, and the runner never widens it: each ctest entry still selects exactly one example's cases.

### 2.3 Chosen collision strategy: shared-globals TU + per-TU namespace wrap

**Chosen:**
- **(a) `EXA_Globals.{hpp,cpp}`** — one TU in the runner defines the deck-visible union from §2.1 at global scope (preserving the dlopen ABI bit-for-bit); a header provides `extern` declarations. Test TUs delete their own definitions and include the header.
- **(b) `namespace exa_<example> { … }`** around each TU's remaining contents (helpers, TU-local state), with all `extern "C"` on helpers removed. Two-line mechanical edit per TU; `TEST_CASE` works inside namespaces (the test *name* is a string and stays global).
- **(c) Rename the 3 duplicate TEST_CASE-name pairs** with an example prefix (e.g. `"Temp_Only_Channel_with_Four_Cylinders_Static_Linear"`, `"Beam_Temperature_Field_write"`). ctest names are unaffected; only direct-binary-invocation filters change (owner sign-off item, §6).

**Rejected:**
- *`-D`-based renaming*: cannot fix `extern "C"` collisions differently than deleting them anyway, breaks the dlopen ABI for deck-visible globals unless exceptions are hand-tracked per TU, makes debugging/grep miserable. Strictly dominated.
- *Generated registry* (codegen): Catch2 already **is** the registry; codegen adds a build step and a new custom seam, and does nothing for the dlopen ABI problem.
- *Pure per-TU namespaces without the shared-globals TU*: silently breaks every deck (§2.1) — unresolved-symbol at dlopen or, worse with lazy binding, reading a wrong symbol. This is the trap the whole design must avoid.

### 2.4 The mechanical per-TU migration recipe (applied once per leaf)

- **R1** Delete file-scope global definitions that appear in `EXA_Globals.hpp`; add `#include "EXA_Globals.hpp"`.
- **R2** Remove `extern "C"` from all helper functions.
- **R3** Wrap everything below the includes in `namespace exa_<example_dir_lowercase> { … }`.
- **R4** For every deleted initializer and every global the TU's deck externs, add explicit assignments at the top of **each** `TEST_CASE` (temporal isolation makes same-name reuse across examples safe — each ctest process runs one example's cases only — but *within* the process nothing may rely on zero-init or file-scope init anymore).
- **R5** Append the canonical tag `[EXA_<ExampleDirName>]` to every `TEST_CASE` tag string.
- **R6** If the TU is one of the 3 colliding pairs, rename the Catch2 test name per §2.3(c).

---

## 3. Runner design

### 3.1 `runner_main.cpp`

A near-copy of `projects/EXA/src/example_main.cpp` (which stays in place for unmigrated leaves during coexistence), plus a guard: refusing to run with no arguments prevents the new failure mode of "all 99 test cases execute sequentially in one working directory" (each case loads `./<X>.so` relative to CWD).

```cpp
// projects/EXA/runner/runner_main.cpp   (new file)
#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include <cstring>
#include "cl_Communication_Manager.hpp"
#include "cl_Logger.hpp"
#include "cl_Performance_Reporter.hpp"
#include "banner.hpp"

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

int main( int argc, char* argv[] )
{
    gMorisComm = moris::Comm_Manager( &argc, &argv );
    gLogger.initialize( 2 );
    moris::print_banner( argc, argv );

    // Guard: every example expects to run in ITS working directory with a
    // Catch2 test spec (ctest supplies both). A bare invocation would run
    // all examples' cases in one directory - refuse instead.
    if ( argc < 2 )
    {
        if ( moris::par_rank() == 0 )
        {
            std::cerr << "EXA-test: pass a Catch2 test spec, e.g.  EXA-test.exe \"[EXA_Laplace_2D]\"\n"
                      << "          (run with --list-tests to see all; ctest supplies the spec and working dir)\n";
        }
        gMorisComm.finalize();
        return 1;
    }

    int tRet = Catch::Session().run( argc, argv );

    gPerfReporter.finalize();
    gMorisComm.finalize();
    return tRet;
}
```

### 3.2 `EXA_Globals.{hpp,cpp}`

```cpp
// projects/EXA/runner/EXA_Globals.hpp   (new file)
// Deck-visible globals: input .so files resolve these BY NAME against the
// executable's -rdynamic symbol table at dlopen time. They must stay at
// global scope with exactly these names/types, defined once (EXA_Globals.cpp).
// Regenerate the union with:
//   grep -rhP '^\s*extern\s+(?!"C")' projects/EXA --include='*.cpp' | grep -v example_test_case | sort -u
#pragma once
#include <string>
#include "moris_typedefs.hpp"    // MRS/COR/src

extern moris::uint  gInterpolationOrder;
extern moris::uint  gTestCaseIndex;
extern moris::uint  gTestIndex;
extern moris::uint  gCaseIndex;
extern moris::uint  gDim;
extern moris::uint  gOrder;
extern moris::uint  tGeoModel;
extern moris::uint  tDim;
extern moris::uint  gLevelSetInterpolationOrder;
extern moris::uint  gFEMInterpolationOrder;
extern moris::uint  gLagrMeshInterpolationOrder;

extern bool         gPrintReferenceValues;
extern bool         gUseBspline;
extern bool         gUseMixedTimeElements;
extern bool         gUseBelosWithILUT;
extern bool         gInletVelocityBCFlag;
extern bool         gInletPressureBCFlag;
extern bool         gHaveStaggeredFA;
extern bool         gHaveStaggeredSA;

extern std::string  tOrder;
extern std::string  tStressType;
extern std::string  tOutputFileName;
extern std::string  gPrecSolver;
```

`EXA_Globals.cpp` defines each with no meaningful initializer (rule R4 forbids relying on it). The comment block is the contract; the regeneration grep is the enforcement tool.

### 3.3 Test selection and ordering

- **Mechanism:** Catch2 v2 positional test-spec — the canonical tag from rule R5. ctest COMMAND: `mpirun -n <P> $<TARGET_FILE:EXA-test> "[EXA_Laplace_2D]"`. Tag-based (not name-wildcard) because several TUs contain multiple cases whose names share no prefix (`Solver_Examples`: `Standard_Monolithic`, `Staggered_FA_and_SA`, …), and because a fresh `EXA_`-prefixed tag is guaranteed collision-free.
- **Ordering:** Catch2 v2 default run order is declaration order (`catch.hpp:5300`, `RunTests::InDeclarationOrder`). All intra-example ordering dependencies live inside single TUs — e.g. `LevelSet_Boxbeam_Adaptive_Refinement_Restart/example_test_case.cpp` runs the create-file case (produces `ADV_Alg_0_Iter_11.hdf5`) before the restart case (deck `:525` reads that hdf5). Declaration order within a TU is preserved under any filter. There are **no cross-example dependencies**: no `set_tests_properties`, `DEPENDS`, or `FIXTURES` exist anywhere under EXA (grep verified), and each example's files live in its own working directory.
- **argv passthrough:** none needed — each TEST_CASE constructs its own deck argv internally (the `.so` path, or `--meshgen Bear_Example.xml` for `mesh_generation/image_bear/example_test_case.cpp:59-66`, which also self-stages its data files into CWD). The Catch2 spec is consumed by `Session::run` and MPI args are stripped by `Comm_Manager(&argc,&argv)` before Catch sees them — the same chain `example_main.cpp` uses today.

### 3.4 Working directories, MPI, valgrind — all preserved verbatim

- **WORKING_DIRECTORY** stays `${CMAKE_CURRENT_BINARY_DIR}/bin` per leaf (today: `Laplace_2D/CMakeLists.txt:49,102`). That is where `dynamic_link_input`'s `.so` lands and where every mesh/exodus/hdf5 artifact is read/written — tests run `./<X>.so` relative to it. The runner binary lives in `projects/EXA/runner/bin`, referenced absolutely via `$<TARGET_FILE:EXA-test>`.
- **MPI:** `${MORIS_EXECUTE_COMMAND}` (`mpirun`, root `CMakeLists.txt:144`) `-n ${PROCS}` unchanged; `par_size()`-switching inside test cases unchanged.
- **Valgrind:** `${VALGRIND} ${VALGRIND_OPTIONS_EXA}` (root `CMakeLists.txt:1126-1133`) interpolated into the COMMAND exactly as today.
- **ctest -j:** concurrent examples become concurrent processes of the same binary in different working directories — no shared mutable state, no change.

---

## 4. `moris_add_example()` and a worked leaf conversion

### 4.1 Why a per-leaf function, not a central manifest

A central manifest was considered and rejected: (i) per-example knowledge (PROCS, irregular test names like `SIMP`/`Comsol_conform`/`Mach_Leading_Edge`/`Fick_Problem`/`single_element`; multiple input decks; the 5 leaves whose `SO_INCLUDES` deviate from the standard block) stays next to the example it describes; (ii) migration diffs are leaf-local and independently revertable; (iii) coexistence is trivial — an unmigrated leaf simply keeps its old file. The one thing that must be central — the source-list aggregation — is a two-line global-property pattern.

### 4.2 Function draft

```cmake
# share/cmake/utilities/moris_add_example.cmake   (new file)
#
# moris_add_example(
#     NAME <ExampleDirId>            # unique; Catch2 tag becomes [EXA_<NAME>]
#     TEST_BASE <name>               # ctest base name == old EXAMPLE_FILE (e.g. Laplace)
#     PROCS <n> [<n> ...]            # registers <TEST_BASE>-<n>-procs per entry
#     [NO_PROCS_SUFFIX]              # single test named exactly <TEST_BASE> (SIMP, Comsol_conform, ...)
#     [SOURCES <files> ...]          # default: example_test_case.cpp
#     [INPUTS <deckbase> ...]        # each: dynamic_link_input(<b> <b> <b>.cpp SO_INCLUDES); default: TEST_BASE
#     [EXTRA_SO_INCLUDES <dirs>...]  # for the 5 leaves whose SO_INCLUDES deviate from the standard block
# )
#
# Preserves, byte-for-byte, the dynamic_link_input() targets, their output
# names/locations, and the ctest names/PROCS of the old per-leaf boilerplate.

function(moris_add_example)
    set(options NO_PROCS_SUFFIX)
    set(oneValueArgs NAME TEST_BASE)
    set(multiValueArgs PROCS SOURCES INPUTS EXTRA_SO_INCLUDES)
    cmake_parse_arguments(EXA "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    if(NOT EXA_SOURCES)
        set(EXA_SOURCES example_test_case.cpp)
    endif()
    if(NOT EXA_INPUTS)
        set(EXA_INPUTS ${EXA_TEST_BASE})
    endif()

    # ---- input .so machinery: identical to the old leaf boilerplate --------
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

    set(SO_TPLS "trilinos" ${ARMADILLO_EIGEN})
    get_property(INT_SRC_LIST GLOBAL PROPERTY INT_SRC_LIST)
    get_property(MTK_SRC_LIST GLOBAL PROPERTY MTK_SRC_LIST)
    set(SO_INCLUDES                                          # == the block 60/65 leaves share
        ${MORIS_BOOST_INCLUDE_DIRS}
        ${MORIS_TRILINOS_INCLUDE_DIRS}
        ${MORIS_PACKAGE_DIR}/ALG/src
        ${INT_SRC_LIST}
        ${MTK_SRC_LIST}
        ${MORIS_PACKAGE_DIR}/FEM/MSI/src
        ${MORIS_PACKAGE_DIR}/FEM/VIS/src
        ${MORIS_PACKAGE_DIR}/SOL/DLA/src
        ${MORIS_PACKAGE_DIR}/SOL/TSA/src
        ${MORIS_PACKAGE_DIR}/SOL/NLA/src
        ${MORIS_PACKAGE_DIR}/SOL/SOL_CORE/src
        ${MORIS_PACKAGE_DIR}/LINALG/src
        ${MORIS_PACKAGE_DIR}/LINALG/src/${LINALG_IMPLEMENTATION_INCLUDES}
        ${MORIS_PACKAGE_DIR}/COM/src
        ${MORIS_PACKAGE_DIR}/PRM/ENM/src
        ${MORIS_PACKAGE_DIR}/MTK/src
        ${MORIS_PACKAGE_DIR}/HMR/src
        ${MORIS_PACKAGE_DIR}/XTK/src
        ${MORIS_PACKAGE_DIR}/PRM/src
        ${MORIS_PACKAGE_DIR}/MRS/COR/src
        ${EXA_EXTRA_SO_INCLUDES})
    foreach(TPL ${SO_TPLS})
        string(TOUPPER ${TPL} TPL)
        list(APPEND SO_INCLUDES ${MORIS_${TPL}_INCLUDE_DIRS})
    endforeach()

    foreach(INPUT ${EXA_INPUTS})
        dynamic_link_input(${INPUT} ${INPUT} ${INPUT}.cpp ${SO_INCLUDES})
    endforeach()

    # ---- register this leaf's test-case TU(s) with the shared runner -------
    foreach(SRC ${EXA_SOURCES})
        get_filename_component(SRC_ABS ${SRC} ABSOLUTE)
        set_property(GLOBAL APPEND PROPERTY EXA_RUNNER_SOURCES ${SRC_ABS})
    endforeach()

    # ---- ctest registrations: EXACT legacy names, PROCS, workdir, valgrind --
    if(MORIS_HAVE_PARALLEL_TESTS)
        foreach(PROCS ${EXA_PROCS})
            if(EXA_NO_PROCS_SUFFIX)
                set(_test_name ${EXA_TEST_BASE})
            else()
                set(_test_name ${EXA_TEST_BASE}-${PROCS}-procs)
            endif()
            add_test(NAME ${_test_name}
                COMMAND ${MORIS_EXECUTE_COMMAND} -n ${PROCS} ${VALGRIND} ${VALGRIND_OPTIONS_EXA}
                        $<TARGET_FILE:EXA-test> "[EXA_${EXA_NAME}]"
                WORKING_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
        endforeach()
    endif()
endfunction()
```

And the aggregation directory, added **last** in `projects/EXA/CMakeLists.txt` (the `$<TARGET_FILE:EXA-test>` generator expression in earlier-processed leaves resolves at generate time, so definition order within configure is fine):

```cmake
# projects/EXA/runner/CMakeLists.txt   (new file)
get_property(EXA_RUNNER_SOURCES GLOBAL PROPERTY EXA_RUNNER_SOURCES)
if(NOT EXA_RUNNER_SOURCES)
    return()   # Step 0 / nothing migrated yet: no runner, no cost
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${BIN})

add_executable(EXA-test
    runner_main.cpp
    EXA_Globals.cpp
    ${EXA_RUNNER_SOURCES})
target_include_directories(EXA-test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})   # EXA_Globals.hpp
target_link_libraries(EXA-test PRIVATE
    ${HMR}-lib ${COM}-lib ${WRK}-lib
    ${MORIS_PETSC_LIBS} ${MORIS_BOOST_LIBS} ${MORIS_ACML_LAPACK_MKL_OPENBLAS_LIBS}
    ${MORIS_ARMADILLO_EIGEN_LIBS} ${MORIS_SUPERLU_LIBS} ${MORIS_LDLIBS}
    ${MORIS_TRILINOS_LIBS} ${MORIS_BASE_LIBS})     # == the identical list all 65 leaves use today
set_target_properties(EXA-test PROPERTIES OUTPUT_NAME EXA-test${EXE_EXT})
```

```diff
--- a/projects/EXA/CMakeLists.txt
+++ b/projects/EXA/CMakeLists.txt
@@
+include(${MORIS_CMAKE_DIR}/utilities/moris_add_example.cmake)
 add_subdirectory(thermal)
 add_subdirectory(structure)
 add_subdirectory(fluid)
 add_subdirectory(optimization)
 add_subdirectory(mesh_generation)
+add_subdirectory(runner)    # must stay last: consumes EXA_RUNNER_SOURCES
```

Note the test-case TUs are listed **directly in `add_executable`** — never moved into an intermediate static library, which would let the linker drop Catch2's self-registration objects.

### 4.3 Worked conversion: `thermal/diffusion/Laplace_2D`

`CMakeLists.txt` — 105 lines → one call (the deleted `EXAMPLE_INCLUDES` loop at `:29-31` was already inert: `include_directories(../COM/src)` resolves relative to the leaf to a nonexistent path; the real headers come from the global `include_directories` at root `:985`):

```cmake
moris_add_example(
    NAME       Laplace_2D
    TEST_BASE  Laplace          # ctest names stay Laplace-1-procs / Laplace-4-procs
    PROCS      1 4
    # INPUTS defaults to TEST_BASE: dynamic_link_input(Laplace Laplace Laplace.cpp ...) as before
)
```

`example_test_case.cpp` — recipe R1–R5 applied (abridged):

```diff
 #include <catch.hpp>
 ...
 #include "fn_norm.hpp"
+#include "EXA_Globals.hpp"                // shared deck-visible globals (dlopen ABI)

 using namespace moris;

-// global variable for interpolation order
-uint gInterpolationOrder;
-
-// flag to print reference values
-bool gPrintReferenceValues = false;
-
-uint gTestIndex;
+// R3: everything below is TU-local; names may repeat across examples.
+namespace exa_laplace_2d
+{

 int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

-extern "C" void
+void                                       // R2: extern "C" removed (no deck references helpers)
 check_results(
         const std::string &aExoFileName,
         uint               aTestCaseIndex,
         bool               aSolvewithPetsc = false )
@@
 TEST_CASE( "Laplace_Anasazi",
-        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Anasazi]" )
+        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Anasazi],[EXA_Laplace_2D]" )   // R5
 {
     int argc = 2;

     gTestIndex          = 0;
     gInterpolationOrder = 1;
+    gPrintReferenceValues = false;         // R4: no reliance on file-scope/zero init
@@
+}    // namespace exa_laplace_2d
```

The generated ctest entry is then byte-comparable to today's: same names `Laplace-1-procs`/`Laplace-4-procs`, same `mpirun -n N`, same `WORKING_DIRECTORY …/Laplace_2D/bin`, command now `…/EXA/runner/bin/EXA-test.exe "[EXA_Laplace_2D]"`.

---

## 5. Migration plan

### 5.1 Steps (each leaves the build green and is independently revertable)

- **Step 0 — infrastructure, zero behavior change.** Add `moris_add_example.cmake`, `runner/` (main, globals, CMakeLists with the empty-property early-return), the `include()`/`add_subdirectory(runner)` hook. No leaf migrated → no runner target, no new edges. *Verify:* `build_metrics.py --diff` vs baseline ≈ 0; `ctest -N` unchanged.
- **Step 1 — pilot: `Laplace_2D`** (the worked diff above; it exercises multi-PROCS, multi-TEST_CASE, `#ifdef`-gated cases, and deck-visible globals). This step *adds* the first runner link (+~150 s, +~35 MB) while removing one example link — metrics ≈ neutral by design; the point is correctness. *Verify per §5.3.*
- **Steps 2–6 — batches by physics directory**, smallest/simplest first, ordered so both TUs of a duplicate-name pair land together:
  1. `thermal/` (16 leaves; includes the `Channel_with_Four_Cylinders_Static_Linear/_Quadratic` rename pair — migrate both affected leaves in the same PR),
  2. `fluid/` (6),
  3. `structure/` (13; includes the `Field_example_write` rename in `Beam_Temperature_Field`, the two `NO_PROCS_SUFFIX` cases, and `Column_Buckling`'s string globals),
  4. `mesh_generation/` (3; no deck `.so` — `INPUTS` empty, `--meshgen` argv),
  5. `optimization/` last (24; the biggest and quirkiest: irregular test names, the Restart two-`.so` leaf, the 4+1 leaves with deviant `SO_INCLUDES` — pass `EXTRA_SO_INCLUDES`).
- **Step 7 — cleanup (separate PR, owner sign-off):** delete `projects/EXA/src/example_main.cpp` (now unused), the two stale `src/` copies, and the per-leaf boilerplate remnants; run the full before/after measurement (§1.5) and archive the new labeled snapshot in `$MORIS_RUNS_DIR/benchmarks/build_metrics/`.

Each batch removes N example links and leaves exactly one runner link — the metrics improvement accrues monotonically and `--diff` after each batch should show `link/count` down by N and `executables/total_mb` down by ≈ 34·N MB.

### 5.2 Coexistence mechanics

Unmigrated leaves are untouched — they still build `Main_<X>.exe` from `../../../src/example_main.cpp` and register their tests against it. Migrated leaves contribute their TU to the runner. Both patterns coexist in one configure; ctest names never collide (each name exists exactly once, on exactly one side of the migration at any commit).

### 5.3 Per-leaf verification (run by moris-build, in a scratch build dir)

```bash
# 1) ctest name/count parity — the hard constraint:
ctest -N | grep -E "Test *#" | sed 's/.*: //' | sort > after.txt   # compare vs the pre-migration list
diff before.txt after.txt                                          # must be empty at every step
# (reference list: 73 names, captured 2026-07-01 via
#  grep -rh "^add_test" build/projects/EXA --include=CTestTestfile.cmake | sed 's/add_test(//;s/ .*//' | sort)

# 2) Catch2 case parity for the migrated leaf:
old/bin/Main_Laplace.exe --list-test-names-only                       > old_cases.txt
new/runner/bin/EXA-test.exe "[EXA_Laplace_2D]" --list-test-names-only > new_cases.txt
diff old_cases.txt new_cases.txt        # identical names, identical order (decl order)

# 3) Behavior: the leaf's tests pass at every PROCS value:
ctest -R "^Laplace-" --output-on-failure

# 4) The pymoris seam is untouched: the leaf's .so edges in build.ninja are unchanged:
grep -A3 "Laplace.so" before/build.ninja > a; grep -A3 "Laplace.so" after/build.ninja > b; diff a b

# 5) Metrics: build_metrics.py --diff shows only the expected -1 exe / -1 link edge.
```

Known-failing tests (`Channel_with_Four_Cylinders_Static*` per `.claude/rules/testing-overview.md`) are compared against their *current* status, not against green.

---

## 6. Risks, open questions, handoffs

### Risks

1. **dlopen ABI regressions** — the whole design hinges on §2.1's union being complete. Mitigation: the regeneration grep is in `EXA_Globals.hpp` itself; per-leaf verification actually runs the deck, which fails loudly on an unresolved symbol.
2. **Initializer semantics (R4)** — the known instances are listed in §2.2(4), but a missed one produces a *silently different simulation* (e.g. interpolation order 0). Mitigation: R4 is applied to every TU, not just known offenders; each leaf's tests run against unchanged reference values, which is exactly what these tests check.
3. **State leakage across TEST_CASEs within one example** — unchanged from today (each exe already ran all its cases in one process; the tag filter selects the same set), but the *bare* runner invocation is a new footgun — mitigated by the `argc < 2` guard.
4. **LTO memory/time for the one link** — same object set plus ~66 small TUs; `-flto=4` partitioning unchanged. If the runner link exceeds ~200 s or RAM becomes an issue, the fallback is `-fno-lto` on the runner target only (examples are correctness tests, not performance baselines) — estimated 20–40 s link, but changes codegen vs today, so proposed only as an opt-in follow-up measurement.
5. **Catch2 v3 migration someday** — this design neither helps nor hurts (one runner main to update instead of 57 copies is strictly easier).

### Open questions for the owner

1. **Catch2-level test-name renames (§2.3c)** — anyone invoking old binaries directly by test name is affected; ctest names are not. Approve the three renames?
2. **The 8 disconnected examples** (`Homogenization_2D/3D`, `Plane_Strain`, `Channel_2D_Compressible`, `Heated_Bubble_2D`, `Channel_with_Four_Cylinders_Transient`, `Two_Channels_with_Separation_Wall_Transient`, `Thermal_Flow_About_Sphere_3D` — have `CMakeLists.txt` but no `add_subdirectory` reaches them) — reconnect during migration (~1.5 s compile each, no extra link) or delete? Reconnecting *adds* ctest entries, which technically violates "exact current set" — needs an explicit decision.
3. **The 2 stale `src/` copies** (`LevelSet_Boxbeam_Adaptive_Refinement/src/`, `Level_Set_Beam_SIMP_Hole_Seeding/src/`) — delete in Step 7?
4. **`build_dbg` policy** — with −20 GB on the table, should dbg keep building examples at all, or is `MORIS_USE_EXAMPLES=OFF` for local dbg (CI keeps it ON) the better default? Orthogonal to this RFC but the numbers belong together.
5. **Correctness-only rider (unquantified, flagged per policy):** the `example_test_case.cpp.template` under `projects/EXA/templates/` should be updated to the post-migration pattern in the same series, or the next new example reintroduces the collisions. Moves no metric; bundled here only as a flag.

### Handoffs

- **Implementer (main session):** apply Step 0 + Step 1 patches from §3/§4 verbatim; then batch PRs per §5.1 using the §2.4 recipe.
- **moris-build:** run the §1.5 measurement protocol (before on `main`, after on the branch, scratch dirs only) and the §5.3 per-leaf verification; owns the env (`spack env activate /home/doble/codes && source ~/.bashrc_moris`).
- **pymoris seam — flagged loudly:** nothing here touches `dynamic_link_input` target names, `.so` output paths, or the `shared_object_file`/`create_shared_object.sh` machinery that pymoris scrapes from `build.ninja`; §5.3 step 4 proves it per leaf. If any future variation of this RFC proposes centralizing `SO_INCLUDES` beyond the byte-identical reproduction in `moris_add_example()`, that crosses the seam and needs a pymoris-side check first.
- **Docs debt (not this RFC):** `share/doc` example-authoring documentation and `.claude/skills/build-moris/` must be updated post-migration; both already contain stale claims tracked in `CMAKE_REVIEW.md`.
