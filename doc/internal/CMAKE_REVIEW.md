# MORIS CMake Build System — Full Review

**Date:** 2026-07-01 · **Tree:** `/home/doble/moris_workspace/moris` (branch `main`, clean) · **Reviewer:** cmake-expert agent · **Verified against:** live caches, `build.ninja`, `.ninja_log`, and two scratch probes. Prior audit `CMAKE_AUDIT.md` (2025-09-28, 52/100) re-checked: its structural findings all still hold, though its line numbers have drifted (e.g. MRS include is now `:1003` not `:1008`; `CMAKE_BUILD_TYPE` force is at `:481`/`:512` not `:736-770`). It remains untracked and unactioned.

## Verdict

This build system works — reliably, on exactly one machine, in exactly one shell environment, and nowhere else. It is a 2015-era recursive-include architecture with a genuinely thoughtful packaging layer bolted on, then years of copy-paste drift. The two things that actually cost you are: (1) a 58-minute full build in which **74% of all CPU time is spent statically re-linking the same code 81 times under `-flto=4`**, and (2) a configure step whose output depends on `LD_LIBRARY_PATH`, env vars, and set-compiler-after-`project()` hacks — which has already silently produced **two live build trees on different compilers** (`build_opt` = system `/usr/bin/g++`, `build_dbg` = Spack `mpic++`). Beyond that there is a family of confirmed typo-bugs that are masked by the global-state architecture — the globals paper over the local mistakes, which is exactly why both are hard to remove. My honest position: the *package layout and export concept are worth keeping*; the *toolchain/flag/configure layer is what I would refuse to build anything new on top of* until it is pinned by presets. None of this requires a rewrite; the top three fixes below are an afternoon.

---

## 1. What is genuinely good (do not refactor this away)

- **The hand-rolled component/export system is ahead of its time for a research code.** Per-package `Config.cmake` + `Targets.cmake`, a component-aware top-level `morisConfig.cmake.in` producing `moris::COMPONENT` imported targets (`share/cmake/morisConfig.cmake.in`), and `configure_package_config_file`/`write_basic_package_version_file` used correctly in every package (e.g. `projects/HMR/src/CMakeLists.txt:186-221`). Most academic codes have nothing like this.
- **The `_new` TPL wrappers create namespaced `INTERFACE IMPORTED` targets** (`moris::gperftools` etc., e.g. `find_wrappers/gperftools_new.cmake:27-31`) — that's the modern pattern, done before it was fashionable.
- **Package libraries are half-modern already**: `target_include_directories(... PUBLIC $<BUILD_INTERFACE:...> $<INSTALL_INTERFACE:...>)` and `target_link_libraries(PUBLIC ...)` (`projects/HMR/src/CMakeLists.txt:141-144`). The migration path is short.
- **`snippets/` is not cruft — it's docs-as-tests, and it's load-bearing.** The 138 `.inc` doxygen example snippets are `#include`d by real unit tests (`projects/MRS/IOS/test/fn_to_stdio.cpp:22`) and fed to Doxygen via `EXAMPLE_PATH`/`INCLUDE_PATH` (`share/doc/Doxyfile.in:963-964,2186`). That's why `"snippets"` is on the include path at root `:130`. Anyone "cleaning up" that directory breaks tests. (The Doxyfile itself is in decent shape: `INPUT`, `IMAGE_PATH`, `HTML_EXTRA_FILES` all point at files that exist.)
- **`test_includes.cmake` is textbook**: a single `test-libs` INTERFACE target (`share/cmake/dependencies/test_includes.cmake:16-17`).
- **Selective package builds actually work.** The despised sentinel/depends graph does deliver: `build_dbg` builds without EXA (`MORIS_USE_EXAMPLES:BOOL=OFF` in its cache) and is compile-dominated, not link-dominated.
- **The in-source-build guard, out-of-source discipline, and a simple honest CI** (`.github/workflows/ctest.yml` — raw flags, dbg + opt, `ctest --output-on-failure`).
- **The doxygen wiring is correct** (`share/doc/CMakeLists.txt:12-24` — `configure_file` + `doxy` custom target). It is merely inert on this machine: `DOXYGEN_EXECUTABLE-NOTFOUND` in `build_opt/CMakeCache.txt`, so no `doxy` edge exists in `build.ninja`. That's a missing tool, not broken CMake.

---

## 2. Findings

### A. Bugs (wrong now)

**A1 — Two live trees, two compilers; compiler selection is out of CMake's control.**
Root `:344-358` defaults `MORIS_C(XX)_COMPILER` to `mpicc`/`mpic++`, then root `:955-956` does `set(CMAKE_C(XX)_COMPILER ...)` — *after* `project()` at `:24`, which official CMake docs say is unsupported (set once, at first configure). The result on disk: `build_opt/CMakeFiles/3.31.8/CMakeCXXCompiler.cmake` says `/usr/bin/g++`; `build_dbg/CMakeFiles/3.31.8/CMakeCXXCompiler.cmake` says `/home/doble/codes/.spack-env/view/bin/mpic++`. So Release binaries and Debug binaries are built by *different toolchains* — system GCC linking Spack-built Boost 1.88/HDF5 in one tree, Spack MPI wrapper in the other. Failure scenario: "works in dbg, crashes in opt" bugs that are ABI artifacts, not code; and any fresh configure lands on whatever compiler the shell happens to expose. This is the single most disqualifying property of the current setup.

**A2 — `ctest -L` label filtering is fiction, and the docs teach it.**
Zero `LABELS` properties exist anywhere: 0 hits across all 182 `CMakeLists.txt` and all of `share/cmake` (94 `add_test` calls, 23 `set_tests_properties` — none set labels). Yet `moris/CLAUDE.md` ("`ctest -L XTK` — by label") and the `run-ctest`/`test-moris` skills document label filtering. `ctest -L XTK` runs **zero tests and exits 0** — the "all green" that validates nothing. This one can actively cost correctness. (Fix is one line; see Patch 2 — verified by scratch probe that a directory-level `LABELS` property inherits into nested `test/` subdirectories.)

**A3 — The typo'd-variable family (option set, guard tests a different name):**
- `option(MORIS_USE_SACADO ...)` but `if (USE_SACADO)` — root `:702-703`. Enabling Sacado does nothing, silently.
- `if(USE_GPERFTOOLS)` vs `MORIS_USE_GPERFTOOLS` — `share/cmake/dependencies/MRS_includes.cmake:46`, `LINALG_Depends.cmake:42`. The proper imported-target link path (`moris::gperftools`) is dead code. gperftools only works in `build_opt` because root `:674` does a *global raw* `link_libraries("-Wl,--no-as-needed -lprofiler -Wl,--as-needed -ltcmalloc")` with a hardcoded `-I/usr/include/gperftools/` (root `:672`) — i.e. the **system** gperftools is silently used, bypassing Spack's, on every one of the 153 link lines that carry `-ltcmalloc` in `build_opt/build.ninja`.
- `if (HAVE_PEDANTIC)` at root `:352` (Intel branch) tests a variable that doesn't exist yet and isn't the option (`MORIS_HAVE_PEDANTIC`). Dead in practice.

**A4 — The pedantic block reuses one cached check variable** (root `:565-584`): `check_cxx_compiler_flag` caches `HAVE_PEDANTIC` on the first call (`-pedantic`), so the checks for `-pedantic-errors`, `-Wall`, `-Werror` are all skipped and the flags appended untested (confirmed in the live `build_opt` FLAGS line: `-pedantic -pedantic-errors -Wall -Werror`). Harmless on GCC today; blind on any other compiler. Related: `MORIS_HAVE_PEDANTIC` defaults **ON** (root `:182-183`), so `-Werror` is the everyday default — contrary to the canon (opt-in only) and contrary to what the `ci-local` skill claims ("`-Werror` is enabled only in CI preset"). Every GCC upgrade will break the build at some random warning.

**A5 — `STK_Depends.cmake` is a vestigial file with live side effects.**
There is no `projects/STK` package and no `set(STK ...)` anywhere in the root (source-dir block `:53-91` has no STK), so `list(APPEND MORIS_SOURCE_DIRS ${STK})` (`STK_Depends.cmake:19`) appends nothing, and `STK_TPL_DEPENDENCIES` is consumed by *no file in the tree*. But the file **does** run on every `BUILD_MAIN` configure (`main_includes.cmake:43`) and its copy-paste remnant at `:36-39` appends STK's TPLs to **`HMR_TPL_DEPENDENCIES`** — which flows straight into HMR's PUBLIC link/include loop (`projects/HMR/src/CMakeLists.txt:147-158`). HMR links superlu/trilinos through pollution, not declaration. Works today; will confound any TPL surgery.

**A6 — `EXA_Depends.cmake:22` writes `MAIN_TPL_DEPENDENCIES`, not `EXA_TPL_DEPENDENCIES`.** Masked only because `main_includes.cmake:22` sets the identical value behind its sentinel.

**A7 — `dynamic_link_input()` works by triple accident** (`share/cmake/utilities/dynamic_link_input.cmake`): declares parameter `so_includes` but reads `${SO_INCLUDES}` (`:37`) — resolves only because every EXA caller happens to have `SO_INCLUDES` set in caller scope; `target_compile_definitions(... INTERFACE ...)` (`:41`) puts the definitions on *consumers* of the `.so`, not the `.so` itself — masked by the global `add_definitions` at root `:982`; and `SO_LIB_REQS` (`:19-32`) is entirely commented out, so the `.so` links nothing and resolves at `dlopen` against the `-rdynamic` executable. This function is the poster child for "the globals mask the bugs."

**A8 — `FindPETSc.cmake:71`** sets `PETSC_LIBRARIES` from `${PETSC_LIBRARY}` five lines before `:76` defines `PETSC_LIBRARY`. Empty on that path. Still present.

**A9 — Install lists diverge from the build.** `find_wrappers/CMakeLists.txt:8-40` omits `lbfgsb_new.cmake`, `gperftools_new.cmake`, `arborx_new.cmake`; `find_modules/CMakeLists.txt:10-25` omits `FindLAPACK.cmake`, `FindLBFGSB.cmake`, `FindGperftools.cmake`. Since every package Config includes its wrappers by name, an installed `find_package(moris)` fails on any package whose TPL list touches the omitted ones. **Honest severity: near-zero for this workspace** — nothing in the workflow runs `cmake --install`; pymoris consumes the build tree. File under "fix if you ever ship."

**A10 — `create_shared_object.sh` mutates the tracked source tree** (moves `projects/mains/input_file.cpp`, symlinks the user file, hand-rolled lockfile — `share/scripts/create_shared_object.sh:108-134`). New find this pass: line 130 is `rm -f "MORISROOT/$builddir/lib/input_file.so"` — missing `$`, so it deletes a nonexistent relative path and the stale-`.so` cleanup silently never happens (benign only because ninja rebuilds the target anyway).

**A11 — PYBIND is a ghost.** `projects/PYBIND/` is untracked and referenced by *nothing* (0 hits in root, `projects/CMakeLists.txt`, or any depends file). The `moris_py.cpython-313-*.so` in `build_opt/lib` is dated Jan 1 — a stale artifact the current tree cannot rebuild. Anyone "rebuilding and testing" the Python bindings is testing a six-month-old binary.

**A12 — The workspace skills document presets that do not exist.** No `CMakePresets.json` anywhere; `.claude/skills/ci-local/SKILL.md:15` says `cmake --preset moris-ci`, `build-moris/build-rules.md:24` names `moris-ninja-relwithdebinfo`. The real CI (`.github/workflows/ctest.yml:20,36`) uses raw `-D` flags. The docs are not just stale — they describe a configuration mechanism that was never built.

### B. Cost / fragility (works, but bleeds time)

**B1 — The link tail is the build.** From `build_opt/.ninja_log` (verified this pass): one contiguous session, **57.9 min wall**, 2540 edges, 396 CPU-min total; **553 link edges = 293.6 CPU-min = 74% of all edge time**; the top 12 edges are all EXA example links at **147–151 s each**. 81 executables (57 EXA + 24 test) each statically link the full PUBLIC closure under a raw `-flto=4` (root `:524` — hand-rolled, not `CheckIPOSupported`/`INTERPROCEDURAL_OPTIMIZATION`, so LTO parallelism is hard-capped at 4 regardless of `-j`). Everything is STATIC (`MORIS_HAVE_SHARED` OFF, `:147-148`; `LIB_LINK_MODE`, `:656-659`) and everything is PUBLIC, so the closure never shrinks.

**B2 — No ccache, anywhere.** Not installed on the machine, no `CMAKE_CXX_COMPILER_LAUNCHER` in either cache. `build_dbg` is **36 GB**.

**B3 — Configure output is a function of the shell.** `LIB_SEARCH_PATH` is built from `$LD_LIBRARY_PATH` + `$SPACK_LINK_DIRS` (root `:377-379`), then eight `find_library` + global `link_libraries()` calls bake absolute paths to `libgfortran`/`ssl`/`crypto`/`z` into every target (root `:382-455`); `$GFORTLIB` steers which gfortran. Two shells → two different dependency sets, cached forever. Combined with A1, this is why one should not build reproducibility tooling (pymoris fingerprinting) on top of the current configure step.

**B4 — `CMAKE_BUILD_TYPE` is force-overwritten** (root `:481`, `:512` — `FORCE`). `RelWithDebInfo` is unreachable; the real switch is the non-standard `MORIS_HAVE_DEBUG`. Any user/preset intent is silently discarded.

**B5 — Global state everywhere, and it's the masking agent.** `add_definitions`/`include_directories` at root `:982-985`; ~190 global `include_directories` calls tree-wide vs 57 `target_include_directories` (counts re-verified); package files still layering globals on top of their own target-scoped includes (`projects/XTK/src/CMakeLists.txt:105-119`, `projects/GEN/src/CMakeLists.txt:23-28`). Note also the TPL loop puts **`LIBRARY_DIRS` on the include path** (`projects/HMR/src/CMakeLists.txt:154-155`). Every bug in section A7 compiles anyway *because* of these globals — removing them is the refactor, and it must be done package-by-package.

**B6 — `CMAKE_MODULE_PATH` shadows official Find modules** (root `:122-124`): `find_modules/` contains `FindHDF5.cmake`, `FindLAPACK.cmake`, `FindEigen3.cmake`. Any `find_package(HDF5)` in module mode — including ones issued from inside TPL logic — gets MORIS's weaker copy instead of CMake's. Canon violation with real teeth on any machine but this one.

**B7 — `BUILD_ALL` is sticky.** `build_all_exe.cmake` force-sets every `BUILD_*` to ON **with `FORCE` in the cache**; setting `BUILD_ALL=OFF` later leaves them all ON. Cache surgery or a fresh build dir required. (The `build_opt` cache shows the fossil: `BUILD_ALL:BOOL=OFF` with everything still enabled.)

**B8 — `/usr/include` explicitly on the global include path** (root `:130`) — reorders system headers ahead of toolchain paths; benign on this box, classic latent breakage elsewhere.

### C. Style / old-fashioned but fine to leave

- 17 dead old-style wrappers alongside the `_new` set — **0 references anywhere** (verified). Delete at leisure; they cost nothing but confusion.
- `share/cmake/utilities/CMakeLists.txt` is empty (header only) yet still `add_subdirectory`'d at root `:1174`. `linear_algebra_lib_fix_obsolete.cmake` announces its own obsolescence.
- Hand-rolled `-std=c++17` with a *silent* fallback to C++14 (root `:604-609`) instead of `target_compile_features` — on GCC it always takes 17, so it's cosmetic, but the fallback would be a nasty silent downgrade on a hypothetical old compiler.
- `cmake_minimum_required(VERSION 3.17)` bare, no range/policy pin, configured under both 3.28.3 and 3.31.8 (both compiler-ID dirs exist in each build tree). Not urgent; CMake 4.0's `<3.5` removal does not affect a 3.17 floor.
- The cyclic package graph (`projects/CMakeLists.txt:26-30` — the comment openly admits it) — with static libs and PUBLIC links this *works*, and untangling it is a research project of its own. Leave it.
- **Explicitly does not matter for this workspace:** Windows/MSVC support, `GNUInstallDirs`, `NAMESPACE` on the export, SOVERSION, sanitizer matrices, install-tree completeness (A9). Single platform, single user, build-tree consumption. Don't let anyone cargo-cult these onto you.

---

## 3. Recommendations, leverage-ordered

1. **Create `CMakePresets.json` pinning the Spack toolchain + real option sets** (Patch 1). Kills A1 (compiler divergence) at the root for all *future* configures, makes A12 (lying skills) true retroactively, and gives moris-build one command per configuration. ~30 minutes including a fresh scratch verification. Note: existing `build_opt`/`build_dbg` keep their current compilers until *deliberately* reconfigured fresh — that's a moris-build decision, not something to do casually (pymoris fingerprints against those caches).
2. **Make `ctest -L` real** (Patch 2). One line in `projects/CMakeLists.txt`, verified mechanism. Also instantly gives `ctest -LE EXA` to skip the slow example tests.
3. **Attack the link tail: `MORIS_USE_LTO=OFF` + `MORIS_USE_EXAMPLES=OFF` in the dev preset, keep both ON in a `full` preset** (folded into Patch 1). Zero source changes; the ninja log says this alone removes the ~294 CPU-min link tail from daily work — expect full passes in the ~15–25 min range instead of 58. Longer-term (optional): replace the raw `-flto=4` with `CheckIPOSupported` + `INTERPROCEDURAL_OPTIMIZATION` so LTO parallelism scales and per-target opt-out becomes possible.
4. **Install ccache and wire `CMAKE_<LANG>_COMPILER_LAUNCHER` into the presets** (already in Patch 1, commented). Biggest payoff on branch switches and flag churn.
5. **Fix the typo-guard family** (Patch 3: `USE_SACADO`, 2× `USE_GPERFTOOLS`, `:352 HAVE_PEDANTIC`, unique vars in the pedantic block). Ten minutes, removes silent misconfiguration.
6. **Decide PYBIND's fate**: either wire `projects/PYBIND` into the depends graph behind `BUILD_PYBINDINGS` and track it in git, or delete the stale `moris_py*.so` from `build_opt/lib` so nobody trusts it. Needs the owner's call.
7. **Delete `STK_Depends.cmake`** (moving its `MOD`/`TOL`/`DLA` includes into `main_includes.cmake` directly), fix `EXA_Depends`, fix `FindPETSc.cmake` ordering, fix the `rm` line in `create_shared_object.sh`. Mechanical.
8. **Only then** consider the structural refactor (PRIVATE demotion, killing root globals, one `moris_base` INTERFACE target for flags/definitions) — package by package, each step leaving the build green. Do **not** attempt this before 1–4; the presets are what make each step verifiable.

**If you only do 3 things:** (1) Presets with pinned Spack compilers; (2) the one-line test-labels fix; (3) LTO+examples OFF for daily builds. Together: reproducible toolchain, honest test filtering, and roughly two-thirds of the build time back — for well under a day of work.

---

## 4. Draft patches

### Patch 1 — `CMakePresets.json` (new file at repo root)

```json
{
  "version": 3,
  "cmakeMinimumRequired": { "major": 3, "minor": 21 },
  "configurePresets": [
    {
      "name": "base",
      "hidden": true,
      "generator": "Ninja",
      "cacheVariables": {
        "CMAKE_C_COMPILER": "/home/doble/codes/.spack-env/view/bin/mpicc",
        "CMAKE_CXX_COMPILER": "/home/doble/codes/.spack-env/view/bin/mpic++",
        "BUILD_ALL": "ON",
        "MORIS_HAVE_SYMBOLIC": "ON",
        "MORIS_HAVE_SYMBOLIC_STRONG": "OFF"
      }
    },
    {
      "name": "opt", "displayName": "Release, dev (no LTO, no examples)",
      "inherits": "base", "binaryDir": "${sourceDir}/build_opt",
      "cacheVariables": {
        "MORIS_HAVE_DEBUG": "OFF",
        "MORIS_USE_LTO": "OFF",
        "MORIS_USE_EXAMPLES": "OFF"
      }
    },
    {
      "name": "opt-full", "displayName": "Release, full (LTO + examples, CI-equivalent)",
      "inherits": "base", "binaryDir": "${sourceDir}/build_opt_full",
      "cacheVariables": {
        "MORIS_HAVE_DEBUG": "OFF",
        "MORIS_USE_LTO": "ON",
        "MORIS_USE_EXAMPLES": "ON"
      }
    },
    {
      "name": "dbg", "displayName": "Debug",
      "inherits": "base", "binaryDir": "${sourceDir}/build_dbg",
      "cacheVariables": {
        "MORIS_HAVE_DEBUG": "ON",
        "MORIS_USE_EXAMPLES": "OFF"
      }
    }
  ],
  "buildPresets": [
    { "name": "opt", "configurePreset": "opt", "jobs": 4 },
    { "name": "dbg", "configurePreset": "dbg", "jobs": 4 }
  ],
  "testPresets": [
    { "name": "opt", "configurePreset": "opt", "output": { "outputOnFailure": true } }
  ]
}
```

Notes: once ccache is installed, add `"CMAKE_C_COMPILER_LAUNCHER": "ccache", "CMAKE_CXX_COMPILER_LAUNCHER": "ccache"` to `base`. Because of A1 (root `:955-956` overriding compilers) and B7 (sticky `BUILD_ALL`), **presets only take effect on a fresh build directory** — do not re-point them at the existing `build_opt`/`build_dbg`. Companion change (recommended, same commit): delete root `:955-956` so the preset/cache is the single source of compiler truth; and update `.claude/skills/ci-local` + `build-moris` to name these presets instead of the fictional `moris-ci`.

### Patch 2 — real ctest labels (`projects/CMakeLists.txt`)

Verified by two scratch probes: a directory-scope `LABELS` property is inherited by tests created in nested subdirectories, and re-setting it per loop iteration correctly labels each package without retro-effect.

```diff
 # Add source directories
 # This is done outside of the previous loop because cyclical dependencies may
 # use the variables from packages that haven't been called yet
 foreach(MORIS_SOURCE_DIR ${MORIS_SOURCE_DIRS})
+    # Label every test created under this package with the package name,
+    # so `ctest -L <PKG>` and `ctest -LE EXA` work as documented.
+    string(REGEX REPLACE "([^\\/]*)\\/" "" DIR_NAME ${MORIS_SOURCE_DIR})
+    set_property(DIRECTORY PROPERTY LABELS ${DIR_NAME})
     add_subdirectory(${MORIS_SOURCE_DIR})
 endforeach()
```

### Patch 3 — the typo-guard family

`CMakeLists.txt`:

```diff
@@ line 703
 option(MORIS_USE_SACADO "UseSacado (check paths in CMakeLists.txt)." OFF)
-if (USE_SACADO)
+if (MORIS_USE_SACADO)
@@ lines 565-584 (give each check its own cache variable)
-    check_cxx_compiler_flag("-pedantic" HAVE_PEDANTIC )
-    if (HAVE_PEDANTIC)
+    check_cxx_compiler_flag("-pedantic" HAVE_FLAG_PEDANTIC )
+    if (HAVE_FLAG_PEDANTIC)
         set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -pedantic" )
     endif()
-    check_cxx_compiler_flag("-pedantic-errors" HAVE_PEDANTIC )
-    if (HAVE_PEDANTIC)
+    check_cxx_compiler_flag("-pedantic-errors" HAVE_FLAG_PEDANTIC_ERRORS )
+    if (HAVE_FLAG_PEDANTIC_ERRORS)
         set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -pedantic-errors" )
     endif()
-    check_cxx_compiler_flag("-Wall" HAVE_PEDANTIC )
-    if (HAVE_PEDANTIC)
+    check_cxx_compiler_flag("-Wall" HAVE_FLAG_WALL )
+    if (HAVE_FLAG_WALL)
         set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Wall" )
     endif()
-    check_cxx_compiler_flag("-Werror" HAVE_PEDANTIC )
-    if (HAVE_PEDANTIC)
+    check_cxx_compiler_flag("-Werror" HAVE_FLAG_WERROR )
+    if (MORIS_USE_WERROR AND HAVE_FLAG_WERROR)
         set(MORIS_CXX_FLAGS "${MORIS_CXX_FLAGS} -Werror" )
     endif()
```

plus a new `option(MORIS_USE_WERROR "Treat warnings as errors" OFF)` near `:182` (set it ON in CI via `-DMORIS_USE_WERROR=ON` so the `ci-local` skill's claim becomes true), and in `share/cmake/dependencies/MRS_includes.cmake:46` and `LINALG_Depends.cmake:42`:

```diff
-if(USE_GPERFTOOLS) #> TEMPORARY SOLUTION
+if(MORIS_USE_GPERFTOOLS) #> TEMPORARY SOLUTION
```

(Behavior note: fixing the `USE_GPERFTOOLS` guard makes the proper `moris::gperftools` imported target active *in addition to* the raw `link_libraries` hack at root `:672-674`; the raw hack should be deleted in the same change, otherwise the system gperftools still wins over Spack's.)

---

## 5. Handoffs

- **moris-build:** owns any actual reconfigure. Specifically: (a) verify Patch 1 presets in a *fresh scratch* build dir before touching the live trees; (b) decide when/whether to retire the current `build_opt` (system g++) in favor of a preset-configured tree — this invalidates pymoris' cache fingerprints, so schedule it; (c) install `ccache` and `doxygen` (the `doxy` target will materialize on next configure once doxygen exists).
- **Implementer / main session:** apply Patches 2 and 3 (safe, no reconfigure semantics change beyond new cache vars), delete `STK_Depends.cmake` per rec 7, fix `create_shared_object.sh:130` (`MORISROOT` → `$MORISROOT`), and update `.claude/skills/ci-local/SKILL.md`, `.claude/skills/build-moris/build-rules.md`, and `moris/CLAUDE.md`'s `ctest -L` line to match reality.
- **User decision required:** (1) PYBIND — resurrect or bury (A11); (2) whether `-Werror` stays the local default or becomes CI-only (A4/Patch 3 assumes CI-only); (3) when to accept the cache invalidation of re-toolchaining `build_opt`.
- **Loud flag:** anything that changes `build_opt/build.ninja` structure (presets, LTO off, examples off) will move the `shared_object_file` edges that **pymoris scrapes** (`moris_toolchain.py`). That seam has no stability contract — whoever reconfigures must re-run the pymoris compile path end-to-end (`moris-run <deck> --compile-only`) before declaring victory.
