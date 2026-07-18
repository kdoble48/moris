# CMake Audit Report — MORIS

Date: 2025-09-28

## Executive Summary

- Score: 52/100
- Overall, the build is functional in its original environment, but relies heavily on global state, legacy CMake constructs, and unconditional dependencies that break fresh, minimal configures. Packaging and per-package exports are in place and a positive, but target‑scoped usage requirements and cross‑platform concerns need significant work.

Key strengths
- Per‑package `Config.cmake` and `Targets.cmake` export flow exists; consumers get `moris::COMPONENT` targets (share/cmake/morisConfig.cmake.in:63).
- Imported targets are created in many custom Find modules (e.g. `OPENBLAS::openblas`).
- Installation exports are present across subpackages.

Key issues
- Unconditional inclusion of MRS components forces Boost/HDF5/MPI on every configure, causing failures in minimal environments (CMakeLists.txt:1008 via share/cmake/dependencies/MRS_includes.cmake:12).
- Extensive use of global commands and variables (`include_directories`, `add_definitions`, `link_libraries`, `CMAKE_CXX_FLAGS`, forced compilers) instead of target‑scoped `target_*`.
- Tooling and presets missing; no `CMakePresets.json`.
- CMake minimum is 3.17; aim 3.20+ and set `cmake_policy(VERSION)`.
- Cross‑platform gaps (hard‑coded GCC/ELF flags, mpic++ defaults) and forced `CMAKE_BUILD_TYPE`.

## Inventory & Structure

- Root: `CMakeLists.txt` (CMake 3.17; `project(MORIS 1.0.0)`) — CMakeLists.txt:16, CMakeLists.txt:24
- CMake files:
  - CMakeLists: 182 files
  - Modules (`*.cmake`): 94 files
- Notable trees
  - `projects/*` — per‑package libs/tests/examples
  - `share/cmake/dependencies/*` — inter‑package dependency wiring
  - `share/cmake/find_modules/*` — custom Find modules (define imported targets)
  - `share/cmake/find_wrappers/*` — wrappers; new variants add imported targets, older ones set globals
  - `share/cmake/utilities/*` — macros, helpers

Root policies
- `cmake_minimum_required(VERSION 3.17)` (CMakeLists.txt:16)
- No `cmake_policy(VERSION ...)` nor explicit `cmake_policy(SET ...)` detected.

Root global state (selected)
- Forced compilers to mpicc/mpic++ if not set (CMakeLists.txt:300–325)
- Global `link_libraries()` adding system libs (CMakeLists.txt:385–453, 438–453)
- Global definitions and includes (CMakeLists.txt:980–983)
- Forced `CMAKE_BUILD_TYPE` (Debug/Release) (CMakeLists.txt:736–770)
- Manual global flags for warnings/opts/LTO (CMakeLists.txt:552–635)
- Unconditional MRS include (CMakeLists.txt:1008 → share/cmake/dependencies/MRS_includes.cmake:12)

## Target Inventory (high‑level)

Library targets (representative; “-lib” outputs use `OUTPUT_NAME` to drop suffix):
- `${LINALG}-lib` (INTERFACE) — projects/LINALG/src/CMakeLists.txt:205
  - Publishes LINALG API; TPLs: `${ARMADILLO_EIGEN}`, `${ACML_LAPACK_MKL_OPENBLAS}`, `superlu`, `arpack` (share/cmake/dependencies/LINALG_Depends.cmake:19–36)
  - Moris deps: `${COM}`, `${ALG}`
- `${COM}-lib` — projects/COM/src/CMakeLists.txt:37
  - Moris deps: `${LINALG}-lib`
  - TPLs: `mpi` via wrapper (share/cmake/dependencies/COM_Depends.cmake:20–31)
- `${ALG}-lib` — projects/ALG/src/CMakeLists.txt:41
  - Moris deps: `${LINALG}-lib`
- `${TOL}-lib` — projects/TOL/src/CMakeLists.txt:32–35
  - Moris deps: `${LINALG}-lib`, `${ALG}-lib`
- `${HMR}-lib` — projects/HMR/src/CMakeLists.txt:125–140
  - Moris deps: `${LINALG}-lib`, `${CNT}-lib`, `${IOS}-lib`, `${ALG}-lib`, `${MTK}-lib`, `${MAP}-lib`
- `${MSI}-lib` — projects/FEM/MSI/src/CMakeLists.txt:61
  - Moris deps: `${DLA}-lib`, `${PRM}-lib`, `${ENM}-lib`
- `${DLA}-lib` — projects/SOL/DLA/src/CMakeLists.txt:94
  - Moris deps: `${LINALG}-lib`, `${COM}-lib`, `${MTK}-lib`, `${SOL_CORE}-lib`
- `${NLA}-lib` — projects/SOL/NLA/src/CMakeLists.txt:63
  - Moris deps: `${LINALG}-lib`, `${DLA}-lib`, `${COM}-lib`, `${MSI}-lib`, `${SOL_CORE}-lib`, `${ENM}-lib`
- `${VIS}-lib` — projects/FEM/VIS/src/CMakeLists.txt:63
  - Moris deps: `${LINALG}-lib`, `${COM}-lib`, `${CNT}-lib`, `${HMR}-lib`, `${MTK}-lib`, `${ENM}-lib`
- `${INT}-lib` — projects/FEM/INT/src/CMakeLists.txt:566
  - Moris deps: `${LINALG}-lib`, `${MTK}-lib`, `${MSI}-lib`, `${VIS}-lib`, `${PRM}-lib`, `${ENM}-lib`
- `${WRK}-lib` — projects/WRK/src/CMakeLists.txt:77
  - Moris deps: `${LINALG}-lib`, `${INT}-lib`, `${XTK}-lib`, `${HMR}-lib`, `${GEN}-lib`, `${MDL}-lib`, `${OPT}-lib`, `${MIG}-lib`
- `${GEN}-lib` — projects/GEN/src/CMakeLists.txt:171
  - Moris deps: `${COM}-lib`, `${LINALG}-lib`, `${MAP}-lib`, `${SDF}-lib`, `${ENM}-lib`
- `${SDF}-lib` — projects/GEN/SDF/src/CMakeLists.txt:70
  - Moris deps: `${LINALG}-lib`, `${COM}-lib`, `${CNT}-lib`, `${MTK}-lib`
- `${SOL_CORE}-lib` — projects/SOL/SOL_CORE/src/CMakeLists.txt:69
  - Moris deps: `${DLA}-lib`, `${ENM}-lib`
- `${MAP}-lib` — projects/MTK/MAP/src/CMakeLists.txt:44
- `${MTK}-lib` — projects/MTK/src/CMakeLists.txt:343 (deps: `${COM}`, `${LINALG}`, `${HMR}`, `${TOL}`, `${ENM}`; projects/MTK/src/CMakeLists.txt:332–341)
- `${XTK}-lib`, `${OPT}-lib`, `${MOD}-lib`, `${MIG}-lib`, `${PRM}-lib`, `${ENM}-lib`, `${IOS}-lib`, `${CNT}-lib`, `${ASR}-lib`, `${COR}-lib` (various)

Executables (selected)
- `moris` main — projects/mains/CMakeLists.txt:59 (links many `${..}-lib` and raw MORIS_* variables)
- Dozens of example/test executables under `projects/EXA/**` and per‑package `test/` trees.

Target usage requirements
- Many targets set some `target_include_directories()` and `target_link_libraries()` with `PUBLIC` (good), but also inject global include/def lists from wrappers: e.g. projects/ALG/src/CMakeLists.txt:57–67, projects/COM/src/CMakeLists.txt:45–57.

## Modern CMake Compliance

Target‑scoped vs global
- Global: `add_definitions(${MORIS_DEFINITIONS})`, `include_directories(${MORIS_INCDIRS})`, `link_libraries(...)` (CMakeLists.txt:980–983, 385–453), numerous `include_directories()` in tests/examples.
- Target‑scoped: used within many package `CMakeLists.txt` for local includes and link deps (good), but mixed with global variables injected by wrappers.

C++ standard
- Only `projects/GUI/src/CMakeLists.txt` sets `CMAKE_CXX_STANDARD 17` (projects/GUI/src/CMakeLists.txt:14–15). Root uses custom `HAVE_STD_CPP17` detection to append `-std=c++17` to global flags (CMakeLists.txt:610–620). Prefer `target_compile_features(<tgt> PUBLIC cxx_std_17)` for transitivity.

Usage requirements
- Many links are `PUBLIC`, but TPLs are not consistently propagated via imported targets; wrappers often add plain include/defs variables.

Generator expressions
- Minimal usage. Recommend `$<BUILD_INTERFACE:...>` and `$<INSTALL_INTERFACE:...>` are used in parts (good), but could be applied consistently and replace many globals.

## Options, Presets, and Tooling

- Options: extensive `MORIS_*` options (good prefixing). Some mutually exclusive sets are enforced (e.g., MKL vs LAPACK vs OPENBLAS).
- Cache pollution: multiple global `set()` and `include_directories()` add to process‑wide state.
- Presets: No `CMakePresets.json`.
- Tooling: No ccache/sccache; LTO done via flags, not `INTERPROCEDURAL_OPTIMIZATION`; no unity build toggle; reproducible build flags absent.

Suggested Presets (skeleton)
```json
{
  "version": 5,
  "configurePresets": [
    {
      "name": "default",
      "generator": "Ninja",
      "binaryDir": "build",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "RelWithDebInfo",
        "MORIS_USE_TESTS": "OFF"
      }
    },
    { "name": "asan", "inherits": "default", "binaryDir": "build-asan",
      "cacheVariables": {
        "CMAKE_BUILD_TYPE": "RelWithDebInfo",
        "CMAKE_CXX_FLAGS": "-fsanitize=address,undefined -fno-omit-frame-pointer",
        "CMAKE_LINKER_FLAGS": "-fsanitize=address,undefined"
      }
    },
    { "name": "release", "inherits": "default", "cacheVariables": {"CMAKE_BUILD_TYPE": "Release"} }
  ],
  "buildPresets": [ {"name":"default","configurePreset":"default","jobs":0} ],
  "testPresets": [ {"name":"default","configurePreset":"default","output": {"outputOnFailure": true}} ]
}
```

## Dependency Management

Wrappers present (selection from `share/cmake/find_wrappers`):
- `acml(_new).cmake`, `openblas(_new).cmake`, `mkl(_new).cmake`, `lapack(_new).cmake`
- `mpi(_new).cmake`, `petsc_new.cmake`, `slepc_new.cmake`, `trilinos(_new|_new_install).cmake`
- `eigen(_new).cmake`, `arpack(_new).cmake`, `superlu(_new).cmake`, `suitesparse(_new).cmake`, `boost(_new).cmake`, `hdf5_new.cmake`, `gcmma(_new).cmake`, `snopt(_new).cmake`, `viennacl(_new).cmake`, `gperftools_new.cmake`, `arborx_new.cmake`

Custom Find modules (selection from `share/cmake/find_modules`):
- Define imported targets (e.g., `OPENBLAS::openblas`, `MKL::mkl`, `LAPACK::lapack`). Good modern pattern.

Concerns
- Mixed old/new wrappers: older wrappers still call `add_definitions()` and `include_directories()` globally (e.g., share/cmake/find_wrappers/openblas.cmake:30–33).
- Some wrapper `find_package()` use REQUIRED (good), others optional; but package missing → later `target_link_libraries(${MORIS}::TPL)` may reference non‑existent imported targets.
- Unconditional MRS includes drag `boost_new.cmake` (REQUIRED) into all configs (see “Run & Validate”).

Recommendations
- Prefer imported targets from Find modules/wrappers consistently, and remove global include/defs variables.
- Add `REQUIRED` with versions for core TPLs (Eigen/OpenBLAS/etc) at the app entry point; guard optional ones.
- Replace ad‑hoc `FindEigen.cmake` with upstream `Eigen3Config` usage where possible: `find_package(Eigen3 3.3 CONFIG REQUIRED)` and use `Eigen3::Eigen`.

## Installation, Packaging, Export

What’s present
- Per‑package `configure_package_config_file()` and `write_basic_package_version_file()` usage (good).
- `install(TARGETS ... EXPORT ...)` + `install(EXPORT ...)` per package (good).
- Top‑level general interface target `${MORIS}_general` exporting compile definitions (CMakeLists.txt:1185–1200).

Gaps
- No `GNUInstallDirs`; custom dirs used throughout.
- No `NAMESPACE` on `install(EXPORT ...)` (consumers get raw target names, but `morisConfig.cmake.in` compensates by creating `moris::COMPONENT` wrapper INTERFACE targets at consume time).
- `SOVERSION`/`VERSION` not set for shared libs (relevant if `MORIS_HAVE_SHARED=ON`).
- Relocatability generally okay via `PACKAGE_PREFIX_DIR`, but hand‑rolled includes in generated configs string‑concat wrapper includes; consider `find_dependency()`.

Snippets
- Prefer
```cmake
include(GNUInstallDirs)
install(TARGETS mylib EXPORT morisTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(EXPORT morisTargets NAMESPACE moris:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/moris)
```

## Testing & CI

- `enable_testing()` present (CMakeLists.txt:1166). Tests are added manually via `add_test()` in per‑package trees.
- Catch2 vendored header exists (`include/catch.hpp`) but no `catch_discover_tests()` usage.

Recommendations
- Add test presets: `ctest --test-dir <build> --output-on-failure`.
- If moving to Catch2 v3, use `catch_discover_tests(<exe>)`.
- CI: matrix over Linux(gcc/clang), macOS, Windows; drive via CMake Presets (see “Validation Commands”).

## Warnings, Sanitizers, Security

- Warnings: `-Wall -Werror -pedantic(-errors)` added globally (CMakeLists.txt:565–583). Suggest `-Werror` only in CI via option.
- Sanitizers: not present; add ASan/UBSan presets (see Presets above).
- `POSITION_INDEPENDENT_CODE`: not set on targets; only global `-fPIC`. Prefer `set_target_properties(<tgt> PROPERTIES POSITION_INDEPENDENT_CODE ON)` when needed.
- RPATH: not managed explicitly; okay to rely on defaults, but consider install RPATH policy if shipping shared libs.
- Windows DLL exports: no strategy in tree.

## Cross‑platform & Toolchains

Issues
- Hard‑coded compilers (mpicc/mpic++) and tool flags (`-rdynamic`, `-m64`, `-xCORE-AVX2`, etc.) applied globally (CMakeLists.txt:300–325, 620–666).
- No MSVC handling (`/W4`, `/EHsc`, `/permissive-`, `/bigobj`), no Windows library guards for `dl`, `rt`, `ssl`, `crypto`.

Recommendations
- Guard flags by compiler and platform; prefer `target_compile_options()` with generator expressions.
- Avoid setting `CMAKE_CXX_FLAGS` globally; use an INTERFACE “compile options” target linked by all libs.
- Do not overwrite `CMAKE_BUILD_TYPE`; let presets drive it.

Example
```cmake
if(MSVC)
  target_compile_options(moris_base INTERFACE /W4 /EHsc /permissive-)
else()
  target_compile_options(moris_base INTERFACE -Wall -Wextra -Wpedantic)
endif()
```

## Findings (by severity)

Critical
- Unconditional MRS inclusion forces Boost/HDF5/MPI and breaks configure in clean envs.
  - CMakeLists.txt:1008 → share/cmake/dependencies/MRS_includes.cmake:12
  - Error (see log): `Could not find a package configuration file provided by "Boost"`.
- Forced compilers to mpicc/mpic++ and forced `CMAKE_BUILD_TYPE`.
  - CMakeLists.txt:300–325, 736–770

High
- Global `link_libraries`, `include_directories`, `add_definitions` at root.
  - CMakeLists.txt:385–453, 980–983
- Global flags for optimization, LTO, warning levels (no per‑target control).
  - CMakeLists.txt:552–666
- Mixed old/new wrapper styles; some wrappers leak includes/defs globally.
  - Example: share/cmake/find_wrappers/openblas.cmake (global) vs openblas_new.cmake (imported target)
- No `cmake_policy(VERSION)` and minimum < 3.20.

Medium
- No `CMakePresets.json` for reproducible developer flow.
- No `GNUInstallDirs`; install dirs hard‑coded.
- No sanitizer presets; `-Werror` globally (prefer CI‑only).
- Incomplete Windows/MSVC support and guards.

Low
- Sparse use of generator expressions beyond install/build interface includes.
- MAIN executable links raw MORIS_* vars instead of imported targets.
  - projects/mains/CMakeLists.txt:31–68

## Quick‑Fix Patchset (top 5)

1) Raise minimum CMake and set policies (root CMakeLists)
```diff
@@
-cmake_minimum_required(VERSION 3.17)
+cmake_minimum_required(VERSION 3.20)
+cmake_policy(VERSION 3.20)
```

2) Make MRS includes optional behind a toggle (default OFF)
```diff
@@
-include(${MORIS_DEPENDS_DIR}/MRS_includes.cmake)
+option(MORIS_ENABLE_MRS "Include core MRS components by default" OFF)
+if(MORIS_ENABLE_MRS)
+  include(${MORIS_DEPENDS_DIR}/MRS_includes.cmake)
+endif()
```

3) Stop forcing compilers/build type; prefer user/preset control
```diff
@@
-if (NOT MORIS_C_COMPILER )
-    set(MORIS_C_COMPILER "mpicc" )
-endif()
-if (NOT MORIS_CXX_COMPILER )
-    set(MORIS_CXX_COMPILER "mpic++")
-endif()
-set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type Release" FORCE)
+# Respect externally provided compilers; do not force MPI wrappers.
+# Toolchain or presets should set compilers as needed.
+# set(MORIS_C_COMPILER "mpicc")
+# set(MORIS_CXX_COMPILER "mpic++")
+# Do not force CMAKE_BUILD_TYPE; use presets instead.
```

4) Replace root‑level global include/defs with an interface target
```cmake
add_library(moris_base INTERFACE)
target_compile_definitions(moris_base INTERFACE ${MORIS_DEFINITIONS})
target_include_directories(moris_base INTERFACE ${MORIS_INCDIRS})
# Link ‘moris_base’ to all package libraries
# target_link_libraries(${ALG}-lib PUBLIC moris_base) # repeat for packages
```

5) Adopt target‑scoped C++ standard and IPO knobs
```cmake
target_compile_features(${ALG}-lib PUBLIC cxx_std_17)
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_ok)
if(ipo_ok)
  set_property(TARGET ${ALG}-lib PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
```

Additional suggested files
- CMakePresets.json — see skeleton above.
- cmake/Conventions.cmake — centralize warnings and sanitizer flags (target‑scoped) and include it in each package CMakeLists.

## Refactor Plan (prioritized)

1. Add presets and raise CMake min + policies
2. Gate MRS includes behind `MORIS_ENABLE_MRS` (default OFF)
3. Remove root `include_directories`, `add_definitions`, `link_libraries`
4. Introduce `moris_base` INTERFACE target and link it in all libs
5. Standardize wrapper usage to “_new” imported‑target variants only
6. Switch TPL usages to imported targets; drop global include/defs variables
7. Replace global flags with `target_compile_options` per lib
8. Introduce `GNUInstallDirs` and `NAMESPACE` on install exports
9. Add sanitizer and CI presets; move `-Werror` to CI
10. Add Windows/MSVC guards and equivalent flags

## Validation Commands

Linux/macOS (Ninja)
```sh
cmake -S . -B build -G Ninja -DMORIS_ENABLE_MRS=OFF -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build -j
ctest --test-dir build --output-on-failure
```

With presets (after adding CMakePresets.json)
```sh
cmake --preset default
cmake --build --preset default
ctest --preset default
```

Windows (Visual Studio generator)
```bat
cmake -S . -B build -G "Visual Studio 17 2022" -A x64 -DMORIS_ENABLE_MRS=OFF
cmake --build build --config RelWithDebInfo
ctest --test-dir build --output-on-failure -C RelWithDebInfo
```

## Run & Validate (this environment)

Commands attempted
```sh
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DMORIS_HAVE_PARALLEL=OFF -DMORIS_USE_TESTS=OFF -DBUILD_ALL=OFF -DBUILD_MAIN=OFF \
  -DMORIS_USE_EXAMPLES=OFF -DMORIS_USE_ARMA=OFF -DMORIS_USE_EIGEN=ON -DMORIS_USE_OPENBLAS=ON
```

Result: configure failed due to unconditional Boost requirement from MRS includes
```
-- CMAKE_MODULE_PATH: .../share/cmake/find_modules;.../share/cmake/find_wrappers
-- MORIS_CXX_COMPILER is set to mpic++
-- MORIS_LIB_Ssl = MORIS_LIB_Ssl-NOTFOUND
-- MORIS_Lib_Crypto = MORIS_Lib_Crypto-NOTFOUND
-- MORIS_LIB_Z = MORIS_LIB_Z-NOTFOUND
-- MORIS recognized the C++ flags ... -std=c++17 -m64 -fextended-identifiers -rdynamic -fPIC.
CMake Error at share/cmake/find_wrappers/boost_new.cmake:18 (find_package):
  Could not find a package configuration file provided by "Boost" (requested
  version 1.54.0) with any of the following names:
    BoostConfig.cmake
    boost-config.cmake
  Add the installation prefix of "Boost" to CMAKE_PREFIX_PATH or set "Boost_DIR"...
Call Stack (most recent call first):
  projects/MRS/CHR/src/CMakeLists.txt:27 (include)
```

Notes captured
- Missing system libs (ssl, crypto, z) are linked globally; should be target‑scoped and optional.
- For a clean configure/build in general CI, the MRS include should be optional or the TPL set provided via presets/toolchain.

---

If you want, I can draft the first set of patches (min CMake + policy, optional MRS, presets skeleton, and a `moris_base` interface target) and/or wire a basic GitHub Actions CI using presets.

