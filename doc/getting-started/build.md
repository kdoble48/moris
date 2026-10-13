# Building MORIS

MORIS uses CMake with the Ninja generator and out-of-source build trees. The convention is two
trees next to the sources: `build_opt/` (optimized) and `build_dbg/` (debug, with assertions).

## With presets

`CMakePresets.json` defines the standard configurations:

| Preset | Tree | Purpose |
|---|---|---|
| `opt` | `build_opt` | Release build with the Python bindings (`moris_py`) |
| `dbg` | `build_dbg` | Debug build with `MORIS_HAVE_DEBUG` and `MORIS_HAVE_SYMBOLIC` |
| `ci-opt`, `ci-dbg` | as above | CI variants: as `opt`/`dbg`, with the example decks built |

```bash
cmake --preset opt
cmake --build build_opt
```

The `base` preset that the others inherit from points at a specific Spack environment. If your
environment lives elsewhere, override the compiler and prefix paths in a personal
`CMakeUserPresets.json` rather than editing the shared file.

## By hand

An equivalent optimized configuration without presets:

```bash
cd $MORISROOT
cmake -S . -B build_opt -G Ninja \
      -DBUILD_ALL=ON \
      -DMORIS_USE_EXAMPLES=OFF \
      -DBUILD_PYBINDINGS=ON \
      -DMORIS_USE_MUMPS=OFF \
      -DMORIS_HAVE_SNOPT=OFF
cmake --build build_opt
```

and a debug configuration:

```bash
cmake -S . -B build_dbg -G Ninja \
      -DBUILD_ALL=ON \
      -DMORIS_USE_EXAMPLES=OFF \
      -DMORIS_HAVE_DEBUG=ON \
      -DMORIS_HAVE_SYMBOLIC=ON \
      -DMORIS_USE_MUMPS=OFF \
      -DMORIS_HAVE_SNOPT=OFF
cmake --build build_dbg
```

On an Intel machine, `-DMORIS_USE_MKL=ON -DMORIS_USE_OPENBLAS=OFF` selects MKL instead of
OpenBLAS.

## Commonly used options

| Option | Default | Effect |
|---|---|---|
| `BUILD_ALL` | `OFF` | Build every module and the `moris` executable |
| `MORIS_USE_EXAMPLES` | `ON` | Build the `projects/EXA` example decks (they double as integration tests) |
| `MORIS_HAVE_DEBUG` | `OFF` | Enable `MORIS_ASSERT` checks and debug output |
| `MORIS_HAVE_SYMBOLIC` | `OFF` | Symbolic-factorization checks in the linear algebra layer |
| `BUILD_PYBINDINGS` | `OFF` | Build the `moris_py` Python module |
| `MORIS_HAVE_ROL` | `ON` | Trilinos ROL optimizer (`Algorithm_ROL`) |

## Keep examples off during development

The ~64 example decks under `projects/EXA/` each link a separate executable. That linking tail
dominates incremental build time: touching a widely included header relinks every example even
when the libraries rebuilt in seconds. For library work and unit tests, configure with
`-DMORIS_USE_EXAMPLES=OFF`, or build only the targets you need:

```bash
ninja -C build_opt INT-lib INT-test projects/mains/moris
```

A named target builds itself and its dependencies, never the example executables.

## Memory

A single MORIS translation unit can take 1–3 GB in the compiler, and the optimized tree uses
link-time optimization. Size the job count (`ninja -j N`) to the memory you have, not the cores:
roughly one job per 3 GB of free memory.

## Verifying the build

```bash
./build_opt/projects/mains/moris --version
cd build_opt && ctest --output-on-failure
```

See {doc}`../developer/testing` for running subsets of the test suite.
