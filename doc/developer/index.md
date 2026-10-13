# Developer guide

For people changing the MORIS C++ core.

- {doc}`testing` — running the CTest/Catch2 suite and individual test binaries; writing tests.
- {doc}`coding-style` — naming, formatting, and comment conventions.
- {doc}`fork-workflow` — branches, pull requests, and the bug ledger of this fork.
- {doc}`documentation` — building this site and documenting code so it appears in the API
  reference.

## Repository layout

| Path | Contents |
|---|---|
| `projects/<MODULE>/src` | module sources: ALG, COM, FEM, GEN, HMR, LINALG, MIG, MOD, MRS, MTK, OPT, PRM, SOL, TOL, WRK, XTK |
| `projects/<MODULE>/test` | the module's Catch2 unit tests |
| `projects/EXA` | example decks, built as integration tests |
| `projects/PRM` | parameter-list definitions and enums for every module |
| `projects/mains` | the `moris` executable |
| `share/cmake`, `share/install` | CMake helpers and the dependency install scripts |
| `share/scripts` | deck compilation and run helpers, `moris_init.py` |
| `doc/` | this documentation site and the design records |

## Adding a physics class

A new IWG, IQI, constitutive model, or stabilization parameter is a three-layer change:

1. **Enum** — add the type to `projects/PRM/ENM/src/cl_FEM_Enums.hpp`.
2. **Class and factory** — implement `cl_FEM_<SEC>_<Name>.{hpp,cpp}` under
   `projects/FEM/INT/src/<SEC>/` and map the enum to it in `cl_FEM_<SEC>_Factory.cpp`.
3. **Test** — add a `UT_FEM_<SEC>_<Name>.cpp` Catch2 test. For residual-based classes, check
   the analytic Jacobian against finite differences.

Build and run the module tests (`INT-test`) before and after the change, so a regression is
attributable to it.
