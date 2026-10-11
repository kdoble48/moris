# Testing

MORIS tests use [Catch2](https://github.com/catchorg/Catch2), driven by CTest. Each module
builds a test executable (`XTK-test.exe`, `INT-test.exe`, …) from its `test/` directory; with
`MORIS_USE_EXAMPLES=ON`, each example deck under `projects/EXA/` is also an integration test.

## Running tests with CTest

```bash
cd build_opt
ctest --output-on-failure            # everything
ctest -R "<regex>" --output-on-failure   # by test name
ctest -L <label> --output-on-failure     # by label, e.g. XTK, MDL
ctest -N                             # list without running
```

Run in parallel with `-j N`, sized to memory: a single example test can take several GB.
The OPT tests share an HDF5 working directory, so run them serially.

## Running a test binary directly

A module's test binary accepts Catch2 filters, which is faster than CTest when iterating on one
test:

```bash
./projects/XTK/test/bin/XTK-test.exe "[ENRICH_1]"          # one tag
./projects/XTK/test/bin/XTK-test.exe "[ENRICH_1],[GHOST]"  # several tags
./projects/XTK/test/bin/XTK-test.exe "[ENRICH_1]" -s       # show every assertion
./projects/XTK/test/bin/XTK-test.exe -t                    # list tags
```

Source the MORIS environment first; without it, the binaries fail to find the MKL or Trilinos
shared libraries.

## Writing a test

- Put unit tests in the module's `test/` directory, named `UT_<MODULE>_<Subject>.cpp`, and add
  the file to that directory's `CMakeLists.txt`.
- Give each `TEST_CASE` a tag so it can be run alone.
- Compare floating-point results with a tolerance (`Catch::Approx` or a relative-error check),
  never with `==`.
- For residual and Jacobian code, check the analytic Jacobian against a finite-difference
  approximation of the residual.
- An example deck becomes an integration test through `moris_add_example` in its directory's
  `CMakeLists.txt`; its `example_test_case.cpp` asserts on the run's results.

## Before a pull request

Build both trees and run the suites affected by the change. Compare against a run of the same
tests on the base branch: tests that fail on the base branch for environment reasons (for
example, solvers absent from the local Trilinos build) are not regressions, but new failures are.
