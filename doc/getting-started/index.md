# Getting started

Three steps take you from a clean machine to a finished MORIS run:

1. {doc}`install` — install the third-party libraries (Trilinos, PETSc, MPI, HDF5, …) with the
   Spack-based install script.
2. {doc}`build` — configure and build MORIS with CMake and Ninja, in an optimized or a debug tree.
3. {doc}`running` — compile an input deck to a shared object and run it with the `moris`
   executable.

MORIS is developed and tested on Linux (openSUSE, Ubuntu, and Ubuntu under WSL2). The
dependency stack is large; expect the first install to take a few hours, most of it compiling
Trilinos.
