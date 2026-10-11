# Installing dependencies

MORIS depends on a substantial scientific stack: MPI, Trilinos (Epetra/Tpetra, Belos, Amesos,
ML, Ifpack, ROL, SEACAS/Exodus), PETSc and SLEPc, HDF5, Boost, Eigen, Armadillo, ArborX, and a
BLAS/LAPACK implementation (MKL or OpenBLAS). The supported route installs all of them with
[Spack](https://spack.io).

## The install script

The repository ships a driver script:

```bash
share/install/InstallScript.sh
```

Edit the parameters at the top of the script first (install location, compiler, developer mode,
optional packages), then run it. In **developer mode** it creates a Spack environment containing
every MORIS dependency and writes a shell setup file, `~/.bashrc_moris` (or `~/.cshrc_moris`),
that exports `MORISROOT`, the compiler wrappers, and the library paths.

Distribution-specific package prerequisites are listed in:

- `share/install/prerequisites_wsl_ubuntu_22.04`
- `share/install/prerequisites_opensuse_15_5`

If the script fails, `share/install/Install_Step_By_Step.txt` walks through the same steps by
hand.

## The shell environment

Every shell that builds or runs MORIS needs the Spack environment and the MORIS variables:

```bash
spack env activate <path-to-your-spack-env>
source ~/.bashrc_moris
```

Add the `source` line to `~/.bashrc` to load it automatically.

:::{admonition} MKL symbol-lookup errors
:class: warning

When MORIS is linked against Intel MKL, running a test binary or `moris` directly can fail with
`undefined symbol: mkl_blas_dgemm` or `libmkl_intel_lp64.so.2: cannot open shared object file`.
The MKL library directory is missing from `LD_LIBRARY_PATH`; add
`$MKL_DIR/mkl/<version>/lib` to it (and, for threading-layer symbol errors, preload
`libmkl_core.so`, `libmkl_gnu_thread.so`, and `libgomp.so.1`).
:::

## Optional packages

Some packages are optional and controlled by CMake options at configure time (see
{doc}`build`). If a package is absent from your Spack environment, switch its option off:

| Package | CMake option |
|---|---|
| PETSc | `-DMORIS_HAVE_PETSC=OFF` |
| SLEPc | `-DMORIS_HAVE_SLEPC=OFF` |
| Pardiso | `-DMORIS_USE_PARDISO=OFF` |
| MUMPS | `-DMORIS_USE_MUMPS=OFF` |
| GCMMA | `-DMORIS_HAVE_GCMMA=OFF` |
| L-BFGS | `-DMORIS_HAVE_LBFGS=OFF` |
| SNOPT (SQP) | `-DMORIS_HAVE_SNOPT=OFF` |
| Trilinos ROL | `-DMORIS_HAVE_ROL=OFF` |
