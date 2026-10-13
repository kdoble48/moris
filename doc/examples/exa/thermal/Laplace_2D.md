# Laplace problem on a cut square with an eigenvalue solve

A 2D Laplace (steady diffusion, no source) problem on the part of the unit square to the left of
an immersed vertical plane at x = 0.52. The temperature is zero on the bottom and top edges and
y(1 − y) on the immersed plane, and the left edge is left natural. The deck compares the
solution with a 199-term series solution and, in the same solve, computes eigenpairs of the
operator. It exercises quadratic immersed elements with SPG-based enrichment and basis
agglomeration, Dirichlet conditions on an immersed boundary, and two linear-algebra back ends:
Trilinos (AztecOO GMRES plus Anasazi eigen-solvers) and PETSc/SLEPc.

```{exa-summary} projects/EXA/thermal/diffusion/Laplace_2D
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square, offset by −0.0003 in x |
| Phases | 0 (x < 0.52, the computational domain), 1 (x > 0.52, no physics) |
| Material | density 1.0, heat capacity 1.0, conductivity 0.1916229656078740 |
| Source | volumetric and surface heat loads set to 0.0 |
| Bottom and top edges | T = 0 (side ordinals 1 and 3), Nitsche-enforced |
| Immersed plane x = 0.52 | T = y(1 − y) (`Func_TempDistro`), Nitsche-enforced |
| Left edge | no condition (natural) |
| Time | static (`tTimeSteps = 0`) |

`Func_Anal_Sol` evaluates the reference solution as a sum over n = 1 … 199 of
c_n sin(nπy) cosh(nπx), with coefficients c_n computed for the domain length a = 0.52.

## Setup in MORIS

**Mesh (HMR).** 4 × 4 elements, no initial refinement. The deck fixes its own interpolation
order to 2 (`tInterpolationOrder`), giving a quadratic Lagrange mesh and quadratic B-spline
meshes on patterns 0 and 1 (`bspline_orders` = "2,2").

**Geometry and phases (GEN).** One `USER_DEFINED` level set, φ = x − 0.52, discretized on
B-spline mesh 0 with multilinear interpolation; `get_phase_index` assigns the positive side to
phase 1.

**Cut and enrichment (XTK).** Conformal decomposition, `use_SPG_based_enrichment`, enrichment of
mesh 0, `activate_basis_agglomeration`, ghost stabilization off. Enrichment fields are written
for a probe sphere (`write_enrichment_fields_probe_spheres` = "1.0,0.5,0.5,1.0"), together with
Exodus output of the integration and interpolation meshes.

**Physics (FEM).** Phase-based input. `DIFF_LIN_ISO` with `SPATIALDIFF_BULK` on phase 0;
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE` (penalty 10.0) on side
ordinals 1 and 3 and on the phase-0 side facing phase 1 (`neighbor_phases`). The same bulk and
Dirichlet IWGs are repeated for a second DOF type, `L2`, with its own constitutive model
`CMDiffusionL2`. IQIs include the analytic temperature (`PROPERTY`), `DOF` values of `TEMP` and
`L2`, `L2_ERROR_ANALYTIC`, the diffusive flux `THERMAL_ENERGY_DIFFUSIVE_FLUX` on each outer edge,
and `EIGEN_VECTOR` (index 0); `EIGEN_VALUE` is added only when `gTestIndex == 1`. MSI requests
5 eigenvectors for `gTestIndex == 0` and 1 otherwise.

**Solver (SOL).** The solver block is chosen by the global `gTestIndex`:

| | `gTestIndex = 0` (Trilinos) | `gTestIndex = 1` (PETSc) |
|---|---|---|
| Linear solve | AztecOO GMRES (`tSolverType = 3`) with `AZ_dom_decomp`/`AZ_ilut` options and preconditioner 0 (ML, `NSSA`); condition numbers computed (`dense`) | PETSc `fgmres`, `PCType` none |
| Eigen solves | two `EIGEN_SOLVER` algorithms: Block Krylov–Schur (`SM`) and Block Krylov–Schur with Amesos (`LM`), 5 eigenvalues each, Ifpack ILU preconditioner | two `SLEPC_SOLVER` power iterations (`LM`), the second with `shift_invert` and a MUMPS preconditioner |
| Nonlinear | one Newton solver over `TEMP` with algorithms 0, 1, 2 (linear solve + both eigen problems) | NLBGS over a `TEMP` sub-solver (algorithms 0 and 1: linear solve + first eigen problem) and an `L2` sub-solver |
| Matrix output | `heat` | `Jacobian_petsc` |

All Newton algorithms use one iteration and a relative residual drop of 10⁻¹².

**Output (VIS).** `Laplace.exo` on the bulk and outer side sets of both phases: nodal `TEMP`,
`HEAT`, `TEMP_A`, `L2Nodal`, `EigenVec`, `EigenVal`, `TEMPL`, and global `L2Glob` (L2 error) and
`DiffLower` (flux through the bottom edge).

## What the test checks

Two Catch2 cases call the workflow, require a return code of 0, and check `Laplace.exo.e-s.0000`
(on rank 1 for the 4-processor run):

| Case | `gTestIndex` | Built when |
|---|---|---|
| `Laplace_Anasazi` | 0 | always |
| `Laplace_Slepc` | 1 | `MORIS_HAVE_SLEPC` |

Checked values (`CHECK`, so a failure does not stop the remaining checks):

- mesh: 2 dimensions, 126 nodes and 52 elements on 1 processor; 25 nodes and 4 elements on rank 1
  of 4;
- global variable 0 (L2 error) against 6.65582 × 10⁻⁷ and global variable 1 (bottom-edge flux)
  against 0.0504293, relative difference below 10⁻⁶;
- `Laplace_Slepc` on 1 processor only: nodal field 7 at node 0 (eigenvalue) against 1.371841,
  relative difference below 10⁻⁶.

## Run it

```bash
ctest -R Laplace-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Laplace-4-procs --output-on-failure
```

## Notes

- The relative differences are computed as `(actual − reference) / reference` without an
  absolute value, so each check bounds the value from above only.
- VIS lists `IQIBulkEigenVal`, which the FEM block defines only when `gTestIndex == 1`.
- The Trilinos branch includes the `L2` IWGs in the FEM block but solves only `TEMP`.
- The `.dox` documentation file is absent from this directory.
