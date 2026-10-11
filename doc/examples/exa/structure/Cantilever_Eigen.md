# Cantilever eigenvalue analysis

A generalized eigenvalue analysis of a slender 2D plane-stress cantilever, with the mass matrix
as the right-hand-side matrix of the eigenproblem (`RHS_Matrix_Type = MassMat`). The example
exercises the MORIS eigen-solver path in SOL and the `EIGEN_VECTOR` IQI that writes mode shapes
to Exodus. One deck is run with three eigen-solver configurations: Trilinos Anasazi block
Davidson with an Amesos Pardiso or Umfpack inner solve, and SLEPc Krylov–Schur.

```{exa-summary} projects/EXA/structure/eigen/Cantilever_Eigen
```

## Problem

| | |
|---|---|
| Domain | 0.01 × 0.30 rectangle (width × height) |
| Material | linear elastic, isotropic, plane stress, E = 1, ν = 0.0, density 1.0 |
| Support | bottom edge (`SideSet_1_n_p0`), displacement (0, 0) by symmetric Nitsche |
| Load | traction (0, 100.0) on the left edge (`SideSet_4_n_p0`) |
| Eigenproblem | 5 eigenpairs; right-hand-side matrix `MassMat` |
| Discretization | 40 × 1200 linear elements |

## Setup in MORIS

**Mesh (HMR).** 40 × 1200 elements on the 0.01 × 0.30 domain, one Lagrange and one B-spline mesh,
both of order 1 (`tOrder = "1"`). The Lagrange mesh is also written to `Cantilever_Eigen_HMR.exo`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Back_Wall`, equal to x − 100. It
does not intersect the domain, so the whole mesh is phase 0 and the FEM sets are the uncut
`HMR_dummy_n_p0` block and the `SideSet_<k>_n_p0` boundary sets
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
off, Exodus output of the integration mesh.

**Physics (FEM).** One `STRUC_LIN_ISO` model with `Model_Type::PLANE_STRESS`. IWGs:
`STRUC_LINEAR_BULK` on the bulk, `STRUC_LINEAR_NEUMANN` on the left edge,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` on the bottom edge (`DIRICHLET_NITSCHE`, penalty 100.0),
and `TIME_CONTINUITY_DOF` on the bulk with current and previous weights equal to the density and a
zero initial condition. IQIs: `VOLUME`, and two `EIGEN_VECTOR` IQIs (`IQIBulkDISPX`,
`IQIBulkDISPY`, `vectorial_field_index` 0 and 1, `function_parameters = "1"`).

**Solver (SOL).** The deck reads the global `gPrecSolver` (set by the test case) and builds one of
two parameter sets:

| | Trilinos (`Amesos_Pardiso`, `Amesos_Umfpack`) | PETSc/SLEPc (`Slepc`) |
|---|---|---|
| Eigen algorithm | `EIGEN_SOLVER`, `EIGALG_BLOCK_DAVIDSON`, `Which = SM` | `SLEPC_SOLVER`, `krylovschur`, `Which = LM`, `STType = shift_invert` |
| Eigenvalues | `Num_Eig_Vals = 5`, `Block_Size = 5`, `Num_Blocks = 2`, `MaxSubSpaceDims = 75`, `MaxRestarts = 20`, tolerance 10⁻⁵ (relative) | `Num_Eig_Vals = 5`, `is_symmetric = true` |
| Inner solve | `IFPACK` preconditioner with `ifpack_prec_type = Amesos`, `amesos: solver type = gPrecSolver` | PETSc `gmres` with `PCType = mumps` |
| Other | `SOL_save_operator_to_matlab = MassMat` | `SOL_TPL_Type = Petsc` |

Both sets use `RHS_Matrix_Type = MassMat`, one Newton iteration (`NLA_rel_res_norm_drop = 1e-10`),
and one time step over a time frame of 1. MSI requests `number_eigen_vectors = 5`.

**Output (VIS).** `STANDARD` mesh on `HMR_dummy_n_p0`, nodal fields `UX, UY` from the two
eigenvector IQIs, written to `Cantilever_Eigen_<gPrecSolver>.exo`.

## What the test checks

The ctest runs one Catch2 binary invocation with three test cases; the SLEPc case is compiled only
when `MORIS_HAVE_SLEPC` is defined.

| Test case | `gPrecSolver` | Output file | Reference node | Reference (UX, UY) |
|---|---|---|---|---|
| `Cantilever_Eigen_Pardiso` (index 0) | `Amesos_Pardiso` | `Cantilever_Eigen_Amesos_Pardiso.exo` | 14 at (0.00175, 0.0) | (−1.42942084733364e−06, 6.46851986039073e−08) |
| `Cantilever_Eigen_Umfpack` (index 1) | `Amesos_Umfpack` | `Cantilever_Eigen_Amesos_Umfpack.exo` | 43 at (0.00525, 0.00025) | (9.93543881481338e−04, −5.32691066062199e−04) |
| `Cantilever_Eigen_Slepc` (index 2) | `Slepc` | `Cantilever_Eigen_Slepc.exo` | 43 at (0.00525, 0.00025) | same as index 1 |

Each case requires a return code of 0, then checks 2 dimensions, 49241 nodes and 48000 elements,
the reference node coordinates (relative tolerance 10⁻⁵), and the nodal eigenvector displacement
(relative tolerance 10⁻⁵). Before the comparison the computed vector is negated if its x component
has the opposite sign to the reference, since an eigenvector is defined only up to sign.
Eigenvalues are not compared.

## Run it

```bash
ctest -R Cantilever_Eigen-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The Pardiso case uses the `Amesos_Pardiso` backend, which the Trilinos build must provide (see
the Amesos note in {doc}`../../../user-guide/modules`); the SLEPc case is registered only with
`MORIS_HAVE_SLEPC`.

## Notes

- The deck declares `extern uint gOrder` but sets the interpolation order from its own constant
  `tOrder = "1"`.
- The `NumFreeDofs = 1000` key in the Trilinos list carries the comment "number of free dofs =
  2*node_x*node_y"; the value is not tied to this mesh.
