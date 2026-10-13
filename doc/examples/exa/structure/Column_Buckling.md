# Column linear buckling analysis

A linear buckling analysis (LBA) of a column clamped at its base and compressed by a traction on
its top face. A static solve provides the prestress, then a generalized eigenproblem with the
geometric stiffness matrix as right-hand-side matrix (`RHS_Matrix_Type = GeomStiffMat`) yields the
five lowest buckling eigenvalues. The example tests the `STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS` IWG
with a linear or a PK2 (Saint Venant–Kirchhoff) stress measure, the `EIGEN_VALUE` IQI, and the
SLEPc eigen solver, on a body-fitted column in 2D and 3D and on an immersed column cut from a
larger background mesh.

```{exa-summary} projects/EXA/structure/eigen/Column_Buckling
```

## Problem

| | |
|---|---|
| Column | 1.0 × 10.0 (2D) or 1.0 × 10.0 × 1.0 (3D); volume 10.0 |
| Material | linear elastic, isotropic, E = 100, ν = 0.3 |
| Support | bottom face y = 0: displacement 0 by symmetric Nitsche |
| Load | top face y = 10: traction (0, −1.0) or (0, −1.0, 0) |
| Analysis | static prestress solve, then 5 eigenpairs of the buckling problem |

The deck reads four globals set by the test case: `tOrder` (interpolation order), `tDim` (2 or
3), `tStressType` (`Linear` or `PK2`) and `tOutputFileName`.

## Setup in MORIS

**Mesh (HMR).** Body-fitted deck: 5 × 50 (2D) or 5 × 50 × 5 (3D) elements on the column domain.
Immersed deck: 10 × 55 (× 10) elements on a 2 × 11 (× 2) background domain with
`domain_offset` −0.5 in every direction, so the background covers [−0.5, 1.5] × [−0.5, 10.5]
(× [−0.5, 1.5]). One Lagrange and one B-spline mesh, both of order `tOrder`.

**Geometry and phases (GEN).** Body-fitted deck: no geometry; FEM uses `HMR_dummy_n_p0`,
`SideSet_1_n_p0` (bottom) and `SideSet_3_n_p0` (top). Immersed deck: four `LINE` level sets in 2D
(six `PLANE` level sets in 3D) at y = 0, y = 10, x = 0, x = 1 (and z = 0, z = 1). The phase
function `get_phase_index_2d` / `get_phase_index_3d` assigns four phases: 3 outside any lateral
face, 1 below y = 0, 2 above y = 10, and 0 for the column itself.

**Cut and enrichment (XTK).** Conformal decomposition and enrichment of mesh 0 in both decks;
ghost stabilization off in the body-fitted deck and on in the immersed deck.

**Physics (FEM).** Two constitutive models on the same properties: `STRUC_LIN_ISO` and
`STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF`. IWGs:

- `STRUC_LINEAR_BULK` with the linear model;
- `STRUC_NON_LINEAR_GEOMETRIC_STIFFNESS` with the SVK model when `tStressType == "PK2"`, else with
  the linear model;
- `STRUC_LINEAR_NEUMANN` with the traction on the top face;
- `STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` on the bottom face (`DIRICHLET_NITSCHE`, penalty 100.0).

The immersed deck declares FEM phases (`PhaseBulk` = 0, `PhaseBottom` = 1, `PhaseTop` = 2,
`PhaseSides` = 3) and places the boundary IWGs on `SIDESET` elements of `PhaseBulk` with
`neighbor_phases` `PhaseTop` (traction) and `PhaseBottom` (support) instead of naming mesh sets.
It adds a `GHOST_NORMAL_FIELD` IWG with a `GHOST_DISPL` parameter (0.001) and `ghost_order` equal
to the interpolation order. IQIs in both decks: `VOLUME`, five `EIGEN_VALUE` IQIs
(`vectorial_field_index` 0–4), and per eigenvector `EIGEN_VECTOR` IQIs for each displacement
component (`function_parameters` = eigenvector index 0–4).

**Solver (SOL).** `tUsePetsc = true` selects the PETSc parameter set; a Trilinos set
(`EIGALG_BLOCK_DAVIDSON` with an Amesos Pardiso preconditioner) is also defined but not used. The
nonlinear solver runs two algorithms in sequence (`NLA_Nonlinear_solver_algorithms = "0,1"`):

| Step | Nonlinear algorithm | Linear solver |
|---|---|---|
| 0 — prestress | `NLA_max_iter = 2`, `NLA_rel_res_norm_drop = 1e-10`, combined residual/Jacobian assembly | PETSc `gmres` with `PCType = mumps` |
| 1 — buckling | `NLA_max_iter = 1`, `NLA_is_eigen_problem = true` | `SLEPC_SOLVER`, `krylovschur`, `Which = LM`, `Num_Eig_Vals = 5`, `STType = shift_invert`, `is_symmetric = false`, `RHS_Matrix_Type = GeomStiffMat` |

`SOL_TPL_Type = Petsc`; MSI requests `number_eigen_vectors = 5`.

**Output (VIS).** `STANDARD` mesh on the bulk set(s): global fields `VOLUME` and `EVAL1`–`EVAL5`,
nodal fields `EVEC<i>X`, `EVEC<i>Y` (and `EVEC<i>Z` in 3D), written to `tOutputFileName`.

## Cases

All test cases live inside `#ifdef MORIS_HAVE_SLEPC`. Each ctest runs every case at its
processor count; the case index is the first value on 1 processor and the second on 4.

| Test case | Deck | Dim | Order | Stress | Index (1 / 4 procs) | Reference eigenvalues 1–5 |
|---|---|---|---|---|---|---|
| `Column_Buckling_2D` | body-fitted | 2 | 1 | Linear | 0 / 1 | 0.208264, 1.76392, 4.39768, 7.51133, 10.6635 |
| `Column_Buckling_2D` | body-fitted | 2 | 1 | PK2 | 2 / 3 | 0.209892, 1.77737, 4.42993, 7.56415, 10.7357 |
| `Column_Buckling_2D` | body-fitted | 2 | 2 | Linear | 4 / 5 | 0.204191, 1.72616, 4.28966, 7.29788, 10.3161 |
| `Column_Buckling_3D` | body-fitted | 3 | 1 | Linear | 6 / 7 | 0.209678, 0.209678, 1.77607, 1.77607, 4.42856 |
| `Column_Buckling_Immersed_2D` | immersed | 2 | 1 | Linear | 8 / 9 | 0.208056, 1.76131, 4.38765, 7.48695, 10.618 |
| `Column_Buckling_Immersed_3D` | immersed | 3 | 1 | Linear | 10 / 11 | 0.209338, 0.209338, 1.77238, 1.77238, 4.41592 |

`Column_Buckling_2D` runs the deck three times in one test case, changing the globals between
runs. Output files are named `Column_Buckling[_Immersed]_<stress>_PolyOrd_<order>_<dim>D.exo`.

## What the test checks

After each run the test requires a return code of 0; rank 0 then opens the output file and checks:

- dimension, node count and element count against per-index references (for example 306 nodes
  and 250 elements for index 0, 104 and 75 for index 1, 10184 and 16168 for index 10);
- the global `VOLUME` against 10.0, relative tolerance 10⁻⁵;
- each of the five global eigenvalue fields divided by the volume against the table above,
  relative tolerance 10⁻⁵.

Eigenvectors are written but not compared.

## Run it

```bash
ctest -R Column_Buckling-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Column_Buckling-4-procs --output-on-failure
```

Both require a build with `MORIS_HAVE_SLEPC`; without it the example compiles no test cases.

## Notes

- Both constitutive models set `Model_Type::PLANE_STRESS`, also when `tDim = 3`.
- The body-fitted deck writes the HMR Lagrange mesh to `LBAProblem.exo`; the immersed deck omits
  this setting.
- The VIS loop appends an extra comma before each eigenpair's field names
  (`Field_Names += ","` followed by `",EVAL..."`), so the generated `Field_Names`,
  `Field_Type` and `IQI_Names` strings contain empty entries between the `VOLUME` field and each
  group.
