# Laplace problem on a sliver cut with cell agglomeration

A 2D Laplace problem on the part of the unit square to the left of an immersed vertical plane
that leaves only 0.1 % of each cut cell as material. The deck shows when to turn on cell
agglomeration (`activate_cell_agglomeration`) and what it does: without it, the system matrix of
this problem is indefinite; with it, the matrix is positive definite and well conditioned, with
no ghost penalty. It exercises SPG-based enrichment, cell agglomeration with linear and with
quadratic B-splines, symmetric Nitsche conditions on an immersed boundary, and an analytic
error check.

```{exa-summary} projects/EXA/thermal/diffusion/Laplace_2D_Cell_Agglomeration
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square, 20 × 20 background elements (h = 0.05) |
| Immersed plane | x_c = (10 + α) h with α = 10⁻³, so x_c = 0.50005 |
| Phases | 0 (x < x_c, the computational domain), 1 (x > x_c, no physics) |
| Material | density 1.0, heat capacity 1.0, conductivity 1.0 |
| Bottom and top edges | T = 0 (side ordinals 1 and 3), symmetric Nitsche |
| Immersed plane | T = sin(πy) (`Func_Immersed_Temp`), symmetric Nitsche |
| Left edge | no condition (natural, zero flux) |
| Exact solution | T = sin(πy) cosh(πx) / cosh(πx_c) |

Every cell of the column 0.5 < x < 0.55 is cut and keeps the fraction α of its area as material.

## Why cell agglomeration

The basis functions that are supported only on the cut column have almost no material under
them. Their bulk contributions scale with α, while the symmetric Nitsche terms on the immersed
plane do not. With the Nitsche penalty of this deck (10 k/h) and no stabilization, the system
matrix loses definiteness.

Cell agglomeration attaches each cut cell to an uncut neighbour, its root cell. Every basis
function supported only on cut cells is removed from the unknowns and replaced by the extension
of the root cell's polynomial. Here each of the 20 cut cells is attached to its left neighbour.

Measured with this deck, with `tUseCellAgglomeration` set to false and to true:

| | unknowns | negative eigenvalues | condition number |
|---|---|---|---|
| linear, no stabilization | 252 | 21 | 3.9 × 10³ |
| linear, cell agglomeration | 231 | 0 | 371 |
| quadratic, no stabilization | 286 | 22 | 2.6 × 10¹⁰ |
| quadratic, cell agglomeration | 264 | 0 | 363 |

The face-oriented ghost penalty (`ghost_stab`) is the alternative: it keeps all unknowns and
adds a penalty term with a user-chosen weight. Cell agglomeration has no weight to tune. The
deck leaves the ghost penalty off, so the result is due to the agglomeration alone.

## Setup in MORIS

**Mesh (HMR).** 20 × 20 elements, no refinement, one pattern. The global `gTestIndex` selects the
order of both the Lagrange mesh and the B-spline mesh: 1 for `gTestIndex = 0`, 2 for
`gTestIndex = 1`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, φ = x − x_c (`Func_Plane`), with
multilinear interpolation; `get_phase_index` assigns the positive side to phase 1.

**Cut and enrichment (XTK).** Conformal decomposition, `use_SPG_based_enrichment`, enrichment of
mesh 0, `activate_cell_agglomeration`, ghost stabilization off. `visualize_cell_association`
writes the subphase-group index (`SPG_index0`) and the root index (`Root_Cell_Index_B0`) of every
cell to `xtk_temp_ip.exo`.

**Physics (FEM).** Phase-based input. `DIFF_LIN_ISO` with `SPATIALDIFF_BULK` on phase 0;
`SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE` (penalty 10.0) on side
ordinals 1 and 3 and on the phase-0 side facing phase 1 (`neighbor_phases`). IQIs: the `DOF`
value of `TEMP`, the exact temperature (`PROPERTY`), `L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC`
and `VOLUME`.

**Solver (SOL).** AztecOO GMRES with ILU(0), relative residual 10⁻¹⁰; one Newton iteration. The
condition number of the system matrix is computed (`dense`) and printed to the log, and the
matrix and residual are written with `SOL_save_operator_to_matlab`.

**Output (VIS).** `Laplace_Cell_Agglomeration_Linear.exo` or
`Laplace_Cell_Agglomeration_Quadratic.exo` on the bulk sets of phase 0: nodal `TEMP` and
`TEMP_EXACT`, and global `L2ErrSq`, `H1ErrSq` and `Volume`.

## What the test checks

Two Catch2 cases call the workflow on one processor and require a return code of 0:

| Case | `gTestIndex` | Elements |
|---|---|---|
| `Laplace_Cell_Agglomeration_Linear` | 0 | linear Lagrange on linear B-splines |
| `Laplace_Cell_Agglomeration_Quadratic` | 1 | quadratic Lagrange on quadratic B-splines |

Checked values (`CHECK`, so a failure does not stop the remaining checks):

| | linear | quadratic |
|---|---|---|
| nodes / elements of the output mesh | 354 / 280 | 984 / 280 |
| unknowns (length of the residual file) | 231 | 264 |
| solid cells / cells with a root other than themselves / cells without a root | 220 / 20 / 0 | 220 / 20 / 0 |
| integral of (T_h − T)² | 9.922945826250896 × 10⁻⁸ | 4.681940942909169 × 10⁻¹¹ |
| integral of \|∇T_h − ∇T\|² | 3.780860907617879 × 10⁻³ | 1.351853643804102 × 10⁻⁶ |
| solid volume | 0.50005 | 0.50005 |

The error integrals are compared with a relative tolerance of 10⁻⁶ and the volume with 10⁻¹⁰.
The unknown counts and the root counts fail if the agglomeration does nothing.

## Run it

```bash
ctest -R Laplace_Cell_Agglomeration-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The cell-agglomeration path agglomerates every cut cell. The XTK parameter `volume_fraction`
  is used only by the basis agglomeration (`activate_basis_agglomeration`).
- With cell agglomeration, `Root_Cell_Index_B0` holds subphase-group indices, so a cell is its
  own root when the field equals `SPG_index0`.
- The example runs on one processor only.
