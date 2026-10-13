# Two-material bar with an immersed interface

A steady 1D heat-conduction problem written as a 2D strip or a 3D bar: two materials with
different conductivities and heat sources meet at an immersed interface at x = 4.1. The left end
carries the exact temperature, the right end a heat flux, and the deck codes the piecewise
quadratic analytic solution to evaluate L2 and H1 semi-norm errors. The same deck runs in 2D and
3D (global `gDim`) with linear or quadratic interpolation (global `gInterpolationOrder`). It
exercises a material interface cut through elements and coupled by symmetric Nitsche's method,
and writes XTK basis-function and enrichment fields that the debug build checks.

```{exa-summary} projects/EXA/thermal/diffusion/Two_Material_Bar
```

## Problem

| | |
|---|---|
| Domain | length L = 10.0; cross-section 1.0 (2D) or 1.0 × 1.0 (3D) |
| Phases | 1 material 1 (x < 4.1), 0 material 2 (x > 4.1) |
| Conductivity | k₁ = 1.0, k₂ = 0.125 |
| Heat source | Q₁ = 1.0, Q₂ = 2.0 |
| Density, heat capacity | 0.0 in both materials |
| Left end (`SideSet_4`) | exact temperature (T = 0 at x = 0), Nitsche-enforced |
| Right end (`SideSet_2`) | heat flux 5.0 (Neumann) |
| Time | one time step, time frame 1.0 |

`Func_Exact_Temperature` gives, for x < L₁ = 4.1,
T = T_pre − Q₁x²/(2k₁) + x (P + L₁Q₁ + L₂Q₂)/k₁ with T_pre = 0, P = 5.0 and L₂ = L − L₁, and a
quadratic in x with the material-2 constants for x > L₁ (the full expression is in the deck).
The directory also holds `HeatedBarBenchmark.m`, `Two_Material_Bar_Linear.m` and
`Two_Material_Bar_Quadratic.m`.

## Setup in MORIS

**Mesh (HMR).** 2 × 1 (2D) or 2 × 1 × 1 (3D) elements, one level of initial refinement,
refinement and staircase buffer 1. Lagrange and B-spline orders follow `gInterpolationOrder`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Inclusion` = 4.1 − x (values with
magnitude below 10⁻⁸ replaced by 10⁻⁸), with the default phase assignment.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
off. Output: integration and interpolation meshes (`exodus_output_XTK_ig_mesh`,
`exodus_output_XTK_ip_mesh`), basis functions (`write_basis_functions`) and enrichment fields with
the probe sphere "5.0,0.0,0.0,0.0".

**Physics (FEM).** Legacy set-name input
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
`DIFF_LIN_ISO` per material with `SPATIALDIFF_BULK` on `HMR_dummy_{n,c}_p1` and
`HMR_dummy_{n,c}_p0`; `SPATIALDIFF_NEUMANN` (flux 5.0) on `SideSet_2_n_p0`;
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, penalty 100.0) with the exact
temperature on `SideSet_4_n_p1`; `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`,
penalty 100.0) on `dbl_iside_p0_1_p1_0`. Ghost IWGs are written under `if ( tUseGhost )` and are
inactive. IQIs: `DOF` temperature, `PROPERTY` analytic temperature, `L2_ERROR_ANALYTIC`,
`H1_ERROR_ANALYTIC`, and four `VOLUME` IQIs on the bulk sets and on `iside_b0_1_b1_0`.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type, Newton with at most 2 iterations
and relative residual drop 10⁻⁸, one time step.

**Output (VIS).** `Two_Material_Bar_<gTestCaseIndex>.exo`: nodal `TEMP` and `TEMP_ANALYTIC`,
global `L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC` and `VOLUME`, and the volume IQIs written as
`ELEMENTAL_AVG`, `ELEMENTAL_INT`, `FACETED_AVG` and `FACETED_INT` fields.

## Test cases

| `gTestCaseIndex` | Catch2 case | `gDim` | Order | Nodes / elements | Reference node | Temperature |
|---|---|---|---|---|---|---|
| 0 | `Two_Material_Bar_Linear` | 2 | 1 | 45 / 26 | 41 at (4.1, 0.5) | 77.31885467530020 |
| 1 | `Two_Material_Bar_Linear` | 3 | 1 | 235 / 332 | 72 at (4.1, 0.5, 0.5) | 77.29432195240925 |
| 2 | `Two_Material_Bar_Quadratic` | 2 | 2 | 70 / 26 | 66 at (4.1, 0.5) | 77.28500000000312 |
| 3 | `Two_Material_Bar_Quadratic` | 3 | 2 | 390 / 332 | 129 at (4.1, 0.5, 0.5) | 77.28500000001044 |

## What the test checks

Each Catch2 case runs the 2D and then the 3D variant; the return codes are not checked. For each
output file the test requires the dimension and the node and element counts above, the reference
node coordinates (relative tolerance 10⁻⁸), a time value of 1.0 (10⁻⁸), and the nodal temperature
at the reference node (relative tolerance 10⁻⁴). It also compares global variables 0 and 1 with
an absolute tolerance of 10⁻⁶:

| Case | Global variable 0 | Global variable 1 |
|---|---|---|
| 0 | 348.1332925575233 | 683.8626208796658 |
| 1 | 347.9932121344261 | 683.8626251218382 |
| 2 | 8.922568376614086 × 10⁻²² | 2.911872014125403 × 10⁻²³ |
| 3 | 5.762418056802439 × 10⁻²⁰ | 2.805634603810048 × 10⁻²¹ |

In a build with `MORIS_HAVE_DEBUG`, case 2 additionally reads `xtk_temp.exo` and `xtk_temp_ip.exo`:
54 nodes / 26 elements and 60 nodes / 10 elements, node coordinates (1.25, 0.25) and (1.25, 0.0),
and exact equality of the enrichment weights `weights_mi_0_ind_4`, `weights_mi_0_ind_22`
(1.0, −10.0) and of `Basis_mi_0_ind_0`, `Basis_mi_0_ind_4`, `weights_mi_0_ind_10`
(0.0625, 0.375, 0.0625).

## Run it

```bash
ctest -R Two_Material_Bar-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The reference values of global variables 0 and 1 for the quadratic cases are of order 10⁻²⁰,
  so the absolute tolerance of 10⁻⁶ there checks only that they are near zero. VIS lists
  `L2_ERROR_ANALYTIC` and `H1_ERROR_ANALYTIC` as the first two global fields.
- The `.dox` documentation file is absent from this directory.
