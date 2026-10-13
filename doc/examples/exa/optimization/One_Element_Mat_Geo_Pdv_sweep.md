# One-element material and geometry PDV sweep

A finite-difference sensitivity check on a single cut element. A straight-line level set cuts
one bilinear element into a material and a void part; the material carries a SIMP-type Young's
modulus and density that depend on a GEN property field through the PDV `LS1`. The OPT
`SWEEP` algorithm compares the analytic gradients of a strain-energy objective and a
volume constraint with respect to the line parameters against forward, backward and central
finite differences. The deck tests that geometric (interface position) and material (property
PDV) contributions to the sensitivities are combined consistently.

```{exa-summary} projects/EXA/optimization/One_Element_Mat_Geo_Pdv_sweep
```

## Problem

| | |
|---|---|
| Domain | 10.0 × 10.0 square, offset (−5, −5), one element |
| Geometry | `LINE` level set through (2.0, 0.0) with normal (1.0, 0.0), `isocontour_threshold = 0.5` |
| Phases | 0 material (bulk sets `HMR_dummy_*_p0`), 1 void |
| Material | linear elastic, ν = 0.0, E = 2.0 · ρ³; ρ from the PDV `LS1` (below) |
| Supports | displacement (0, 0) on `SideSet_4_n_p0,SideSet_4_c_p0`, symmetric Nitsche (penalty 100) |
| Load | traction (0, 1.0) on the material–void interface `iside_b0_0_b1_1` (`Func_Neumann_U`; its y < 20 condition holds everywhere in the domain) |
| Objective | strain energy / 1437.01 |
| Constraint | density-weighted volume / 40 − 1 ≤ 0 |

The PDV `LS1` comes from a `CONSTANT` GEN property `LvL_Set_Field` with value 0.8 on the
material bulk sets. In `tYoungsFunc` and `tDensityFunc` it is mapped to
ρ = (LS1 − 0.5)/(1 − 0.5), clamped to [10⁻⁴, 1].

## Setup in MORIS

**Mesh (HMR).** 1 × 1 element, two linear Lagrange and two linear B-spline meshes on patterns
0 and 1 (`lagrange_to_bspline = "0,1;-1"`), no initial refinement. MSI places UX, UY on
B-spline mesh 0.

**Geometry and phases (GEN).** One `LINE` geometry whose four parameters are given as
(lower, initial, upper) triples: `center_x` (−5, 2, 5), `center_y` (−1, 0, 1), `normal_x`
(−2, 1, 2), `normal_y` (−1, 0, 1); these are the ADVs. `PDV_types = "LS1"`. The `CONSTANT`
property is given as a single value 0.8 with `discretization_lower_bound = 0.001` and
`discretization_upper_bound = 1.0`, and supplies `LS1` on `HMR_dummy_n_p0`, `HMR_dummy_c_p0`
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types`: `IQIBulkStrainEnergy`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization off (`tUseGhost = false`), Exodus output of the integration mesh.

**Physics (FEM).** `STRUC_LIN_ISO` with `PropYoungs` (`tYoungsFunc`,
`function_parameters = "2.0,3.0"`, `dv_dependencies = "LS1"`), `STRUC_LINEAR_BULK` on the
material, `STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 100) and
`STRUC_LINEAR_NEUMANN` on the interface. IQIs: `STRAIN_ENERGY` and `VOLUME` weighted by
`PropDensity` (`tDensityFunc`, `dv_dependencies = "LS1"`) on the material, plus `DOF` IQIs for
UX and UY.

**Solver (SOL).** `AMESOS_IMPL` with its default solver type, one Newton solver (up to 5
iterations, relative drop 10⁻⁷) on UX,UY, one time solver. The solver warehouse sets
`SOL_save_operator_to_matlab = "Mat.dat"`.

**Optimization (OPT).** `user_defined` problem with the `SWEEP` algorithm:
`hdf5_path = "shape_opt_test.hdf5"`, `num_evaluations_per_adv = "1"`,
`finite_difference_type = "all"`, `finite_difference_epsilons = "1e-6"`. The objective and
constraint callbacks provide `dcriteria` derivatives; the `dadv` derivatives are zero
([coupling C2](../../../user-guide/interlinks.md#c2-problem-declaration-opt)).

**Output (VIS).** `One_Element_Mat_Geo_Pdv_sweep.exo` on the material bulk sets with nodal UX,
UY and VOL, `Time_Offset = 10.0`.

## What the test checks

The test runs `One_Element_Mat_Geo_Pdv_sweep.so`, requires a return code of 0, and opens
`shape_opt_test.hdf5`. On 1 processor it reads the analytic objective and constraint gradients
(`objective_gradients eval_1-1 analytical`, `constraint_gradients eval_1-1 analytical`) and the
gradients for `fd_forward`, `fd_backward` and `fd_central` at `epsilon_1-1`. For every ADV and
each FD scheme it requires, separately for objective and constraint:

- relative difference |analytic − FD| / (|FD| + 10⁻¹⁴) < 10⁻⁵ when |FD| > 10⁻⁶;
- otherwise agreement within an absolute margin of 10⁻⁶.

It also requires that the objective and constraint gradient vectors have the same length as
each other and as their FD counterparts.

## Run it

```bash
ctest -R One_Element_Mat_Geo_Pdv_sweep-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- `tDerYoungsFunc` and `tDerDensityFunc`, registered as `dv_derivative_functions`, abort with
  `MORIS_ERROR( false, "Do not need this one" )` if called.
- The deck defines `Const_Geometry` (|x| − 3) but no parameter list references it.
- The `.dox` page contains only a "TODO: Add features here" placeholder.
