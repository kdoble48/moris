# Property sensitivity sweep

A finite-difference check of adjoint sensitivities with respect to a design-dependent material
property. A uniform 2 × 2 element plate in plane linear elasticity has a SIMP-type Young's
modulus and density driven by a GEN density field whose B-spline coefficients are the ADVs. No
interface is cut, so all sensitivities pass through the PDV `DENSITY`. The OPT `SWEEP`
algorithm compares the adjoint gradients of a strain-energy objective and a volume constraint
against forward, backward and central finite differences. See also
{doc}`One_Element_Mat_Geo_Pdv_sweep <One_Element_Mat_Geo_Pdv_sweep>`.

```{exa-summary} projects/EXA/optimization/Property_Sensitivity_Sweep
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square, 2 × 2 elements |
| Geometry | `Func_Right_Plane`, level set x − 1000; it does not cut the domain, so the whole mesh is phase 0 |
| Material | linear elastic, ν = 0.0, E = 1.0 · ρ³ with ρ the PDV `DENSITY` |
| Supports | displacement (0, 0) on `SideSet_4_n_p0`, symmetric Nitsche (penalty 100) |
| Load | traction (1.0, 0.0) on `SideSet_2_n_p0` |
| Design field | `CONSTANT` 0.4 discretized on B-spline mesh 0, bounds [0.001, 1.0] |
| Objective | strain energy / 2.7 |
| Constraint | density-weighted volume − 0.8 |
| Time | static, one time step over a time frame of 1.0 |

## Setup in MORIS

**Mesh (HMR).** 2 × 2 elements, one linear Lagrange and one linear B-spline mesh on pattern 0
(`lagrange_to_bspline = "0"`), `adaptive_refinement_level = 1`. The MSI parameter list is
empty.

**Geometry and phases (GEN).** One user-defined level set (`Func_Right_Plane`) and one
`CONSTANT` property `density` (value 0.4) with `discretization_mesh_index = 0`; its B-spline
coefficients on mesh 0 are the ADVs. The property supplies the PDV `DENSITY`
(`PDV_types = "DENSITY"`) on `HMR_dummy_n_p0`, `SideSet_4_n_p0` and `SideSet_2_n_p0`
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types`: `IQIBulkStrainEnergy`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, Exodus output of the integration mesh.

**Physics (FEM).** `STRUC_LIN_ISO` with `PropYoungs` (`tYoungsFunc`,
`function_parameters = "1.0,3.0"`) and ν = 0.0; `STRUC_LINEAR_BULK` on `HMR_dummy_n_p0`,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 100) and
`STRUC_LINEAR_NEUMANN`. `PropYoungs` and `PropDensity` both declare
`dv_dependencies = "DENSITY"` with the derivative functions `tDerYoungsFunc` and
`tDerDensityFunc`. IQIs: `STRAIN_ENERGY`, `VOLUME` weighted by `PropDensity`, and `DOF` IQIs
for UX and UY. The computation list sets `finite_difference_scheme = POINT_3_CENTRAL` with
`finite_difference_perturbation_size = 1.0e-7`.

**Solver (SOL).** `AMESOS_IMPL` with default settings, one nonlinear algorithm with
`NLA_combined_res_jac_assembly = true`, one nonlinear solver on UX,UY, and
`Sensitivity_Analysis_Type = ADJOINT`.

**Optimization (OPT).** `user_defined` problem with the `SWEEP` algorithm:
`hdf5_path = "property_opt_test.hdf5"`, objective and constraint gradients evaluated,
`num_evaluations_per_adv = "1"`, `include_bounds = false`, `finite_difference_type = "all"`,
`finite_difference_epsilons = "1.0e-7"`. The callbacks provide `dcriteria` derivatives
(1/2.7 for the objective, 1.0 for the constraint); the `dadv` derivatives are zero
([coupling C2](../../../user-guide/interlinks.md#c2-problem-declaration-opt)).

**Output (VIS).** `Property_Sensitivity_Sweep.exo` on `HMR_dummy_n_p0` with nodal UX and UY.

## What the test checks

The test runs `Property_Sensitivity_Sweep.so`, requires a return code of 0, and opens
`property_opt_test.hdf5`. On 1 processor it reads the analytic gradients
(`objective_gradients eval_1-1 analytical`, `constraint_gradients eval_1-1 analytical`) and the
`fd_forward`, `fd_backward` and `fd_central` gradients at `epsilon_1-1`. For every ADV and each
FD scheme it checks `analytical == Approx( fd ).margin( 1E-16 )` for the objective and for the
constraint, i.e. Catch2's default relative tolerance of `Approx` with an absolute margin of
10⁻¹⁶. It also requires matching gradient vector lengths.

## Run it

```bash
ctest -R Property_Sensitivity_Sweep-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The directory has no `example_documentation.dox`.
- `tDerYoungsFunc` returns β · p · ρ^(β−1) (p the parameter vector) without the shape-function
  factor `N()` that `tDerDensityFunc` returns.
- IQIs of type `NODE_SENSITIVITY` and `DOF_SENSITIVITY` are present in the deck but commented
  out.
