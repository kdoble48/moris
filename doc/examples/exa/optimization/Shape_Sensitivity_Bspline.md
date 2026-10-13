# Shape sensitivities of analytic and B-spline level sets

A finite-difference check of shape sensitivities for two straight-interface geometries on a
coarse elastic block, in 2D and 3D. Six geometry configurations switch each interface between
fixed, parameterized by an analytic ADV (its center x-coordinate), and discretized on the
B-spline mesh (its B-spline coefficients become the ADVs). For each configuration the SWEEP
algorithm compares the analytical gradients of a strain-energy objective and a volume constraint
with forward, backward and central finite differences.

```{exa-summary} projects/EXA/optimization/Shape_Sensitivity_Bspline
```

## Problem

| | |
|---|---|
| Domain | [−1, 1]² (2D) or [−1, 1]³ (3D), 2 elements per direction |
| Interface 1 | vertical line (2D) or plane (3D) x = 0.2, normal (1, 0[, 0]) |
| Interface 2 | oblique line or plane through (0.65, −1[, −1]) with normal (0.7071, 0.7071[, 0]), i.e. x + y = −0.35 |
| Material phase | phase 0, where both level sets are negative (x < 0.2 and x + y < −0.35) |
| Material | linear elastic, E = 1.0, ν = 0.0 |
| Supports | side set 4 (x = −1), zero displacement, symmetric Nitsche |
| Analysis | static, one time step |
| Objective | strain energy of phase 0 |
| Constraint | volume of phase 0 (density property 1.0) |

## Geometry configurations

The global `tGeoModel` selects the configuration; the deck sets the parameters in
`GENParameterList`.

| `tGeoModel` | Interface 1 (vertical) | Interface 2 (oblique) |
|---|---|---|
| 0 | analytic ADV: `center_x` (0.0, 0.2, 1.0) | fixed |
| 1 | fixed | analytic ADV: `center_x` (0.0, 0.65, 1.0) |
| 2 | analytic ADV | analytic ADV |
| 3 | B-spline coefficients (`discretization_mesh_index = 0`) | fixed |
| 4 | fixed | B-spline coefficients (`discretization_mesh_index = 0`) |
| 5 | B-spline coefficients | B-spline coefficients |

Analytic ADVs are given as (lower bound, initial value, upper bound). The B-spline cases set no
`discretization_lower/upper_bound`, so the parameter-list defaults apply. The comment block at
the top of each deck lists a different case numbering (including a case 6); the table above
follows the code.

## Setup in MORIS

**Mesh (HMR).** 2 × 2 (2D) or 2 × 2 × 2 (3D) elements on a domain of edge 2.0 offset by −1.0,
linear Lagrange and B-spline meshes (pattern 0), refinement and staircase buffers 3,
`adaptive_refinement_level` 1.

**Geometry and phases (GEN).** Two `LINE` geometries in 2D, two `PLANE` geometries in 3D, with
parameters per the table above. The default phase table applies; the physics lives on phase 0.
`IQI_types` is `IQIBulkStrainEnergy`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).
GEN writes its field mesh to `gen_shape_sensitivities.exo` (2D) or
`gen_shape_sensitivities_3D.exo` (3D).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, and Exodus output of the integration mesh.

**Physics (FEM).** A `STRUC_LIN_ISO` constitutive model (E = 1.0, ν = 0.0) and a
`STRUC_LINEAR_BULK` IWG on `HMR_dummy_{n,c}_p0`; a
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` IWG with a `DIRICHLET_NITSCHE` parameter of 100.0 on
`SideSet_4_n_p0` and `SideSet_4_c_p0`. IQIs: `STRAIN_ENERGY`, `VOLUME` weighted by `PropDensity`
(1.0), and nodal displacement `DOF` IQIs (one in 2D, one per component in 3D).

**Solver (SOL).** `AMESOS_IMPL` direct solve, no preconditioner, one Newton iteration, one time
step.

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1`, `include_bounds = false` and
`finite_difference_type = "all"`. The 2D deck uses the default `finite_difference_epsilons`
(10⁻⁸); the 3D deck sets 10⁻⁶. Results go to `shape_opt_test_2D.hdf5` or
`shape_opt_test_3D.hdf5`.

**Output (VIS).** `shape_sensitivities.exo` (2D, field `U`) or `shape_sensitivities_3D.exo`
(3D, fields `UX`, `UY`, `UZ`) on phase 0, written every step.

## What the test checks

Two test cases, `Shape_Sensitivity_Bspline_2D` and `Shape_Sensitivity_Bspline_3D`, each first
delete `*exo*` and `*hdf5*` files in the working directory, then loop over `tGeoModel` = 0 … 5.
Each run must return 0. For every ADV and each FD type (forward, backward, central), the
analytical objective and constraint gradients must equal the FD values under Catch2 `Approx`
with its default relative tolerance and an absolute margin of `MORIS_REAL_EPS` (`CHECK`). The
loops run only when the test is on 1 processor.

## Run it

```bash
ctest -R Shape_Sensitivity_Bspline-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The deck reads the global `tGeoModel`, which the test sets; there are no environment variables.

## Notes

- The traction IWG (`STRUC_LINEAR_NEUMANN`, property `Func_Traction_U`) has an empty
  `mesh_set_names`, so no load set is assigned to it in either deck.
- `PropFlux` (10.0) is defined but not referenced.
