# SIMP density optimization of a 2D cantilever

A density-based (SIMP) topology optimization of a 2D cantilever under a tip load, run for two
GCMMA iterations. The design variables are the B-spline coefficients of a density field that GEN
publishes as a `DENSITY` PDV; FEM properties read that PDV to penalize the Young's modulus and
weight the volume. The deck exercises the property pathway by which ADVs reach the physics
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)) on an
uncut mesh.

```{exa-summary} projects/EXA/optimization/SIMP
```

## Problem

| | |
|---|---|
| Domain | 2.0 × 1.0 rectangle |
| Phases | one material phase (phase 0); the level set x − 1000 is negative everywhere, so nothing is cut |
| Material | linear elastic, isotropic, E(ρ) = ρ³ (`tYoungsFunc`, parameters 1.0, 3.0), ν = 0.0 |
| Design field | density ρ, initial value 0.4, bounds [0.001, 1.0] |
| Supports | left edge (side set 4), zero displacement, symmetric Nitsche |
| Load | traction (0, 1.0) on the right edge (side set 2), applied only below y = 0.2 |
| Objective | strain energy / 2.7 |
| Constraint | density-weighted volume − 0.8 (constraint type 1) |
| Optimizer | GCMMA, 2 iterations (`max_its = 2`) |

## Setup in MORIS

**Mesh (HMR).** 40 × 20 elements, linear Lagrange and linear B-spline meshes (pattern 0),
refinement and staircase buffers 3, `adaptive_refinement_level` 1, HMR `severity_level` 1.

**Geometry and phases (GEN).** One `USER_DEFINED` level-set geometry, `Func_Right_Plane`
(x − 1000), which places the whole domain in phase 0. The design enters as a GEN property: a
`CONSTANT` field named `density` with value 0.4, discretized on B-spline mesh 0
(`discretization_mesh_index = 0`, bounds 0.001 and 1.0), so its coefficients are the ADVs. It
publishes `pdv_type = DENSITY` on `HMR_dummy_n_p0`, `SideSet_4_n_p0` and `SideSet_2_n_p0`
(`PDV_types = "DENSITY"`). `IQI_types` is `IQIBulkStrainEnergy`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, and Exodus output of the integration mesh. No interface lies in the domain.

**Physics (FEM).** Two design-dependent properties (`dv_dependencies = "DENSITY"`):
`PropYoungs`, value function `tYoungsFunc` (parameters times ρ raised to the second parameter,
3.0) with derivative `tDerYoungsFunc`, and `PropDensity`, value ρ (`tDensityFunc`) with
derivative `tDerDensityFunc`. The `STRUC_LIN_ISO` constitutive model takes `PropYoungs` and
`PropPoisson` (0.0). IWGs: `STRUC_LINEAR_BULK` on `HMR_dummy_n_p0`,
`STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE` on `SideSet_4_n_p0` with a `DIRICHLET_NITSCHE`
parameter of 100.0, and `STRUC_LINEAR_NEUMANN` on `SideSet_2_n_p0` with the traction
`Func_Neumann_U`. IQIs: `STRAIN_ENERGY`, `VOLUME` weighted by `PropDensity`, and a `DOF` IQI.

**Solver (SOL).** `BELOS_IMPL` with an `IFPACK` ILU preconditioner, default Newton settings, one
time step.

**Optimization (OPT).** `GCMMA` with `max_its = 2` and otherwise default parameters. The
objective is criterion 0 divided by 2.7; the constraint is criterion 1 minus 0.8.

**Output (VIS).** `SIMP.exo` on `HMR_dummy_n_p0` with the nodal field `VOLUME` from
`IQIBulkVolume`.

## What the test checks

The test requires a return code of 0, then reads `SIMP.exo` and checks the third nodal field
(index 2) at the first time step for nodes 0, 100, …, 800 against reference values with Catch2
`Approx` default tolerance:

| Node | 0 | 100 | 200 | 300 | 400 | 500 | 600 | 700 | 800 |
|---|---|---|---|---|---|---|---|---|---|
| Value | 0.419977 | 0.419974 | 0.380038 | 0.397622 | 0.380049 | 0.380093 | 0.380112 | 0.419978 | 0.419973 |

It runs on 1 processor.

## Run it

```bash
ctest -R '^SIMP$' --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is exactly `SIMP` (no processor suffix); the anchored pattern avoids matching
other tests whose names contain "SIMP".

## Notes

- `tYoungsFunc` multiplies the whole parameter row (1.0, 3.0) by ρ³; the exponent is read from
  the second entry.
- `PropFlux` (10.0) is defined but not referenced by any IWG or IQI.
- The CMake comment records that the legacy test was also gated on `MORIS_HAVE_GCMMA`;
  `moris_add_example` registers it under `MORIS_HAVE_PARALLEL_TESTS` only.
