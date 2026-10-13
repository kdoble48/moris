# Contact sensitivity test

A finite-difference sensitivity check for linear elasticity with a contact
interface between two materials. A tilted line splits a unit square (cube) with a central hole
into two bodies; one is supported, the other is pushed toward it by a prescribed
displacement. The interface position is the design. A `SWEEP` compares adjoint gradients of
strain energies, volume, and interface perimeter with central finite differences, for an
analytic interface parameter and for a B-spline-discretized interface field.

```{exa-summary} projects/EXA/optimization/Contact_Sensitivity_Test
```

## Problem

| | |
|---|---|
| Domain | unit square (2D) or unit cube (3D), selected by `gDim` |
| Hole | circle (2D) or sphere (3D) of radius 0.166 centred at (0.5, 0.5[, 0.5]); void |
| Interface | `LINE` level set through (0.638, 0) with normal (−sin 75°, −cos 75°), i.e. tilted 75° from horizontal |
| Phases | phase table `[0, 0, 1, 2]`: phase 0 void, phase 1 material 2, phase 2 material 1 |
| Material 1 (p2) | linear elastic, E = 2.0, ν = 0.0 |
| Material 2 (p1) | linear elastic, E = 1.0, ν = 0.3 |
| Bedding | 10⁻⁶ on both materials |
| Supports | material 1 on side set 4 (left edge), zero displacement, Nitsche-enforced |
| Load | material 2 on side set 2 (right edge), prescribed displacement (−0.1, 0[, 0]), Nitsche-enforced |
| Interface condition | `STRUC_LINEAR_CONTACT_UNSYMMETRIC_NITSCHE` between material 1 and material 2 |
| Criteria | strain energy of material 1, strain energy of material 2, volume of material 1, interface perimeter |

## Setup in MORIS

**Mesh (HMR).** One element per direction on a domain of edge 1.0, with two refinement
patterns: pattern 0 (level set) refined once and pattern 1 (displacement) refined three times.
B-spline orders come from `gLevelSetInterpolationOrder` and `gFEMInterpolationOrder`, the
Lagrange order from `gLagrMeshInterpolationOrder`; the test sets all three to 1. The Lagrange
mesh uses pattern 1. Displacements live on B-spline mesh 1 (`UX`, `UY`[, `UZ`] set to 1 in MSI).

**Geometry and phases (GEN).** Geometry 0 is the hole (`CIRCLE` or `SPHERE`, no ADVs). Geometry
1 is the interface `LINE`. Its parameterization depends on `gUseBspline`:

| `gUseBspline` | Interface ADVs | Level-set property `LevelsetField` | FD ADV indices |
|---|---|---|---|
| `false` (analytic) | `center_x` as (0.319, 0.638, 1.276) | `CONSTANT` 1.0, discretized on mesh 0 (bounds ±2) | `1` |
| `true` (B-spline) | line field discretized on mesh 0 (bounds ±2) | `SCALED_FIELD` of the line field, factor 1.0 | `3,4,5` (2D), `7,15,21` (3D) |

In both cases `LevelsetField` publishes PDV type `LS1` on all bulk sets
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types` is `IQIBulkStrainEnergy_Material1`, `IQIBulkStrainEnergy_Material2`,
`IQIBulkVolume_Material1`, `IQIPerimeter_InterfaceMaterial12`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization with ghost visualization, Exodus output of the integration mesh.

**Physics (FEM).** Two `STRUC_LIN_ISO` models (material 1 on `HMR_dummy_*_p2`, material 2 on
`HMR_dummy_*_p1`), each with a `STRUC_LINEAR_BULK` IWG carrying the bedding property. Support
and load use `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` with a `DIRICHLET_NITSCHE` penalty of
10.0. The contact IWG acts on `dbl_iside_p0_2_p1_1` with a `NITSCHE_INTERFACE` penalty of 10.0
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
Ghost terms use `GHOST_NORMAL_FIELD` with `GHOST_DISPL` 0.005 on `ghost_p1` and `ghost_p2`.
The FEM property `PropLevelSet` returns the `LS1` field value and has a DV derivative; the
volume IQIs use it as their `Density`, and a `PROPERTY` IQI outputs it. The perimeter is a
`VOLUME` IQI on `iside_b0_2_b1_1`. FEM-internal finite differences use `POINT_3_CENTRAL` with a
perturbation of 10⁻⁵.

**Solver (SOL).** Direct `AMESOS_IMPL` solve with `Amesos_Mumps` when `MORIS_USE_MUMPS` is
defined and `Amesos_Superludist` otherwise, no preconditioner. Newton runs up to 20 iterations
to a relative residual drop of 10⁻⁹.

**Optimization (OPT).** `SWEEP` with `finite_difference_type = "all"`,
`finite_difference_epsilons = 1.0e-5`, one evaluation per ADV, and
`finite_difference_adv_indices` restricted as in the table above. The objective is the sum of
the four criteria, and each criterion is also a constraint. Results go to
`Contact_Sensitivity_Test_<case>_SEN.hdf5`.

**Output (VIS).** `Contact_Sensitivity_Test_<case>.exo` on the bulk, boundary, and interface
sets, with nodal displacements, the level-set property (2D), and the criteria as global fields.

## What the test checks

For each case the test reads `constraint_gradients eval_1-1 analytical` and
`constraint_gradients eval_1-1 epsilon_1-1 fd_central`, requires equal sizes, and checks that
every entry has a relative difference below 10⁻³. No Exodus or criterion values are compared.

| Case (1 proc) | Case (4 procs) | Dimension | Interface |
|---|---|---|---|
| 0 | 10 | 2D | analytic |
| 1 | 11 | 2D | B-spline |
| 3 | 13 | 3D | analytic |
| 4 | 14 | 3D | B-spline |

CMake registers only the 1-processor run, so cases 10–14 are not exercised by ctest.

## Run it

```bash
ctest -R Contact_Sensitivity_Test --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is `Contact_Sensitivity_Test-1-procs`.

## Notes

- `tPrescribedTraction = false` and `tUseContact = true` are hard-coded. Setting them the other
  way switches the load to a `STRUC_LINEAR_NEUMANN` traction (1, 0[, 0]) and the interface to
  `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE`.
- In 3D the interface is still a `LINE` field with only `normal_x` and `normal_y` set.
- The 2D VIS `Field_Names` string contains an empty entry (`"UX,UY,,Levelset,..."`), giving
  eight names for seven field types and IQIs; the field named `VolumeMaterial2` is wired to
  `IQIBulkVolume_Material1`.
- The user function `Plane_2D3D` is defined but not referenced by any parameter list.
- {doc}`Bedding_Sensitivity_Test <Bedding_Sensitivity_Test>` uses the same check for a bedding
  term.
