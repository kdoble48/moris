# Homogenization of a 3D unit cell with a spherical inclusion

A periodic 3D unit cell — a stiff spherical inclusion embedded in a softer matrix — is loaded by a
unit eigenstrain to compute a homogenized constitutive coefficient. The deck is the 3D version of
{doc}`Homogenization_2D <Homogenization_2D>`: the inclusion is an immersed level-set geometry cut
by XTK, periodicity is enforced weakly across the three pairs of opposite cell faces through MIG,
and the homogenized value is integrated by the `HOMOGENIZED_CONSTITUTIVE` IQI.

```{exa-summary} projects/EXA/structure/linear/Homogenization_3D
```

## Problem

| | |
|---|---|
| Domain | unit cube 1.0 × 1.0 × 1.0 |
| Phases | `p1` inclusion (sphere, center (0.5, 0.5, 0.5), radius 0.25), `p0` matrix |
| Inclusion | E = 2.0, ν = 0.1, full 3D model |
| Matrix | E = 1.0, ν = 0.1, full 3D model |
| Load | eigenstrain (1.0, 0, 0, 0, 0, 0) in both phases |
| Supports | displacement (0, 0, 0) on `SideSet_1_n_p0`, `SideSet_4_n_p0`, `SideSet_6_n_p0`, active only where x, y and z are all < 0.1 |
| Periodicity | side-set pairs `4,2;1,3;5,6` |
| Analysis | static, linear |

## Setup in MORIS

The module configuration matches {doc}`Homogenization_2D <Homogenization_2D>` block for block; the
differences are:

| | 2D | 3D |
|---|---|---|
| HMR mesh | 10 × 10 on 1 × 1 | 5 × 5 × 5 on 1 × 1 × 1 |
| Inclusion level set | circle, r = 0.25 | sphere, r = 0.25 |
| Constitutive `model_type` | `PLANE_STRESS` | `FULL` |
| Eigenstrain | 3 components | 6 components |
| Periodic pairs (MIG) | `4,2;1,3` | `4,2;1,3;5,6` |
| Dirichlet side sets | 1, 4 | 1, 4, 6 |
| DOFs | `UX,UY` | `UX,UY,UZ` |

**Mesh (HMR).** Linear Lagrange and linear B-spline meshes (pattern 0), refinement and staircase
buffer 1.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Inclusion`, returning
0.25 − ‖x − (0.5, 0.5, 0.5)‖, clamped to 10⁻⁸ near zero.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`. FEM consumes the phase bulk sets, the interface
`dbl_iside_p0_1_p1_0`, and `ghost_p0`, `ghost_p1`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Physics (FEM).** Per-phase `STRUC_LINEAR_BULK`; a Nitsche Dirichlet IWG
(`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE`, `DIRICHLET_NITSCHE` 100.0) whose `Select` property
(`Func_Dirichlet`) is the identity only for x, y, z < 0.1; the inclusion–matrix interface and the
periodic set `P00` with `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE` 100.0
and 10000.0); `GHOST_NORMAL_FIELD` on both ghost sets (`GHOST_DISPL` 0.01). IQIs: `DOF` for UX
and UY, `PROPERTY` for the two moduli, `HOMOGENIZED_CONSTITUTIVE` per phase.

**Solver (SOL).** `AMESOS_IMPL`, default Newton, DOFs `UX,UY,UZ`; the operator is written to
`jp2.dat` (`SOL_save_operator_to_matlab`).

**Output (VIS).** `Homogenization_3D.exo` on both phases: nodal UX, UY and the two moduli, global
`IQIIn` and `IQIOut`.

## What the test checks

The test requires a return code of 0 and, in serial only, opens `Homogenization_3D.exo` and
requires 3 dimensions, 976 nodes and 1731 elements; nodes 33 and 584 at (0.4, 0.6, 0.4) and
(0.2, 0.8, 1.0) (relative error below 10⁻²); and the sum of global variables 0 and 1 equal to
1.069 within an absolute tolerance of 10⁻².

## Run it

The ctest name is exactly `Homogenization_3D` (`NO_PROCS_SUFFIX`). The example is not currently
built: `add_subdirectory(Homogenization_3D)` is commented out in
`projects/EXA/structure/linear/CMakeLists.txt`. After re-enabling it, in a build with
`MORIS_USE_EXAMPLES=ON`:

```bash
ctest -R Homogenization_3D --output-on-failure
```

## Notes

- As in the 2D deck, the parameter lists use the pre-`Module_Parameter_Lists` API, and
  `SPPeriodicNitscheP11` (1000.0) is defined but not referenced by any IWG.
- The `.dox` page for this directory carries the 2D page identifier and title
  (`EXA_Structure_Linear_Homogenization_2D`) and has no content beyond the feature name.
