# Homogenization of a 2D unit cell with a circular inclusion

A periodic 2D unit cell — a stiff circular inclusion embedded in a softer matrix — is loaded by a
unit eigenstrain to compute a homogenized constitutive coefficient. The inclusion is an immersed
level-set geometry cut by XTK, periodicity is enforced weakly across opposite cell faces through
MIG side-set pairs, and the homogenized value is integrated by the `HOMOGENIZED_CONSTITUTIVE` IQI.
The 3D counterpart is {doc}`Homogenization_3D <Homogenization_3D>`.

```{exa-summary} projects/EXA/structure/linear/Homogenization_2D
```

## Problem

| | |
|---|---|
| Domain | unit square 1.0 × 1.0 |
| Phases | `p1` inclusion (circle, center (0.5, 0.5), radius 0.25), `p0` matrix |
| Inclusion | E = 2.0, ν = 0.1, plane stress |
| Matrix | E = 1.0, ν = 0.1, plane stress |
| Load | eigenstrain (1.0, 0.0, 0.0) in both phases |
| Supports | displacement (0, 0) on `SideSet_1_n_p0` and `SideSet_4_n_p0`, active only where x < 0.1 and y < 0.1 |
| Periodicity | side-set pairs `4,2;1,3` |
| Analysis | static, linear |

## Setup in MORIS

**Mesh (HMR).** 10 × 10 elements on the unit square, linear Lagrange and linear B-spline meshes
(pattern 0), refinement and staircase buffer 1.

**Geometry and phases (GEN).** One `USER_DEFINED` level-set geometry, `Inclusion`, returning
0.25 − ‖x − (0.5, 0.5)‖ with values smaller than 10⁻⁸ in magnitude clamped to 10⁻⁸. No design
variables and no phase table are set.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`, Exodus output of the integration mesh. The deck consumes
the bulk sets `HMR_dummy_{n,c}_p0/p1`, the interface `dbl_iside_p0_1_p1_0` and the ghost sets
`ghost_p0`, `ghost_p1`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Periodic boundaries (MIG).** `periodic_side_set_pair = "4,2;1,3"` pairs the left/right and
bottom/top faces; the periodic IWG acts on the resulting set `P00`.

**Physics (FEM).** Four `STRUC_LIN_ISO` constitutive models (`PLANE_STRESS`), each carrying
`YoungsModulus`, `PoissonRatio` and `EigenStrain`: one per phase plus two periodic copies. IWGs:

- `STRUC_LINEAR_BULK` on each phase;
- `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on the support sets, with a `Select` property
  (`Func_Dirichlet`) that switches the constraint on only for x < 0.1 and y < 0.1;
  `DIRICHLET_NITSCHE` parameter 100.0;
- `STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE` on the inclusion–matrix interface
  (`NITSCHE_INTERFACE`, 100.0) and on `P00` (`NITSCHE_INTERFACE`, 10000.0, matrix on both sides);
- `GHOST_NORMAL_FIELD` on `ghost_p0` and `ghost_p1` (`GHOST_DISPL`, 0.01).

IQIs: `DOF` for UX and UY, `PROPERTY` for the two Young's moduli, and
`HOMOGENIZED_CONSTITUTIVE` on each phase (`IQIHomInner`, `IQIHomOuter`). A third interface
parameter, `SPPeriodicNitscheP11` (1000.0), is defined but not referenced by any IWG.

**Solver (SOL).** `AMESOS_IMPL` direct solver, default Newton settings, DOFs `UX,UY`. The solver
warehouse sets `SOL_save_operator_to_matlab = "jp2.dat"`.

**Output (VIS).** `Homogenization_2D.exo` on both phases: nodal UX, UY and the two moduli, and the
two homogenized IQIs as global variables `IQIIn` and `IQIOut`.

## What the test checks

The test requires a return code of 0 and, in serial only, opens `Homogenization_2D.exo` and
requires 2 dimensions, 337 nodes and 256 elements; nodes 10 and 140 at (0.5, 0.5) and (0.5, 0.1)
(relative error below 10⁻²); and the sum of global variables 0 and 1 — the two
`HOMOGENIZED_CONSTITUTIVE` contributions — equal to 1.14 within an absolute tolerance of 10⁻².

## Run it

The example is registered with `NO_PROCS_SUFFIX`, so its ctest name is exactly
`Homogenization_2D`. It is not currently built: `add_subdirectory(Homogenization_2D)` is commented
out in `projects/EXA/structure/linear/CMakeLists.txt`. After re-enabling that line, in a build
with `MORIS_USE_EXAMPLES=ON`:

```bash
ctest -R Homogenization_2D --output-on-failure
```

## Notes

- The leaf `CMakeLists.txt` records that the deck predates the current parameter-list API (it
  calls `aParameterLists( i ).push_back( prm::create_..._parameter_list() )` directly), which is
  why the parent leaves it disconnected.
- `tOuterSurface` (side sets 1–4) and the density properties (0.0) are defined in the deck but
  not used by any IWG.
