# Thick-walled spherical pressure vessel with a thermal gradient

A thick-walled spherical shell, modeled with three symmetry planes, is loaded by internal and external pressure and by
different temperatures on its inner and outer surfaces; displacement and temperature are solved
together with thermal expansion coupling the two. The shell is an immersed geometry (two sphere
level sets) cut by XTK in both decks. They differ in how the three symmetry planes are
represented: on the HMR domain boundary in `Pressure_Vessel_3D`, or as three additional plane
level sets that cut the background mesh in `Pressure_Vessel_3D_Immersed`. The immersed deck is
also run with quadratic interpolation.

```{exa-summary} projects/EXA/structure/linear/Thick_Pressure_Vessel
```

## Problem

| | |
|---|---|
| Geometry | spherical shell between radii 0.425 and 0.45 (GEN `SPHERE` fields) |
| Material | E = 1.0, ν = 0.3, CTE = 1.0e−5, reference temperature 0.0; conductivity 1.0, density 1.0, heat capacity 0.0 |
| Pressure | −1.0 on the inner surface, −2.0 on the outer surface (`STRUC_LINEAR_NEUMANN`, `Pressure`) |
| Temperature | 100.0 on the inner surface, 200.0 on the outer surface (Nitsche) |
| Symmetry | normal displacement zero on the x, y and z symmetry planes (Nitsche with a `Select` matrix) |
| Analysis | one time step over a time frame of 1.0 |

## Setup in MORIS

The two decks share the FEM, SOL and VIS blocks and differ in mesh, geometry and set names:

| | `Pressure_Vessel_3D` | `Pressure_Vessel_3D_Immersed` |
|---|---|---|
| HMR mesh | 20³ elements on 0.5³ | 24³ elements on 0.6³, offset (−0.1, −0.1, −0.1) |
| GEN geometries | 2 spheres | 2 spheres + 3 planes (normals x, y, z) |
| Vessel bulk set | `HMR_dummy_{n,c}_p1` | `HMR_dummy_{n,c}_p15` |
| x / y / z supports | `SideSet_4`, `SideSet_1`, `SideSet_5` (`_n_p1`, `_c_p1`) | `iside_b0_15_b1_11`, `_13`, `_14` |
| Inner / outer surface | `iside_b0_1_b1_0` / `iside_b0_1_b1_3` | `iside_b0_15_b1_7` / `iside_b0_15_b1_31` |
| Ghost set | `ghost_p1` | `ghost_p15` |
| Level-set clamp | 10⁻⁸ | 10⁻¹⁴ |

**Mesh (HMR).** One Lagrange and one B-spline mesh (pattern 0), both of order
`gInterpolationOrder`, which the test sets to 1 or 2. Refinement and staircase buffer 1.

**Geometry and phases (GEN).** Two `SPHERE` level sets with radii 0.45 and 0.425; the immersed
deck adds three `PLANE` level sets with normals (1, 0, 0), (0, 1, 0), (0, 0, 1). No phase table is
set; the decks address the vessel as phase 1 (two geometries) and phase 15 (five geometries).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`. Supports and loads are applied on XTK interface side sets
(`iside_b0_…`)
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Physics (FEM).** `STRUC_LIN_ISO` (`FULL`) with dependencies `UX,UY,UZ;TEMP` and properties
`YoungsModulus`, `PoissonRatio`, `CTE`, `ReferenceTemperature`; `DIFF_LIN_ISO` on `TEMP`. IWGs:
`STRUC_LINEAR_BULK` and `SPATIALDIFF_BULK` on the vessel;
three `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` symmetry conditions (`DIRICHLET_NITSCHE`
100.0); two `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` temperature conditions (100.0); two
`STRUC_LINEAR_NEUMANN` pressures; `GHOST_NORMAL_FIELD` for displacement and temperature on the
ghost set (`GHOST_DISPL` 0.01 each). IQIs: `DOF` for UX, UY, UZ and TEMP, `STRAIN_ENERGY`, and
`VOLUME` (weighted by density).

**Solver (SOL).** `AMESOS_IMPL`, combined residual/Jacobian assembly, one nonlinear solver on
`UX,UY,UZ;TEMP`, one time step.

**Output (VIS).** `Pressure_Vessel_3D_Case_<gCaseIndex>.exo` on the vessel: nodal UX, UY, UZ,
TEMP and global `STRAIN_ENERGY`, `VOLUME`.

## What the test checks

Three Catch2 cases select the deck, the interpolation order and `gCaseIndex`:

| Case | Deck | Order | Procs | Reference node | Nodes / elements |
|---|---|---|---|---|---|
| 0 | `Pressure_Vessel_3D` | 1 | 1 | 11115 | 17824 / 39168 |
| 1 | `Pressure_Vessel_3D` | 1 | 4 | 8350 | 6061 / 13034 |
| 2 | `Pressure_Vessel_3D_Immersed` | 1 | 1 | 11115 | 17824 / 39168 |
| 3 | `Pressure_Vessel_3D_Immersed` | 1 | 4 | 3250 | 5716 / 12288 |
| 4 | `Pressure_Vessel_3D_Immersed` | 2 | 4 only (skipped in serial) | 3250 | 5716 / 12288 |

For each case the test opens the case's Exodus file (on 4 processors, only rank 1 performs the
checks) and requires 3 dimensions and the node and element counts above (exact), the
reference node's coordinates, the time value 1.0, the displacement vector, the temperature, and
global variable 0 (the strain energy). All value comparisons are relative with tolerance 10⁻⁵
(time: 10⁻⁸). For example, case 0 expects a displacement of
(−1.619383350258344, −1.619383350258375, −1.801581709307490), temperature 200.0945133270288 and
strain energy 0.4911276929614884; case 4 expects temperature 200.3526463545086 and strain energy
0.4916643196542578. The return code of the workflow call is not checked.

## Run it

```bash
ctest -R Thick_Pressure_Vessel-1-procs --output-on-failure     # cases 0 and 2
ctest -R Thick_Pressure_Vessel-4-procs --output-on-failure     # cases 1, 3 and 4
```

Both require a build with `MORIS_USE_EXAMPLES=ON`.

## Notes

- The `.dox` page for this directory has no content beyond its title.
