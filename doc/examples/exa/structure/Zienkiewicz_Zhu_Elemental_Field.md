# Zienkiewicz–Zhu von Mises stress field on a plate with a hole

A 2D plane-stress plate with a circular hole is clamped on one edge and pulled on the opposite
edge. A nodal von Mises stress field is solved for alongside the displacements, and the
`ZIENKIEWICZ_ZHU_VON_MISES_STRESS` IQI evaluates an elemental Zienkiewicz–Zhu quantity from it.
The deck writes that quantity as an elemental FEM field to HDF5, and the test compares it entry
by entry with a stored reference. The example tests elemental field output from an IQI and the
ZZ stress IQI.

```{exa-summary} projects/EXA/structure/linear/Zienkiewicz_Zhu_Elemental_Field
```

## Problem

| | |
|---|---|
| Domain | square [−2, 2] × [−2, 2] |
| Geometry | level set ‖x‖ − 1.1: material phase `p1` outside a circle of radius 1.1 centered at the origin |
| Material | E = 1.0, ν = 0.0, plane stress |
| Support | `SideSet_4_n_p1` (left edge), displacement (0, 0) by Nitsche |
| Load | traction (0.0, 1.0) on `SideSet_2_n_p1` (right edge) |
| Analysis | one time step over a time frame of 1.0 |

## Setup in MORIS

**Mesh (HMR).** 4 × 4 elements on the 4 × 4 domain (offset (−2, −2)), one level of initial
refinement (`pattern_initial_refinement = 1`), linear Lagrange and B-spline meshes.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `circle`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
off (`tUseGhost = false`).

**Physics (FEM).** `STRUC_LIN_ISO` (`PLANE_STRESS`) on `UX,UY`. IWGs: `STRUC_LINEAR_BULK` on `p1`;
`STRUC_VON_MISES_STRESS` with residual `STRESS_DOF` on `p1`;
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on the left edge (`DIRICHLET_NITSCHE` 100.0);
`STRUC_LINEAR_NEUMANN` on the right edge. IQIs: `DOF` for UX, UY and `STRESS_DOF`, and
`ZIENKIEWICZ_ZHU_VON_MISES_STRESS` (`IQIBulkZienkiewiczZhu`) on `STRESS_DOF` with the structural
constitutive model. A `FEM::FIELDS` entry declares `FieldZienkiewiczZhu` with
`field_entity_type = ELEMENTAL`, `field_type = FIELD_1`, `IQI_Name = IQIBulkZienkiewiczZhu`, and
`field_output_to_file = "Field_Zienkiewicz_Zhu.hdf5"`.

**Solver (SOL).** `AMESOS_IMPL`; Newton with relative residual drop 10⁻¹², relaxation 1.0, at
most 3 iterations; DOFs `STRESS_DOF,UX,UY`, all initialized to 0.0; one time step
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Output (VIS).** `Zienkiewicz_Zhu_Elemental_Field.exo` on `p1`: nodal UX, UY, `Stress` and the
ZZ quantity as `ELEMENTAL_AVG`
([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)).

## What the test checks

In serial only (other processor counts skip the body), the test runs the deck, requires a return
code of 0, and reads the dataset `FieldZienkiewiczZhu` from the output
`Field_Zienkiewicz_Zhu.hdf5` and from the reference
`projects/EXA/structure/linear/Zienkiewicz_Zhu_Elemental_Field/Field_Zienkiewicz_Zhu_Ref.hdf5`.
It checks that both have the same number of entries and that every entry agrees within an
absolute tolerance of 10⁻⁹ (`CHECK`, so all mismatches are reported).

## Run it

```bash
ctest -R Zienkiewicz_Zhu_Elemental_Field-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The OPT list declares `problem = "user_defined"` and the deck defines objective and constraint
  callbacks, but `is_optimization_problem = false`, so they are not used.
- Several deck variables (lengths `sL`, `sL1`, flux `sP2`, conductivity `sK1`, capacity, CTE) are
  defined but not referenced by any parameter list.
- The `.dox` page lists "Imported Elemental Field"; the deck exports the elemental field to HDF5
  and does not import one.
