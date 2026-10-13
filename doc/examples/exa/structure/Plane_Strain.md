# Plane strain block under uniaxial load

A four-element first-order mesh, read from an Exodus file, is compressed by a uniform
traction on its top edge while roller-type supports on the left and bottom edges leave the
in-plane Poisson contraction unrestrained. The example tests the `PLANE_STRAIN` model type of the
linear isotropic constitutive model, including the out-of-plane normal stress it produces.

```{exa-summary} projects/EXA/structure/linear/Plane_Strain
```

## Problem

| | |
|---|---|
| Mesh | `4elements.e`: 2D, 9 nodes, 4 elements, block `block_main` |
| Material | linear elastic, isotropic: E = 1.0, ν = 0.3, plane strain |
| Supports | `sideset_left`: UX = 0; `sideset_bottom`: UY = 0 (Nitsche, component selected) |
| Load | traction (0.0, −1.0) on `sideset_top` |
| Thickness | property value 0.1 on all IWGs |
| Analysis | static, one solve |

The `.dox` page cites a plane-strain reference text for the solution. The test compares against
σ_x = 0, σ_y = −1.0 and an out-of-plane stress σ_z = −0.3.

## Setup in MORIS

**Mesh (STK).** Workflow `STK_FEM`; the STK parameter list reads
`projects/EXA/structure/linear/Plane_Strain/4elements.e`.

**Physics (FEM).** One `STRUC_LIN_ISO` constitutive model with
`model_type = fem::Model_Type::PLANE_STRAIN` on `UX,UY`. IWGs:

- `STRUC_LINEAR_BULK` on `block_main`, bedding 10⁻⁶, thickness 0.1;
- two `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` IWGs with zero prescribed displacement and a
  `Select` matrix: `Func_Select_Y` (diag(0, 1)) on `sideset_bottom`, `Func_Select_X`
  (diag(1, 0)) on `sideset_left`; `DIRICHLET_NITSCHE` parameter 100.0;
- `STRUC_LINEAR_NEUMANN` with traction (0.0, −1.0) on `sideset_top`.

IQIs: `DOF` for UX and UY, `VON_MISES_STRESS`, and three `NORMAL_STRESS` IQIs with
`vectorial_field_index` 0, 1 and 2 (x, y and out-of-plane).

**Solver (SOL).** `AMESOS_IMPL` direct solver; Newton with relative residual drop 10⁻⁹, at most 5
iterations; DOFs `UX,UY`.

**Output (VIS).** `Plane_Strain_Problem.exo` on `block_main` with nodal `UX`, `UY`, `STRESSVM`,
`STRESSX`, `STRESSY`, `STRESSOOP`.

## What the test checks

The test requires a return code of 0 and, in serial only, opens `Plane_Strain_Problem.exo`. It
requires 2 dimensions, 9 nodes and 4 elements; node 7 at (1.0, 1.0); the displacement at that
node within a norm of 10⁻⁵ of (0.78, −1.82); and the nodal fields 5, 6 and 7 (σ_x, σ_y, σ_z)
compared with 0.0, −1.0 and −0.3, tolerance 10⁻⁵.

## Run it

The example is registered with `TEST_BASE Plane_Strain_Problem` and `NO_PROCS_SUFFIX`, so the
ctest name is exactly `Plane_Strain_Problem`. It is not currently built: no `add_subdirectory`
in `projects/EXA/structure/linear/CMakeLists.txt` reaches this directory. Once connected, in a
build with `MORIS_USE_EXAMPLES=ON`:

```bash
ctest -R Plane_Strain_Problem --output-on-failure
```

## Notes

- The deck uses the pre-`Module_Parameter_Lists` API (`aParameterLists( i ).push_back( ... )`);
  the leaf `CMakeLists.txt` names this as the reason it is disconnected.
- The stress comparisons are written as `(value − reference) < tolerance`, without an absolute
  value.
- The `.dox` table lists the mesh as `cylinder_geometry.e` and the top-edge pressure as
  (0.0, 1.0); the deck reads `4elements.e` and applies (0.0, −1.0).
- See {doc}`Axisymmetric <Axisymmetric>` for the axisymmetric model type on a body-fitted mesh.
