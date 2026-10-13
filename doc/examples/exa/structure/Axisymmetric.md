# Axisymmetric cylinder under radial pressure

A linear elastic, static analysis of a solid cylinder loaded by a radial pressure, modeled as a
2D axisymmetric problem on a body-fitted mesh read from an Exodus file. The example tests the
axisymmetric model type of the linear isotropic constitutive model and the axisymmetric
weighting of the bulk, Dirichlet and Neumann IWGs. A second deck repeats the analysis with
thermal-expansion properties added to the constitutive model.

```{exa-summary} projects/EXA/structure/linear/Axisymmetric
```

## Problem

| | |
|---|---|
| Domain | square (x, y) ∈ [−1, 1] × [−1, 1] (size 2.0 × 2.0), from `cylinder_geometry.e` |
| Mesh | 2 × 2 first-order quadrilaterals, 9 nodes, one block `block_main` |
| Axis of rotation | the bottom edge, the line through (0, −1) and (1, −1) |
| Material | linear elastic, isotropic, E = 1.0, ν = 0.3 (density property 1.0 is defined but not used by any model) |
| Support | `sideset_left`: Nitsche-enforced Dirichlet condition with prescribed value (0, 0) on the x component only (selection matrix diag(1, 0)) |
| Load | `sideset_top`: traction (0.0, −1.0) |
| Analysis | static, linear |

The radial coordinate is the distance from the rotation axis, so y is the radial direction and x
the axial direction. The `example_documentation.dox` file states the reference solution: for this
solid cylinder under radial pressure only, the radial stress equals the applied traction, and the
remaining stress components follow from the constitutive law with a radial displacement linear in
*r*. `AxisymmetricCylinder.cae` holds the reference model.

## Setup in MORIS

**Mesh (STK).** The workflow is `STK_FEM` (set in the OPT list with
`is_optimization_problem = false`); the STK list reads `cylinder_geometry.e` from the MORIS
source tree via `get_base_moris_dir()`. There is no HMR, XTK or GEN block, and the FEM sets are
the block and side-set names stored in the Exodus file (`block_main`, `sideset_left`,
`sideset_top`; the file also carries `sideset_bottom` and `sideset_right`).

**Physics (FEM).** One `STRUC_LIN_ISO` constitutive model with
`model_type = fem::Model_Type::AXISYMMETRIC`. Besides Young's modulus and Poisson's ratio it takes
the property `PropRotAxis` in the `AxisymRotationAxis` slot. That property's value function,
`Func_AxiRotation`, receives the axis as two points `{{0.0, −1.0}, {1.0, −1.0}}` and returns the
4-vector {2πr, r, n₁, n₂}: the circumferential length, the radius, and the outward radial unit
normal at the evaluation point. The same property is passed as the `Thickness` property of all
three IWGs, so each residual is integrated with the 2πr weight:

- `STRUC_LINEAR_BULK` on `block_main`, with a bedding property of 10⁻⁶;
- `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on `sideset_left`, with a `DIRICHLET_NITSCHE`
  stabilization parameter (penalty 100.0) and the `Select` property `PropSelectX`;
- `STRUC_LINEAR_NEUMANN` on `sideset_top` with the traction (0.0, −1.0).

IQIs: `DOF` for UX and UY, `VON_MISES_STRESS`, and three `NORMAL_STRESS` IQIs with
`vectorial_field_index` 0, 1 and 2 (axial, radial and out-of-plane/azimuthal stress). The links
follow [coupling C4](../../../user-guide/interlinks.md#c4-inside-fem-bottom-up-references).

**Solver (SOL).** Direct solve with `AMESOS_IMPL` (default backend), Newton with
`NLA_rel_res_norm_drop = 1e-9`, at most 5 iterations, separate residual and Jacobian assembly,
no preconditioner, default time-solver algorithm settings. DOFs `UX,UY`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)); the MSI list is
left empty.

**Output (VIS).** `STANDARD` mesh on `block_main` with nodal fields `UX, UY, STRESSVM, STRESSX,
STRESSRAD, STRESSOOP`, written to `Axisymmetric_Problem.exo`.

## Decks

| Deck | Difference | Output |
|---|---|---|
| `Axisymmetric_Problem.cpp` | as described above | `Axisymmetric_Problem.exo` |
| `Axisymmetric_Problem_ConstantThermalField.cpp` | adds `PropCTE` (0.01), `PropTemp` (100.0, slot `PropertyTemperature`) and `PropRefTemp` (0.0) to the constitutive model; no temperature DOF is solved | `Axisymmetric_Problem_ConstantThermalField.exo` |

## What the test checks

The single Catch2 case `Axisymmetric_Problem_Linear` runs both decks in sequence and requires a
return code of 0 for each. It then aborts unless run on one processor and checks, at node 7:

- mesh of the first deck: 2 dimensions, 9 nodes, 4 elements;
- node coordinates (1.0, 1.0), relative tolerance 10⁻⁵ (both decks);
- first deck only: displacement (UX, UY) = (1.2, −1.4), absolute norm tolerance 10⁻⁵;
- both decks: stress fields 5, 6, 7 (axial, radial, azimuthal) against 0.0, −1.0 and −1.0.

## Run it

```bash
ctest -R Axisymmetric_Problem-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The stress checks are written as `value - reference < 1e-5` without an absolute value, so they
  bound the stresses from above only.
- The deck sets the mesh path from the MORIS source directory, so it runs from any working
  directory as long as `get_base_moris_dir()` resolves.
- The `.dox` file's "Reference Solution" section opens with a sentence about a cantilever beam
  under a temperature gradient, carried over from another example; the figure and remaining text
  describe the cylinder.
