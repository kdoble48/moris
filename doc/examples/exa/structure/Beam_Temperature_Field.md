# Beam deflection under an imported temperature field

A linear elastic cantilever beam is loaded only by thermal expansion: a nodal temperature
field, read from an Exodus file, varies linearly from 0 at the bottom of the beam to 1 at the
top and bends the beam through the coefficient of thermal expansion (CTE). The example tests the
CTE path of the linear isotropic constitutive model on a body-fitted STK mesh with quadratic
elements, and the import of a nodal field through the FEM `FIELDS` parameter list.

```{exa-summary} projects/EXA/structure/linear/Beam_Temperature_Field
```

## Problem

| | |
|---|---|
| Domain | 2D beam, 20 × 1 (inches, per the `.dox`), body-fitted mesh `beam_geometry_quadratic.e` |
| Phases | one block, `block_main` |
| Material | linear elastic, isotropic: E = 1.0, ν = 0.3, CTE = 0.01, reference temperature 0.0 |
| Temperature | nodal field from `beam_temperature_quadratic.e`, 0 at the bottom and 1 at the top |
| Supports | `sideset_left`, displacement (0, 0) imposed weakly by Nitsche |
| Load | thermal strain only; no traction or body force |
| Analysis | static, one solve |

The `.dox` page derives the beam-theory reference: the curvature is α(T_top − T_bottom)/d, so the
tip deflection is v(L) = α(T_top − T_bottom)L²/(2d). The MATLAB script
`DeflectionDueToTemperature.m` evaluates this with α = 0.01, L = 20, d = 1, ΔT = 1, giving a
magnitude of 2.0; `BeamTemperatureField.cae` holds a companion Abaqus model.

## Setup in MORIS

**Mesh (STK).** The workflow is `STK_FEM` (no HMR, GEN or XTK). The mesh is read from
`projects/EXA/structure/linear/Beam_Temperature_Field/beam_geometry_quadratic.e` through the
`input_file` key of the STK parameter list; the discretization is the second-order one stored in
the file.

**Physics (FEM).** One `STRUC_LIN_ISO` constitutive model (`CMStrucLinIso1`) on `UX,UY` takes the
properties `YoungsModulus`, `PoissonRatio`, `CTE`, `ReferenceTemperature` and
`PropertyTemperature`. The temperature property `PropField` is evaluated by `Func_Field`, which
returns the value of the field interpolator for `FIELD_1`. That field is declared in
`FEM::FIELDS` as `FieldNodalTEMP` (`field_entity_type = NODAL`,
`field_create_from_file = beam_temperature_quadratic.e`). The residual is a
`STRUC_LINEAR_BULK` IWG on `block_main` with a bedding stiffness of 10⁻⁸ and a
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` IWG on `sideset_left` with a `DIRICHLET_NITSCHE`
stabilization parameter of 100.0. IQIs: `DOF` for UX and UY, and `PROPERTY` for the imported
temperature (`IQITempField`).

**Solver (SOL).** `AMESOS_IMPL` direct solver, Newton with relative residual drop 10⁻⁹, at most
5 iterations, no preconditioner. The DOF set is `UX,UY`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Output (VIS).** `Beam_Temperature_Field_Problem.exo` on `block_main` with nodal `UX`, `UY` and
`Temp`, written every step.

## What the test checks

The test runs the deck through the workflow interface, requires a return code of 0, and then
opens `Beam_Temperature_Field_Problem.exo` (serial only; any other processor count raises an
error). It requires:

- 2 dimensions, 3975 nodes and 924 elements;
- node 3973 at coordinates (10.0, 0.5), relative error below 10⁻⁶;
- UY at that node, and UY interpolated at those coordinates, compared with −2.0 with tolerance
  10⁻³;
- the temperature at that node and at those coordinates compared with 1.0 with tolerance 10⁻⁶.

## Run it

```bash
ctest -R Beam_Temperature_Field-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The deck locates its mesh and field files relative to the MORIS source root returned by
`get_base_moris_dir()`.

## Notes

- The displacement and temperature comparisons are written as `(value − reference) < tolerance`,
  without an absolute value, so they bound the difference from one side only.
- The `.dox` feature table lists a "simple" Dirichlet condition and a nodal Neumann condition;
  the deck contains only the Nitsche condition on `sideset_left` and no Neumann IWG.
- The Catch2 case is named `Beam_Temperature_Field_write`.
