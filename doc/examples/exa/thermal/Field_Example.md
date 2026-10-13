# Nodal field write, read, and compare

Three decks that exercise the FEM field I/O: a static heat-conduction problem on an immersed
domain writes its nodal temperature to an HDF5 field file, a second run reads a stored reference
field back as a property and evaluates the difference to the computed temperature, and a third
run does the same while loading a stored solution vector. The example tests the `FEM::FIELDS`
parameter list (`field_output_to_file`, `field_create_from_file`), field-dependent properties,
and `SOL_load_sol_vec_from_file`.

```{exa-summary} projects/EXA/thermal/diffusion/Field_Example
```

## Problem

| | |
|---|---|
| Background domain | 2.0 × 2.0 square, 20 × 20 elements |
| Geometry | level set φ = x² + y²/2 − 0.331; FEM acts on phase 1 only (sets `HMR_dummy_n_p1,HMR_dummy_c_p1`) |
| Material | linear isotropic diffusion, conductivity 1.0 |
| Interface | heat flux 1.0 on `iside_b0_1_b1_0` |
| Right edge | temperature 1.0 on `SideSet_2_n_p1` |
| Analysis | static, one time step (time frame 1.0) |

## Setup in MORIS

**Mesh (HMR).** Quadratic Lagrange mesh on pattern 0 and quadratic B-spline mesh on pattern 1,
refinement and staircase buffers 1.

**Geometry and phases (GEN).** One `USER_DEFINED` level set (`LevelSetFunction`) discretized on
B-spline mesh 0 (`discretization_mesh_index = 0`) with `use_multilinear_interpolation = true`,
no refinement.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on, `high_to_low_dbl_side_sets = true`, integration-mesh Exodus output.

**Physics (FEM).** `DIFF_LIN_ISO` on `HMR_dummy_n_p1,HMR_dummy_c_p1`; IWGs `SPATIALDIFF_BULK`,
`SPATIALDIFF_NEUMANN` on the interface, `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE`
(`DIRICHLET_NITSCHE`, 100.0) on the right edge, and `GHOST_NORMAL_FIELD` with `GHOST_DISPL`
(0.01) on `ghost_p1`. A `DOF` IQI `IQIBulkTEMP` reports TEMP. One `FEM::FIELDS` entry declares
the nodal field `FieldNodalTEMP` of type `FIELD_1`.

**Solver (SOL).** `AMESOS_IMPL`, Newton with relative drop 10⁻⁸, relaxation 1.0, at most 2
iterations, no preconditioner; MSI sets `order_adofs_by_host = false`.

**Output (VIS).** `<deck>.exo` on the phase-1 sets.

### The three decks

| | `Field_example_write` | `Field_example_read` | `Field_example_compare` |
|---|---|---|---|
| `pattern_initial_refinement` | 1, 1 | 1, 1 | 1, 0 |
| `FieldNodalTEMP` | filled from `IQIBulkTEMP`, written to `Field_example_write.hdf5` | read from `Field_example_ref.hdf5` | read from `Field_example_ref.hdf5` |
| Extra property | — | `PropField` = value of `FIELD_1` (`Func_Field`, `field_dependencies = FIELD_1`) | same |
| Extra IQI | — | `IQIBulkL2Error`: `H1_ERROR` on TEMP, `function_parameters` `1.0/0.0`, `PropField` as `L2_Reference` | same |
| Solution vector | solved | solved | loaded from `Solution_Vector_Ref.hdf5` |
| VIS fields | TEMP | TEMP, DIFFERENCE (`IQIBulkL2Error`) | TEMP, DIFFERENCE |

Both reference HDF5 files live in the source directory and are located through
`moris::get_base_moris_dir()`.

## What the test checks

One ctest runs three Catch2 cases in the order write, read, compare. Each checks its Exodus file
for 2 dimensions, 6067 nodes, and 1570 elements, then:

- **write** (`gTestCaseIndex = 0`): loads `FieldNodalTEMP` from `Field_example_ref.hdf5` and
  from `./Field_example_write.hdf5`, checks equal entry counts, and checks
  `value − reference < 10⁻¹²` entry by entry (a signed difference, not an absolute value).
- **read** (`gTestCaseIndex = 1`): node 5992 at (0.175, 0.775) to relative 10⁻¹²; nodal field
  index 3 at step 0 equals 0.0 to 10⁻¹².
- **compare** (`gTestCaseIndex = 2`): same node; nodal field index 3 equals 5.30216 · 10⁻⁶ to
  10⁻¹⁰.

None of the cases checks the workflow return code.

## Run it

```bash
ctest -R Field_Example-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The leaf builds three input libraries (`INPUTS Field_example_write Field_example_read
Field_example_compare`).

## Notes

- The write deck defines `Func_Select_X` and `Func_Select_Y` and a commented-out
  `L2_ERROR_ANALYTIC` IQI; neither is used.
- The compare deck contains a commented-out `SOL_save_final_sol_vec_to_file` line pointing at
  the same `Solution_Vector_Ref.hdf5`.
