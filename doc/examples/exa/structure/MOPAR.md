# Axisymmetric thermoelastic analysis of a two-material aeroshell section

A two-material aeroshell section — an aluminum core with an ablator layer on its exterior,
according to the `.dox` page — is analyzed as an axisymmetric linear elastic body. The loads are a
temperature field imported from an Exodus file (acting through thermal expansion) and a wall
traction interpolated from a tabulated data file. The deck runs the `STK_FEM` workflow on a
body-fitted triangle mesh and demonstrates the `AXISYMMETRIC` model type, temperature-dependent
material properties driven by an imported nodal field, and von Mises stress recovery through a
staggered (NLBGS) solve.

```{exa-summary} projects/EXA/structure/thermo_elastic/MOPAR
```

## Problem

| | |
|---|---|
| Mesh | `exomesh.e`, 2D, blocks `block_1` and `block_2`, side sets `side_1`–`side_3` |
| Rotation axis | through (0, 0) and (1, 0), i.e. the x-axis |
| `block_1` | E = 50 GPa, ν = 0.10, CTE = 2.0e−6 (deck variables `*_c`) |
| `block_2` | E = 6.89 GPa, ν = 0.33, CTE = 23.6e−6 (deck variables `*_al`) |
| Reference temperature | 293.0 |
| Temperature | nodal field `FIELD_1` read from `traction.e` |
| Traction | on `side_1`, linearly interpolated in y from columns of `wall_data.txt` |
| Supports | UX = 0 on `side_2` and `side_3` (Nitsche with `Select` = diag(1, 0)) |
| Analysis | static, single time step |

The Young's modulus is evaluated by `Func_Thermal_Modulus_Block`, which interpolates in the
temperature table (0, 293, 1000, 10000); both moduli tables hold the same value at every
temperature, so E is constant in each block. The densities (1879 and 2700) are defined but enter
no IWG.

## Setup in MORIS

**Mesh (STK).** Workflow `STK_FEM`; the mesh is
`projects/EXA/structure/thermo_elastic/MOPAR/exomesh.e`. The `.dox` page and `readme.txt` record
how it was built (Tecplot data → ParaView → Cubit, `tri3` elements); `mapper.m` and the `misc/`
MATLAB scripts belong to that preprocessing chain.

**Physics (FEM).** Two `STRUC_LIN_ISO` constitutive models with
`model_type = fem::Model_Type::AXISYMMETRIC`, each taking `YoungsModulus`, `PoissonRatio`,
`AxisymRotationAxis` (`Func_AxiRotation`), `CTE`, `PropertyTemperature` (`Func_Field`, the
imported field) and `ReferenceTemperature`. With `byGeometry = false` the properties are chosen
per block; the alternative path selects them by the line y = 0.0160091/0.09 (x − 0.01). IWGs, each
with the rotation-axis property passed as `Thickness`:

- `STRUC_LINEAR_BULK` on `block_1` and `block_2`;
- `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on `side_2` and `side_3` (`DIRICHLET_NITSCHE`
  100.0);
- `STRUC_LINEAR_NEUMANN` on `side_1` with traction `Func_Import_Traction`;
- `STRUC_VON_MISES_STRESS` projecting the von Mises stress of `block_1` into `STRESS_DOF` and that
  of `block_2` into `VX`; on the opposite block each of these fields has an `L2` IWG with zero
  source, zero L2 coefficient and diffusion 10⁻⁵.

IQIs: `DOF` for UX, UY, `STRESS_DOF` (`block_1`) and `VX` (`block_2`); `PROPERTY` for the
temperature and the two moduli. The temperature field is declared in `FEM::FIELDS`
(`FieldNodalTEMP`, `NODAL`, `field_create_from_file = traction.e`).

**Solver (SOL).** `AMESOS_IMPL`. Three Newton sub-solvers on `UX,UY`, `STRESS_DOF` and `VX`
(relative drop 10⁻⁹, 1 iteration) are chained by an `NLBGS_SOLVER` (algorithm 1, one iteration)
over `UX,UY;STRESS_DOF;VX`; the time solver uses nonlinear solver 3.

**Output (VIS).** `AxisymmetricProblem.exo` on both blocks: nodal UX, UY, `STRESS1`, `STRESS2`,
`Temp` and elemental-average `EMOD1`, `EMOD2`.

## What the test checks

The test calls the workflow interface (return code not checked) and, in serial only, opens
`AxisymmetricProblem.exo` and requires 2 dimensions, 7059 nodes and 13228 elements (exact), then
at node 729, with relative tolerance 10⁻⁵:

| Quantity | Reference |
|---|---|
| coordinates | (1.882561482489109e−02, 2.055238932371140e−03) |
| time value | 1 (tolerance 10⁻⁸) |
| displacement (UX, UY) | (−1.128321042350804e−04, 9.853289535187406e−06) |
| temperature (nodal field 6) | 1.038658447265625e+03 |
| stress (nodal field 4, `STRESS1`) | 2.397992760916115e+08 |

## Run it

```bash
ctest -R AxisymmetricProblem-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The deck reads `exomesh.e`, `traction.e` and `wall_data.txt` from the MORIS source tree via
`get_base_moris_dir()`.

## Notes

- The `.dox` table states a Dirichlet condition "V = 0.5 at inlet" and describes side set 1 as a
  heat-flux Dirichlet condition; in the deck `side_1` carries the structural traction and the
  Dirichlet conditions are UX = 0 on `side_2` and `side_3`. No heat equation is solved.
- By the `.dox` property tables, the `block_1` values (ν = 0.10, CTE = 2.0e−6) are the ablator's
  and the `block_2` values (ν = 0.33, CTE = 23.6e−6) are aluminum's; the deck itself labels them
  only by the suffixes `_c` and `_al`.
- For another axisymmetric example see {doc}`Axisymmetric <Axisymmetric>`.
