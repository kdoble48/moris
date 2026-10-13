# Heat method normalization

A two-field diffusion problem on a domain cut by a circular (spherical) level set, used to test
the sensitivities of the heat-method normalization of a design level-set field. A first field
`THETA` is diffused from the interface; a second field `PHID` is solved with `THETA` as a
coupled input. Error measures compare `PHID` with the absolute value of the design level set.
A `SWEEP` compares the adjoint gradients of all criteria with forward, backward, and central
finite differences. The .dox lists the feature as "Heat Method Normalization".

```{exa-summary} projects/EXA/optimization/Heat_Method
```

## Problem

| | |
|---|---|
| Domain | [−1, 1]² (2D) or [−1, 1]³ (3D), selected by `gDim` |
| Geometry | user-defined signed distance `Inclusion`: 0.471 − ‖x‖ (circle or sphere of radius 0.471 at the origin); values with magnitude below 10⁻⁸ are replaced by 10⁻⁸ |
| Phases | default two-phase split by the single geometry: `p1` ("inner") and `p0` ("outer"), both active |
| `THETA` | linear isotropic diffusion: conductivity 1.0, density 1.0, heat capacity 0.1; `THETA` = 1.0 on both sides of the interface (Nitsche) |
| `PHID` | linear isotropic diffusion: conductivity 1.0, density 0.0, heat capacity 0.0, `EigenStrain` property 1.0, `THETA` as a dependency; `PHID` = 0.0 on both sides of the interface (Nitsche) |
| Time | one time step, time-continuity term on `THETA` with weights 10.0 and initial condition 0.0 |
| Criteria | H1 errors of `PHID` against a constant and against the design level set, inner and outer volumes, maxima of `THETA` and `PHID` |

## Setup in MORIS

**Mesh (HMR).** 4 elements per direction on a domain of edge 2.0, offset to −1.0, no initial
refinement. One Lagrange and one B-spline mesh, both of order `gInterpolationOrder` (1 or 2,
set by the test).

**Geometry and phases (GEN).** One `USER_DEFINED` level set (`field_function_name =
"Inclusion"`, name `ADVfield`) discretized on B-spline mesh 0 with bounds ±2.0; its
coefficients are the ADVs. A `SCALED_FIELD` property (`LevelsetField`, factor 1.0) publishes
it as PDV `LS1` on all four bulk sets
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
`IQI_types` lists six criteria in this order: `IQIH1ErrorConst`, `IQIH1Error`,
`IQIBulkVolume_Inner`, `IQIBulkVolume_Outer`, `IQIThetaMeasure`, `IQIPhidMeasure`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, Exodus output of the integration mesh.

**Physics (FEM).** Two `DIFF_LIN_ISO` models, `CMDiffusionTheta` and `CMDiffusionPhi`. Each
field has a `SPATIALDIFF_BULK` IWG over the whole domain,
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` IWGs on `iside_b0_1_b1_0` and `iside_b0_0_b1_1`
(`DIRICHLET_NITSCHE` penalty 100.0), and `GHOST_NORMAL_FIELD` IWGs on `ghost_p0` and
`ghost_p1` (`GHOST_DISPL` 0.01). `THETA` also has a `TIME_CONTINUITY_DOF` IWG. The design
level set enters through FEM properties with user value and DV-derivative functions:
`PropLevelSet` returns |φ| and `PropLevelSetGradx` its spatial gradient (sign-corrected,
`tUseAbsoluteValue = true`). The criteria are:

| IQI | Type | Definition in the deck |
|---|---|---|
| `IQIH1ErrorConst` | `H1_ERROR` on `PHID` | references `PropLevelSetConst` (1.0) and `PropLevelSetGradxConst` (1.0 per component); `function_parameters "1.0 / 1.0 / 0.0"` |
| `IQIH1Error` | `H1_ERROR` on `PHID` | references `PropLevelSet` and `PropLevelSetGradx`; same parameters |
| `IQIBulkVolume_Inner` / `_Outer` | `VOLUME` | `p1` / `p0` bulk sets |
| `IQIThetaMeasure` / `IQIPhidMeasure` | `MAX_DOF` on `THETA` / `PHID` | `function_parameters "1.0/1.0"` |

FEM-internal finite differences use `POINT_3_CENTRAL` with a perturbation of 10⁻⁷.

**Solver (SOL).** Direct `AMESOS_IMPL` solve, no preconditioner. Two Newton sub-solvers (one
for `THETA` with `PHID` secondary, one for `PHID` with `THETA` secondary) are combined in an
`NLBGS_SOLVER` (nonlinear solver 2), which the time solver uses
(`TSA_Nonlinear_Solver = 2`). Each algorithm allows 2 iterations to a relative residual drop of
10⁻⁹. A fourth nonlinear solver (Newton on `THETA;PHID`, algorithm 2) is also defined. Both
fields are initialized to 0.0.

**Optimization (OPT).** `SWEEP` with `finite_difference_type = "all"`,
`finite_difference_epsilons = 1.0e-7`, one evaluation per ADV, and
`finite_difference_adv_indices` = `7,8,11,12,13` (2D) or `37,61,62,67,93` (3D). The objective
is the sum of the first four criteria; all six criteria are returned as constraints. Results
go to `SEN_HeatMethod_<case>.hdf5`.

**Output (VIS).** `HeatMethod_<case>.exo` on the bulk and interface sets with nodal `THETA`,
`PHID`, and `PHIdesign` (the `PROPERTY` IQI of `PropLevelSet`).

## What the test checks

Two test cases run four configurations on one processor:

| Test case | Case index | Dimension | Order | Nodes | Elements | `THETA` | `PHID` | `PHIdesign` |
|---|---|---|---|---|---|---|---|---|
| `HeatMethod_Linear` | 0 | 2D | 1 | 82 | 60 | 5.917299018e-01 | 3.872167974e-01 | 4.71e-01 |
| `HeatMethod_Linear` | 1 | 3D | 1 | 414 | 632 | 7.194331940e-01 | 3.587919624e-01 | 4.71e-01 |
| `HeatMethod_Quadratic` | 2 | 2D | 2 | 130 | 60 | 6.349840630e-01 | 2.740283249e-01 | 2.862624742e-01 |
| `HeatMethod_Quadratic` | 3 | 3D | 2 | 1088 | 920 | 7.508299325e-01 | 2.459025775e-01 | 2.506857165e-01 |

For each case the test requires the Exodus dimension, node count, and element count to match
exactly; the reference node (ID 4 in 2D, 13 in 3D) to lie at the origin (relative tolerance
10⁻⁸); the time value to equal 1.0 (10⁻⁸); and `THETA`, `PHID`, and `PHIdesign` at the
reference node to match the table within a relative error of 10⁻⁴. It then reads the analytical
objective and constraint gradients and the `fd_forward`, `fd_backward`, and `fd_central`
gradients from the sweep file and checks each entry to a relative error below 10⁻³.

## Run it

```bash
ctest -R HeatMethod --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is `HeatMethod-1-procs`. The deck reads the globals `gInterpolationOrder`,
`gDim`, and `gTestCaseIndex`, which the test sets for each case.

## Notes

- The Catch2 tags for both test cases are `[structure]`, not `[optimization]`.
- The constants `tNLA_rel_res_norm_drop`, `tNLA_max_iter`, `tTSA_Num_Time_Steps`, and
  `tTSA_Time_Frame` are declared but the SOL parameter list uses literal values.
- In cases 0 and 1 the reference `PHIdesign` at the origin, 0.471, equals the inclusion radius.
