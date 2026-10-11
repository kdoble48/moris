# Cluster-measure sensitivity sweep (two-material diffusion)

A finite-difference check of design sensitivities for criteria built from cluster measures. A
single straight interface splits a unit square (2D) or unit cube (3D) into two diffusive
materials; the four parameters of the interface line are the ADVs. The criteria include volumes,
the interface length, two "strain energies" of the diffusion constitutive models, and four
`STABILIZATION` IQIs that integrate `MEASURE` stabilization parameters (cell, cell-side and
cell-length measures). The SWEEP algorithm compares their analytical gradients against forward,
backward and central finite differences.

```{exa-summary} projects/EXA/optimization/SA_Cluster_Measure
```

## Problem

| | |
|---|---|
| Domain | unit square [0, 1]² (2D) or unit cube [0, 1]³ (3D), one background element per direction |
| Interface | line through (1.35, 0) with normal (1, 1), i.e. x + y = 1.35; it cuts the boundary at (1.0, 0.35) and (0.35, 1.0) |
| Phases | phase 1 on the negative side of the line (x + y < 1.35), phase 0 on the positive side |
| Material (both phases) | conductivity 1.0, density 0.0, heat capacity 0.0, body heat load 1.0 |
| Boundary conditions | side ordinal 4 (x = 0): prescribed temperature 2.0 (Nitsche, both phases) and a surface flux of 5.0 (phase 1) |
| Interface condition | temperature continuity across the material interface, symmetric Nitsche |
| Analysis | steady, one time step, time frame 1.0 |
| Criteria | 9: objective = interface length, constraints = the remaining 8 criteria |

The 3D case uses the same `LINE` geometry (only x and y enter its field), so the interface is a
plane extruded in z.

## Setup in MORIS

**Mesh (HMR).** `number_of_elements_per_dimension` is one element per direction, with
`domain_dimensions` 1.0 in each direction; both are sized from the global `gDim`. Lagrange and
B-spline orders follow `gInterpolationOrder` (set to 1 by the test). Refinement and staircase
buffers are 1.

**Geometry and phases (GEN).** One `LINE` level-set geometry. `center_x`, `center_y`,
`normal_x` and `normal_y` are each given as (lower bound, initial value, upper bound) with all
three equal (1.35, 0.0, 1.0, 1.0), which makes them four ADVs that the sweep perturbs. Phases come
from a user function (`get_phase_index_2d` or `get_phase_index_3d`, chosen by `gDim`) instead of
a phase table. `IQI_types` lists nine criteria in this order: `IQIVolumeInterface1`,
`IQIVolumeMat0`, `IQIVolumeMat1`, `IQISidesetMeasure`, `IQISidesetLengthMeasure`,
`IQIBulkMeasure`, `IQIBulkLengthMeasure`, `IQIBulkStrainEnergy0`, `IQIBulkStrainEnergy1`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on (`tUseGhost = true`), `high_to_low_dbl_side_sets`, and Exodus output of the integration mesh.

**Physics (FEM).** The deck uses the phase-based FEM input (`FEM::PHASES`: `PhaseMaterial0`,
`PhaseMaterial1`, `PhaseAll`) instead of explicit mesh-set names. Per phase: a `DIFF_LIN_ISO`
constitutive model and a `SPATIALDIFF_BULK` IWG. The interface uses
`SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` with a `NITSCHE_INTERFACE` parameter (100.0); the
boundary temperature uses `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` with `DIRICHLET_NITSCHE`
(100.0) on each phase; the flux uses `SPATIALDIFF_NEUMANN`. Ghost penalty: `GHOST_NORMAL_FIELD`
IWGs with `GHOST_DISPL` parameters (0.001) on both phases. Four `MEASURE` stabilization
parameters on phase 1 carry the cluster measures `CELL_SIDE_MEASURE`, `CELL_LENGTH_MEASURE`
(side set and bulk) and `CELL_MEASURE`, all as `PRIMARY,LEADER`; each is read out by a
`STABILIZATION` IQI. The remaining IQIs are `VOLUME` (per phase, and on the phase-1/phase-0
interface side set), `STRAIN_ENERGY` per phase, and a nodal `DOF` IQI for TEMP. The
`COMPUTATION` list sets `is_analytical_sensitivity = false` with `POINT_3_CENTRAL` differencing,
perturbation 10⁻⁶ and the `ABSOLUTE` perturbation strategy.

**Solver (SOL).** `AMESOS_IMPL` direct solve, no preconditioner; Newton with relative residual
drop 10⁻⁸, relaxation 1.0, at most 2 iterations; one time step.

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1`, `include_bounds = false`,
`finite_difference_type = "all"` and `finite_difference_epsilons = 1e-6`; objective and
constraint gradients are written to `Cluster_Measure_2Mat_SA<dim>.hdf5`. The objective is
criterion 0; constraints 0–7 are criteria 1–8, all with constraint type 1.

**Output (VIS).** `Cluster_Measure_2Mat_SA<dim>.exo` on both phases: nodal TEMP and the global
fields `VOLUME_MAT1`, `VOLUME_MAT0`, `VOLUME_INT1`, written every step.

## What the test checks

One test case runs the deck twice, first with `gDim = 2` (test case 0), then with `gDim = 3`
(test case 1), both with interpolation order 1. For each run it reads the Exodus and sweep files
and requires:

| Quantity | 2D | 3D | Tolerance |
|---|---|---|---|
| dimensions / nodes / elements | 2 / 16 / 8 | 3 / 58 / 68 | exact |
| coordinates of node 5 (2D) or 20 (3D) | (1.0, 0.35) | (1.0, 0.35, 1.0) | relative 10⁻⁶ |
| time value | 1.0 | 1.0 | relative 10⁻⁸ |
| temperature at that node (3rd nodal field) | 2.641125244043063 | 2.655482875127814 | relative 10⁻⁴ |

It then compares the analytical objective and constraint gradients with the forward, backward
and central FD gradients for every ADV and every constraint, requiring a relative difference
below 10⁻³ (`REQUIRE`).

## Run it

```bash
ctest -R SA_Cluster_Measure-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The deck reads the globals `gDim` and `gInterpolationOrder`, which the test case sets; there are
no environment variables.

## Notes

- The deck basename (`Cluster_Measure_2Mat_SA`) differs from the directory and ctest name.
- `IQIBulkStrainEnergy0` and `IQIBulkStrainEnergy1` apply `STRAIN_ENERGY` to the diffusion
  constitutive models (`CMDiffusion0`, `CMDiffusion1`).
- The Neumann flux and both Dirichlet conditions are assigned to the same side ordinal (4).
