# Heat conduction around a circular inclusion (2D)

A steady 2D heat-conduction problem on the unit square with a circular inclusion cut by XTK.
Both phases carry a volumetric heat source, the inclusion surface carries an additional surface
heat flux, and the left edge is held at zero temperature. The deck demonstrates a two-phase
thermal model with an immersed material interface coupled by Nitsche's method, ghost
stabilization in both phases, and a level-set field stored on a quadratic B-spline mesh while
the temperature is interpolated on a refined linear B-spline mesh. The same physics on an STK
input mesh is {doc}`STK_XTK_HeatConduction <STK_XTK_HeatConduction>`.

```{exa-summary} projects/EXA/thermal/diffusion/HeatConduction_2D
```

## Problem

| | |
|---|---|
| Domain | 1.0 × 1.0 square |
| Phases | 0 "solid" (outside the circle), 1 "void" (inside the circle) |
| Inclusion | circle centred at (0.5, 0.5), radius 0.31 |
| Material, phase 0 | density 1.0, heat capacity 1.0, conductivity 0.1916229656078740 |
| Material, phase 1 | density 1.0, heat capacity 1.0, conductivity 10.0 |
| Volumetric heat source | 100.0 in both phases (`Func_HeatLoad`) |
| Surface heat flux | 100.0 on the phase-0 side of the circle (Neumann) |
| Dirichlet condition | T = 0 on the left edge (side ordinal 4), Nitsche-enforced |
| Time | static: `tTimeSteps = 0`, so the time-continuity IWG and the time-step settings are not added |

Despite its name, phase 1 is a conducting material with a higher conductivity than phase 0.

## Setup in MORIS

**Mesh (HMR).** 4 × 4 elements on the unit square. Two refinement patterns: pattern 0 (no
initial refinement) carries a quadratic B-spline mesh for the level set; pattern 1 (two levels
of initial refinement, element edge 1/16) carries a linear B-spline mesh for the temperature.
The single Lagrange mesh is quadratic on pattern 1 (`lagrange_orders` = 2,
`bspline_orders` = "2,1", `lagrange_to_bspline` = "0,1"). MSI assigns `TEMP` to B-spline mesh 1
([coupling C6](../../../user-guide/interlinks.md#c6-mesh-indices-hmr-xtk-gen)).

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Func_Sphere`
(φ = 0.31 − ‖x − (0.5, 0.5)‖), discretized on B-spline mesh 0 (`discretization_mesh_index` = 0,
bounds ±2.0). `get_phase_index` maps the positive side of the level set to phase 1. GEN also
lists seven `IQI_types` for the user-defined OPT callbacks
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on, `high_to_low_dbl_side_sets`, and Exodus output of the integration mesh.

**Physics (FEM).** The deck uses the phase-based FEM input (`PhaseSolid`, `PhaseVoid`,
`PhaseAll`) instead of explicit mesh-set names. Each phase has a `DIFF_LIN_ISO` constitutive
model and a `SPATIALDIFF_BULK` IWG with the volumetric load. The solid–void interface is coupled
by `SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, penalty 100.0). The left-edge
temperature uses `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, penalty 100.0);
the surface flux uses `SPATIALDIFF_NEUMANN` on the phase-0 side set adjacent to phase 1
(`neighbor_phases` = `PhaseVoid`). Ghost penalty: `GHOST_NORMAL_FIELD` with `GHOST_DISPL`
(0.05) in each phase. IQIs: nodal temperature (`DOF`), volumetric and surface input energy
(`PROPERTY`), diffusive flux `THERMAL_ENERGY_DIFFUSIVE_FLUX` on each of the four outer edges, and
stored thermal energy (`PROPERTY` of ρ c_p ∂T/∂t).

**Solver (SOL).** Direct solve with `Amesos_Mumps` if MORIS is built with MUMPS, else
`Amesos_Superludist`; Newton with at most 20 iterations, relative residual drop 10⁻¹², no
relaxation. The Jacobian is written to MATLAB format (`SOL_save_operator_to_matlab`).

**Optimization (OPT).** `is_optimization_problem` is false, so the SWEEP algorithm block and the
OPT callbacks are not exercised; the deck runs one forward solve.

**Output (VIS).** `HeatConduction.exo` with nodal `TEMP` and `HEAT` on the bulk sets of both
phases, time offset 10.0.

## What the test checks

The Catch2 case `HeatConduction` calls the workflow once and requires a return code of 0, then
calls it a second time without checking the return code, and compares
`HeatConduction.exo.e-s.0000`:

- 2 dimensions, 1425 nodes, 556 elements (exact);
- node 641 at (0.6875, 0.84375), relative tolerance 10⁻⁵;
- time value 1.0, relative tolerance 10⁻⁸;
- nodal temperature 407.93805018193 at node 641, relative tolerance 10⁻⁵.

The test source also holds a 4-processor reference (node 107, temperature 384.267960788395),
but CMake registers only the 1-processor run.

## Run it

```bash
ctest -R HeatConduction-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The `.dox` documentation file is absent from this directory.
- The analytic-geometry branch (`tUseAnalyticGeometry = true`, radius as a design variable) is
  present in the deck but switched off.
