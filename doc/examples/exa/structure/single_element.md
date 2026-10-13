# Single-element transient thermoelastic problem

One plane-stress element is heated through its top edge and pressed on its left edge while its
bottom edge is held fixed. Temperature is transient (three time steps) and enters the
structural problem through thermal expansion. The deck solves the coupled problem with a
staggered nonlinear block Gauss–Seidel (NLBGS) solver over a temperature and a displacement
Newton sub-solver, uses GGLS stabilization on the diffusion equation, and assigns one time level
to the displacements and two to the temperature.

```{exa-summary} projects/EXA/structure/thermo_elastic/single_element
```

## Problem

| | |
|---|---|
| Domain | 0.01 × 0.05 rectangle, one element |
| Phases | one (`HMR_dummy_n_p0`); the level set x − 100 does not cut the domain |
| Thermal | density 8960, heat capacity 1460, conductivity 370; initial temperature 1.0 |
| Structural | E = 120e9, ν = 0.34, plane stress; CTE 0.1, reference temperature 0.0; bedding 10⁻⁵ is defined |
| Heat load | flux 1.0e4 on the top edge (`SideSet_3_n_p0`) |
| Mechanical load | pressure 100000.0 on the left edge (`SideSet_4_n_p0`) |
| Support | displacement (0, 0) on the bottom edge (`SideSet_1_n_p0`), Nitsche |
| Time | 3 time steps over a time frame of 0.3 |

## Setup in MORIS

**Mesh (HMR).** 1 × 1 element on 0.01 × 0.05; Lagrange and B-spline order `tOrder = "1"`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Back_Wall` (x − 100, clamped to
10⁻⁸ near zero); the GEN output mesh is written to `GEN_single_element.exo`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, no ghost
stabilization, enriched and Exodus integration-mesh output on.

**Physics (FEM).** The deck is organized by the flags `tHaveThermo`, `tHaveStruct`,
`tHaveStaggered`, `tHaveMixedSolverConfig`, `tUseGGLS` (all `true`) and `tUseBelos` (`false`).
With these settings:

- thermal: `DIFF_LIN_ISO` (conductivity, density, heat capacity); `SPATIALDIFF_BULK` with the
  `GGLS_DIFFUSION` stabilization parameter (fed dummy phase-change properties: latent heat 0.0,
  phase-state function 2.0, phase-change constant 100.0); `SPATIALDIFF_NEUMANN` on the top edge;
  `TIME_CONTINUITY_DOF` with weights 100.0 / 100.0 and the initial condition 1.0
  (`time_continuity = true`);
- structural: `STRUC_LIN_ISO` (`PLANE_STRESS`) with dependencies `UX,UY;TEMP` and properties
  `YoungsModulus`, `PoissonRatio`, `CTE`, `ReferenceTemperature`; `STRUC_LINEAR_BULK`;
  `STRUC_LINEAR_NEUMANN` with the `Pressure` property on the left edge;
  `STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` on the bottom edge (`DIRICHLET_NITSCHE` 100.0);
- IQIs: `DOF` for TEMP, UX and UY, and `VOLUME` over the bulk.

**Solver (SOL).** `AMESOS_IMPL`. Nonlinear algorithm 0 is Newton (relative drop 10⁻¹⁰, relaxation
1.0, at most 3 iterations, combined residual/Jacobian assembly); algorithm 1 is NLBGS. Nonlinear
solver 0 (Newton) acts on `UX,UY`, solver 1 (Newton) on `TEMP`, and solver 2 (NLBGS, algorithm 1)
on `UX,UY;TEMP` with `NLA_Sub_Nonlinear_Solver = "1,0"` (the temperature sub-solver listed first). The time solver
uses nonlinear solver 2, `TSA_Initialize_Sol_Vec = "UX,0.0;UY,0.0;TEMP,1.0"` and
`TSA_time_level_per_type = "UX,1;UY,1;TEMP,2"`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). The Jacobian is
written to `jacobian.dat` (`SOL_save_operator_to_matlab`).

**Output (VIS).** `single_element.exo` on the bulk: nodal TEMP, UX, UY and global VOLUME, written
every step.

## What the test checks

The test runs `single_element.so`, requires a return code of 0, and in serial only (other
processor counts raise an error) opens `single_element.exo`. It requires 2 dimensions, 4 nodes and
1 element, then at node 3:

| Quantity | Reference | Relative tolerance |
|---|---|---|
| coordinates | (0.0, 0.05) | 10⁻⁸ |
| time value | 0.10 | 10⁻⁸ |
| TEMP (time step index 2) | 1.334186887938226 | 10⁻⁸ |
| UX (time step index 2) | −8.310495618961509e−04 | 10⁻⁶ |
| UY (time step index 2) | 5.243129489358523e−03 | 10⁻⁶ |
| global variable 0 (VOLUME) | 5.00e−05 | 10⁻⁸ |

## Run it

The example is registered with `NO_PROCS_SUFFIX`, so the ctest name is exactly `single_element`:

```bash
ctest -R '^single_element$' --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The SOL list contains complete configurations for thermal-only, structural-only and monolithic
  coupled runs; editing the flags at the top of the deck switches between them.
- `check_quadratic_results_serial` (expecting 9 nodes) is defined in the test but never called;
  only the linear check runs.
- `PropTraction` (`Func_Neumann_U`) and `PropBedding` are defined but not referenced by any IWG.
