# Leading edge with phase-change material

A 2D transient thermo-mechanical topology optimization of a wedge-shaped leading edge: a metal
shell (skin and internal fins) encloses a phase-change material (PCM), the outer surface
receives an angularly distributed heat flux, radiation and pressure, and the fin layout is
optimized for peak skin temperature, strain energy and fin perimeter under a fin-volume
constraint. The deck couples temperature and displacement through thermal expansion, uses the
phase-change diffusion model `DIFF_LIN_ISO_PC` in the PCM, and solves the linear systems with
Belos (the feature named in the `.dox` page).

```{exa-summary} projects/EXA/optimization/Mach_Leading_Edge
```

## Problem

| | |
|---|---|
| Domain | 0.016 × 0.008 rectangle, offset (−0.004, −0.004) |
| Geometry | wedge of length 12 mm, half-angle 5°, inner radius 1 mm, wall thickness 2 mm; back wall at x = 12 − 1 − 2 mm |
| Phases | 1 back wall, 2 skin, 3 fins, 4 PCM, 0 elsewhere (16-entry phase table) |
| Shell (skin, fins) | density 1.0e−2, heat capacity 200, conductivity 1.5e−4, E = 2.1e5, ν = 0.31, CTE 1.3e−5, plane stress |
| PCM | density 6.0e−3, heat capacity 500, conductivity 1.5e−5, latent heat 220000, phase-change temperature 750, phase-change constant 30, E = 1.0e3, ν = 0.3, CTE 0 |
| Thermal loads | imposed flux 3.75 × `Func_Heat_Load_Distribution` on the outer shell surface; radiation (emissivity 5.0e−7, ambient 0.0, absolute zero −273.15) on the same surface |
| Mechanical loads | pressure 20.265 (200 atm × 1.01325e5 × 10⁻⁶) on the outer shell surface; bedding 2.1 in all bulk phases; displacement (0, 0) on the skin–back-wall interface |
| Initial condition | temperature 700 (`Func_Initial_Condition`); reference temperature for expansion 20.0 |
| Time | transient, 3 time steps over a time frame of 0.2 |
| Objective | 5 · SE / 1.74024e−5 + max-temperature IQI / 0.0458574 + fin perimeter / 0.0288672 |
| Constraint | fin volume / 8.69684e−6 − 1 ≤ 0 |
| Optimizer | GCMMA, step size 0.1, 1 iteration |

Material values are given in the deck's scaled unit system; its comments state the physical
equivalents (for example "10,000 kg/m^3" for the shell density and "375 W/cm^2" for the heat
load). `Func_Heat_Load_Distribution` scales the flux by a function of the polar angle about the
wedge tip, capped at 90° − 5°.

## Setup in MORIS

**Mesh (HMR).** 36 × 18 linear elements, one linear Lagrange mesh and one linear B-spline mesh
on pattern 0, no initial refinement. The MSI parameter list is empty.

**Geometry and phases (GEN).** Three user-defined level sets — `Outer_Wedge`, `Inner_Wedge`,
`Back_Wall` — and a `SUPERELLIPSE` field array for the fins: 3 × 2 superellipses (exponent 6,
semi-diameters 2.0 mm × ≈0.750 mm), centres between x = 1 mm and 12 mm and
y = ±0.9 · 12 mm · tan 5°. The field array is discretized on B-spline mesh 0
(`discretization_mesh_index = 0`) with bounds ±0.0025; these coefficients are the ADVs
([coupling C3](../../../user-guide/interlinks.md#c3-how-advs-reach-the-physics-gen-fem)).
The phase table maps the 16 sign combinations to the five phases. `IQI_types`: `IQIMaxTemp`,
`IQIBSplineGeometryVolume`, `IQIBSplinesPerimeter`, `IQIStrainEnergy`; the objective uses
`aCriteria` 0, 2, 3 and the constraint `aCriteria` 1
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, `high_to_low_dbl_side_sets = true`, `keep_all_opt_iters = true`, and output
of the enriched and the Exodus integration mesh.

**Physics (FEM).** Heat: `DIFF_LIN_ISO` in the shell and `DIFF_LIN_ISO_PC` in the PCM, each
`SPATIALDIFF_BULK` with a `GGLS_DIFFUSION` stabilization; `SPATIALDIFF_NEUMANN` flux and
`SPATIALDIFF_RADIATION` on `iside_b0_2_b1_0`; `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` on the
skin–fin set `dbl_iside_p0_2_p1_3` and the shell–PCM sets `dbl_iside_p0_2_p1_4`,
`dbl_iside_p0_3_p1_4` (`NITSCHE_INTERFACE`, 100); `TIME_CONTINUITY_DOF` (weights 100) on all
bulk sets. Structure: `STRUC_LIN_ISO` with temperature dependence (CTE, reference temperature)
in shell and PCM, `STRUC_LINEAR_BULK` with bedding, `STRUC_LINEAR_NEUMANN` pressure,
`STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE` on the skin–fin interface, and
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 100) on `iside_b0_2_b1_1`.
`GHOST_NORMAL_FIELD` acts on TEMP and UX,UY over `ghost_p2,ghost_p3` and `ghost_p4`
(`GHOST_DISPL`, 0.01). Criteria: `MAX_DOF` on TEMP over the skin
(`function_parameters = "350.0/10.0/0.0"`), `VOLUME` over the fins, `VOLUME` over the
PCM–solid interfaces `iside_b0_4_b1_3,iside_b0_4_b1_2,iside_b0_4_b1_1` as perimeter, and
`STRAIN_ENERGY` over the shell. `VON_MISES_STRESS`, `DOF` and total `VOLUME` IQIs are defined
for output.

**Solver (SOL).** `BELOS_IMPL` with convergence tolerance 10⁻¹² and an `IFPACK` ILU
preconditioner. Nonlinear algorithm 0 is Newton (relative drop 10⁻⁵, up to 10 iterations),
algorithm 1 is NLBGS with default settings, algorithm 2 sets a drop of 10⁻⁷ and 1 iteration.
Nonlinear solver 2 is an `NLBGS_SOLVER` that runs the TEMP Newton solver and then the UX,UY
Newton solver (`NLA_Sub_Nonlinear_Solver = "1,0"`); it is the forward solver of the time
solver. Solver 3, a monolithic solve over `UX,UY,TEMP` on algorithm 2, is the
`TSA_Nonlinear_Sensitivity_Solver`. `TSA_time_level_per_type = "UX,2;UY,2;TEMP,2"`.

**Optimization (OPT).** `user_defined` GCMMA problem with hand-written objective, constraint
and `dcriteria` callbacks (the `dadv` derivatives are zero).

**Output (VIS).** `Mach_Leading_Edge.exo` on all bulk sets: nodal TEMP, UX, UY, STRESS and
global MAX_DOF and two VOLUME fields, `Save_Frequency = 5`, `Time_Offset = 10.0`
([coupling C8](../../../user-guide/interlinks.md#c8-output-vis-fem-and-xtk)).

## What the test checks

The test runs `Mach_Leading_Edge.so` on 1 processor and requires a return code of 0. It does
not compare any output values.

## Run it

```bash
ctest -R Mach_Leading_Edge --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name has no `-procs` suffix (`NO_PROCS_SUFFIX`).

## Notes

- `CMakeLists.txt` notes that the legacy test was also gated on `MORIS_HAVE_GCMMA`;
  `moris_add_example` registers it under `MORIS_HAVE_PARALLEL_TESTS` only.
- The `.dox` page lists "Belos Solver" and a "TODO: Add features here" placeholder.
