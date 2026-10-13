# Channel with four cylinders, static, temperature only

A diffusion-only version of the channel-with-four-cylinders problems
({doc}`Channel_with_Four_Cylinders_Static <Channel_with_Four_Cylinders_Static>`,
{doc}`Channel_with_Four_Cylinders_Transient <Channel_with_Four_Cylinders_Transient>`); the
documentation describes it as an easy setup for debugging the advection problems. A steady
linear heat-conduction problem is solved on an immersed channel with four circular voids, with a
Nitsche-enforced inlet temperature, a heat flux on the cylinder surfaces, and ghost
stabilization. The test solves the linear case with two linear-solver back ends (PETSc when
available, then Amesos) and the quadratic case with Amesos.

```{exa-summary} projects/EXA/thermal/diffusion/Channel_with_Four_Cylinders_Static_Temp_Only
```

## Problem

| | |
|---|---|
| Background domain | 4.0 × 2.0 rectangle, offset (−1.24, −0.86) |
| Physical domain | 2.2 × 0.41 rectangle, 0 ≤ x ≤ 2.2, 0 ≤ y ≤ 0.41 |
| Voids | four circles of radius 0.05 at (0.4, 0.1), (0.8, 0.3), (1.2, 0.1), (1.6, 0.3) |
| Material | density 1.0, heat capacity 1.0, conductivity 0.00005 |
| Left edge (x = 0) | temperature 0.0, Nitsche-enforced |
| Cylinder surfaces | heat flux 2.0 |
| Other boundaries | no condition written (adiabatic per the documentation) |
| Analysis | static, one solve |

The `.dox` page states that no analytical solution exists; mesh size, a nodal temperature, and
the temperature IQI serve as reference values.

## Setup in MORIS

**Mesh (HMR).** 176 × 88 elements, `processor_decomposition_method = 1` with
`processor_dimensions` 2 × 1. Lagrange and B-spline orders come from `gInterpolationOrder`
(1 or 2); the refinement buffer equals the order, the staircase buffer is 1. Refinement comes
from the geometries (one level each), not from an initial pattern.

**Geometry and phases (GEN).** Eight `USER_DEFINED` level sets — four planes and four circles —
each clamped to a minimum magnitude of 10⁻⁶ and refined once on mesh 0. The phase function
`get_phase_index` maps the geometry signs to five phases: 0 material, 1 void inlet, 2 cylinder,
3 void walls, 4 void outlet. GEN output: `GEN_Channel_with_Four_Cylinders_Static_Temp_Only.exo`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on (`isGhost = true`), integration-mesh Exodus output.

**Physics (FEM).** Named FEM phases and `neighbor_phases` replace explicit set names
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).
One `DIFF_LIN_ISO` model on `PhaseFluid`; IWGs `SPATIALDIFF_BULK`,
`SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` on the boundary to `PhaseInlet` (`DIRICHLET_NITSCHE`,
100.0), `SPATIALDIFF_NEUMANN` on the boundary to `PhaseCylinder`, and `GHOST_NORMAL_FIELD` with
`GHOST_DISPL` (0.005). One `DOF` IQI for TEMP.

**Solver (SOL).** The linear solver depends on `gTestCaseIndex`: 0 selects `PETSC` and sets the
solver warehouse `SOL_TPL_Type` to `Petsc`; 1 selects `AMESOS_IMPL`. Newton with relative drop
10⁻⁴, relaxation 1.0, at most 20 iterations, one static step, no preconditioner. MSI sets
`order_adofs_by_host = false`.

**Output (VIS).** `Channel_with_Four_Cylinders_Static_Temp_Only.exo` on
`HMR_dummy_n_p0,HMR_dummy_c_p0`: nodal TEMP and the global `IQIBulkTEMP`.

## What the test checks

All checks open the Exodus file; values are compared to relative tolerance 10⁻⁴ (marked FIXME as
loose because of serial/parallel differences), coordinates and the time value (1.0) to 10⁻⁸.
On 2 processors only rank 1 checks.

| Case | Procs | Nodes | Elements | Node | TEMP | IQI |
|---|---|---|---|---|---|---|
| linear | 1 | 6728 | 5754 | 1798 | 8.671607168734002e+04 | 8.326070772043242e+04 |
| linear | 2 | 4359 | 3714 | 3519 | same | same |
| quadratic | 1 | 18711 | 6573 | 8583 | 8.673830875995508e+04 | 8.329556290215004e+04 |
| quadratic | 2 | 12065 | 4233 | 12034 | same | same |

The reference node lies at (0.76, 0.24227) in every case. The Catch2 case
`Temp_Only_Channel_with_Four_Cylinders_Static_Linear` runs the PETSc solve and checks it only
when `MORIS_HAVE_PETSC` is defined, then always runs and checks the Amesos solve.
`Temp_Only_Channel_with_Four_Cylinders_Static_Quadratic` runs once with Amesos
(`gTestCaseIndex = 1`). The `Temp_Only_` prefix avoids a Catch2 name collision with the
advection example in the shared EXA runner.

## Run it

```bash
ctest -R Channel_with_Four_Cylinders_Static_Temp_Only-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
ctest -R Channel_with_Four_Cylinders_Static_Temp_Only-2-procs --output-on-failure
```

## Notes

- The `.dox` feature table lists "1 and 4 subdomains"; the CMakeLists registers 1 and 2
  processors, and the test accepts only those two counts.
- The directory contains `domain.png`, the domain sketch referenced by the `.dox` page.
