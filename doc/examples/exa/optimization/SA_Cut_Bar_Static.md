# Cut bar, steady heat conduction: sensitivity sweep

A finite-difference check of shape sensitivities on a steady two-material heat-conduction
problem. A thin vertical bar is split by a horizontal immersed interface into a conductive "fin"
on top and a heat-storage material below; the vertical position of the interface is the single
ADV. The SWEEP algorithm compares the analytical gradients of a maximum-temperature objective and
a volume constraint with forward, backward and central finite differences. The transient
counterpart with phase change is {doc}`SA_Cut_PCBar_Transient <SA_Cut_PCBar_Transient>`.

```{exa-summary} projects/EXA/optimization/SA_Cut_Bar_Static
```

## Problem

| | |
|---|---|
| Background domain | 0.16 × 0.6, offset (−0.0342356, −0.031345) |
| Interface | horizontal line y = s, with s the ADV: initial value and lower bound 0.425 (0.5 × 0.85), upper bound 0.475 (0.5 × 0.95) |
| Phases | phase 3 above the interface (fin), phase 2 below (storage material) |
| Fin material | conductivity 50.0, density 1.0, heat capacity 50.0 |
| Storage material | conductivity 1.0, density 1.0, heat capacity 50.0 |
| Heat flux | 1000.0 on side set 3 (top) of the fin |
| Prescribed temperature | 10.0 on side set 1 (bottom) of the storage material |
| Analysis | steady (default time-solver settings) |
| Objective | `MAX_DOF` of TEMP, parameters reference 10.0 / exponent 2.0 / shift 0.0 |
| Constraint | volume of the storage material |

## Setup in MORIS

**Mesh (HMR).** 1 × 40 linear elements, one Lagrange and one linear B-spline mesh (pattern 0),
refinement and staircase buffers 0, `adaptive_refinement_level` 0.

**Geometry and phases (GEN).** Two `LINE` geometries. The first (center x = −1.0, normal (1, 0))
lies outside the domain, so every point is on its positive side. The second has center
(0.0, s), normal (0, 1), and `center_y` set as (0.425, 0.425, 0.475), which makes it the only ADV.
With the default phase table, phase 3 is positive for both geometries and phase 2 is positive for
the first and negative for the second. `IQI_types` is `IQIMaxTemp`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization,
`high_to_low_dbl_side_sets`, and Exodus output of the enriched and integration meshes.

**Physics (FEM).** Legacy set-name input with explicit mesh sets
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)):
bulk `HMR_dummy_{n,c}_p3` (fin) and `HMR_dummy_{n,c}_p2` (storage), interface
`dbl_iside_p0_3_p1_2`, ghost sets `ghost_p3` and `ghost_p2`. Each material has a `DIFF_LIN_ISO`
constitutive model and a `SPATIALDIFF_BULK` IWG stabilized by a `GGLS_DIFFUSION` parameter.
The interface uses `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` (`NITSCHE_INTERFACE`, 100.0), the
heat flux `SPATIALDIFF_NEUMANN`, and the bottom temperature
`SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, 10.0). Ghost penalty is
`GHOST_NORMAL_FIELD` with `GHOST_DISPL` parameters of 0.01 per material. IQIs: `MAX_DOF` and a
nodal `DOF` IQI on both materials, and `VOLUME` on the storage material.

**Solver (SOL).** `BELOS_IMPL` with an `IFPACK` ILU preconditioner; Newton with relative residual
drop 10⁻⁷, relaxation 1.0, at most 10 iterations.

**Optimization (OPT).** `SWEEP` with `num_evaluations_per_adv = 1`, `include_bounds = false` and
`finite_difference_type = "all"`; `finite_difference_epsilons` is not set, so the sweep default
(10⁻⁸) applies. Results go to `SA_Cut_Bar_Static.hdf5`. Objective = criterion 0, constraint =
criterion 1.

**Output (VIS).** `SA_Cut_Bar_Static.exo` on both materials: nodal TEMP and global `MAX_DOF`.

## What the test checks

The test requires a return code of 0, then reads the analytical and FD gradients (forward,
backward, central) from `SA_Cut_Bar_Static.hdf5`. For every ADV it checks (`CHECK`) that the
relative difference between analytical and FD gradients is below 0.01 for the objective and
below 0.001 for the constraint. It runs on 1 processor.

## Run it

```bash
ctest -R SA_Cut_Bar_Static-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The GGLS parameters receive a latent heat of 0.0 and a phase-change temperature of 10000.0
  through dummy properties; the deck also defines `PropLatentHeat` (1000.0) and `PropPCTemp`
  (20.0), but no list references them.
- `tHaveGhost` and `tOutputCriterion` at the top of the deck toggle the ghost IWGs and the
  time-solver output for debugging.
