# Fickian problem for linear diffusion

A transient heat-diffusion problem following Fick's second law on a thin bar: an initially
uniform bar at 313 is heated by a wall temperature of 350 imposed at one end, all other
boundaries are adiabatic. The 1D semi-infinite version of this problem has a closed-form
error-function solution, which the `.dox` page and two MATLAB scripts in the directory provide
for comparison. The deck exercises time-continuity IWGs, GGLS stabilization of the transient
diffusion operator, and a Nitsche-enforced, time-ramped Dirichlet condition, in linear and
quadratic interpolation.

```{exa-summary} projects/EXA/thermal/diffusion/Fickian_Problem
```

## Problem

| | |
|---|---|
| Domain | 0.168 × 0.0005 bar, 300 × 1 elements |
| Material | density ρ = 0.75, heat capacity c_p = 2.4, conductivity κ = 2.1 · 10⁻⁷ |
| Initial temperature | T₀ = 313 |
| Wall temperature | T_w = 350 on side set 4, ramped linearly from 313 over the first three time steps |
| Other boundaries | no condition written (adiabatic) |
| Time | 20 time steps over a time frame of 2000 (Δt = 100) |

The ramp is `min(350, 313 + 37 t / (3 · 2000 / 20))`; the `.dox` developer notes give the
reason: imposing the full Dirichlet value from the first step leads to instabilities.

The analytical solution quoted on the `.dox` page is
T(x, t) = T_w + (T₀ − T_w) erf(x / (2 √(α t))) with α = κ / (ρ c_p), x measured from the heated
end. `solve_fick_problem_1D.m` evaluates it; `Fick_Problem.m` plots it for a length of 0.25 and
t_max = 16/3 h.

## Setup in MORIS

**Mesh (HMR).** 0.168 × 0.0005 domain with 300 × 1 elements, `processor_decomposition_method = 1`.
Lagrange and B-spline orders are 1 or 2 from `gInterpolationOrder` (any other value is an error);
refinement and staircase buffers 3, `adaptive_refinement_level = 1`.

**Geometry and phases (GEN).** One `USER_DEFINED` level set, φ = y − 10000, refined once on
mesh 0. It is negative on the whole domain, so nothing is cut and all elements belong to phase 0
(the deck's comment shows the alternative `tPlaneBottom - 0.01`).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
and `high_to_low_dbl_side_sets` on, integration-mesh Exodus output.

**Physics (FEM).** `DIFF_LIN_ISO` with the three material properties. IWGs on
`HMR_dummy_n_p0`: `SPATIALDIFF_BULK` with the `GGLS_DIFFUSION` parameter, and
`TIME_CONTINUITY_DOF` (weights 100.0 / 100.0, initial condition 313, `time_continuity = true`).
The wall temperature is `SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE` on `SideSet_4_n_p0` with
`DIRICHLET_NITSCHE` penalty 1000.0. A `GHOST_DISPL` parameter (`SPGPTemp`, 0.01) is defined
but no IWG references it. One `DOF` IQI for TEMP.

**Solver (SOL).** `AMESOS_IMPL`, Newton with relative drop 10⁻⁶, relaxation 1.0, at most 5
iterations; 20 time steps over 2000 with `TSA_time_level_per_type = TEMP,2`.

**Output (VIS).** `Fick_Problem_<gTestCaseIndex>.exo` (`_0` linear, `_1` quadratic) on
`HMR_dummy_n_p0` with nodal TEMP.

## What the test checks

Two Catch2 cases, serial only (any other processor count raises an error). Each requires a
return code of 0, then reads the Exodus file:

| Case | File | Nodes | Elements | Node | TEMP at time step index 19 |
|---|---|---|---|---|---|
| `Fickian_Problem_Linear` (order 1) | `Fick_Problem_0.exo` | 602 | 300 | 14 | 3.442925667877556e+02 |
| `Fickian_Problem_Quadratic` (order 2) | `Fick_Problem_1.exo` | 1803 | 300 | 39 | 3.442923942621712e+02 |

Both check 2 dimensions, the reference node at (3.92 · 10⁻³, 0) to 10⁻⁸, the time value 100 to
10⁻⁸, and the temperature to relative tolerance 10⁻⁵. The references are stored MORIS results;
the test does not evaluate the error-function solution.

## Run it

```bash
ctest -R Fick_Problem --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest name is `Fick_Problem` (`TEST_BASE` with `NO_PROCS_SUFFIX`), not the directory name.

## Notes

- The `.dox` page gives a time frame of 960 split into 20 time slabs; the deck uses 2000.
- The linear test case calls `gLogger.initialize( "Log.log" )` before running.
