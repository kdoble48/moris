# Level-set design variable filter

How the linear filter of the level-set design variables (Kreissl and Maute 2012) works, why
level-set designs need it, and what it changes in an optimization. The level-set geometry
parameter `adv_filter_radius` turns it on. Cases 0 and 1 check the filtered sensitivities
against finite differences on a plate with a hole. Cases 2 and 3 optimize a cantilever without
and with the filter.

```{exa-summary} projects/EXA/optimization/LevelSet_ADV_Filter
```

## What the filter does

A level-set geometry discretized on B-splines (`discretization_mesh_index` ≥ 0) has one design
variable per B-spline coefficient. Without a filter, each optimizer variable is a coefficient.
With `adv_filter_radius` = $r$ > 0, the optimizer keeps its variables $\psi$ and the design
field uses filtered coefficients

$$
\hat\phi_i = \sum_j F_{ij}\,\psi_j, \qquad
F_{ij} = \frac{\max(0,\; r - \lVert x_i - x_j \rVert)}{\sum_k \max(0,\; r - \lVert x_i - x_k \rVert)},
$$

where $x_i$ are the Greville points of the coefficients (the cell centres for uniform quadratic
B-splines). The weights fall off linearly with distance, reach zero at the radius, and each row
sums to one, so a constant field stays constant and a linear field is reproduced wherever the
filter's neighbourhood is complete. The criteria sensitivities with respect to the optimizer's
variables are $F^\mathsf{T}$ times the sensitivities with respect to the coefficients. The
optimizer's bounds stay on $\psi$; because the weights are non-negative and sum to one,
$\hat\phi$ stays inside the same bounds.

```{figure} img/LevelSet_ADV_Filter_influence.png
:alt: One design variable's share given to its neighbours, and its effect on the design field, without the filter and with filter radii of 1.5 and 2 cells.

How far one design variable reaches. Top row: the design variables sit one per cell centre
(dots); the dashed circle is the filter radius; each coloured dot is the share of the centre
variable (red +) it receives. At radius 1.5 cells the variable keeps 39 % and passes 13 % to
each side neighbour and 2 % to each diagonal one; at 2 cells the shares are 24 %, 12 % and 7 %.
Bottom row: the change of the design field when that one variable changes by 1. Without the
filter it is one quadratic B-spline, reaching 1.4 cells from the centre. With the filter the
same change spreads over 2.2 cells with a half (1.5 cells) or a third (2 cells) of the peak.
```

## Why level-set designs need it

The optimizer moves every coefficient on its own. Even a uniform normal velocity, such as the
volume gradient that shrinks every member equally, is then applied unevenly: each coefficient's
gradient is the integral of its own basis function along the boundary, which depends on where
the boundary crosses the coefficient's support. On a boundary inclined to the grid this
modulation repeats every $h/\sin\theta$, so coefficient-wise steps write ripples into the
design field. A boundary aligned with the grid shows none.

Straight bars at angle $\theta$ to the grid, shrunk by pure volume descent in coefficient space
for one cell of boundary travel (RMS ripple of the boundary, in cells):

| Bar angle | Ripple wavelength | No filter | Filter 1.5 cells | Filter 2.5 cells |
|---|---|---|---|---|
| 0° | none | 0 | 0 | 0 |
| 15° | 4–8 cells | 0.030 | 0.007 | 0.0012 |
| 30° | 2–4 cells | 0.015 | 0.003 | 0.0008 |
| 45° | 1–2 cells | 0.002 | 0.002 | 0.0004 |

With the filter, one variable cannot make a bump narrower than about four cells, so the step
is smooth along the boundary. Barrera, Geiss and Maute (2022) use radius 1.5 cells in 2D and
1.8 cells in 3D.

## Problem

| | Cases 0, 1: sensitivity check | Cases 2, 3: optimization |
|---|---|---|
| Domain | 2 × 1 plate on [0, 2] × [0, 1] | same |
| Mesh | 12 × 6 cells (h = 1/6) | 40 × 20 cells (h = 1/20) |
| Design | one circular hole, centre (1.0, 0.5), radius 0.3; 112 design variables | staggered grid of 13 holes, radius 0.12; 924 design variables |
| Material | linear elastic, E = 1, ν = 0.3; material where the level set is positive | same |
| Supports | left edge x = 0 clamped (Nitsche) | same |
| Load | traction (0, −1) on the right edge | traction (0, −10) on the right edge, y ∈ [0.45, 0.55] |
| Criteria | strain energy and volume | strain energy and volume |
| Optimization | none (SWEEP at the initial design) | minimum strain energy, volume ≤ 1.0 (half the domain) |
| Filter | case 0: off; case 1: radius 1.5 h | case 2: off; case 3: radius 1.5 h |

## Setup in MORIS

**Mesh (HMR).** One quadratic Lagrange mesh and one quadratic B-spline mesh (index 0), on 12 × 6
or 40 × 20 cells. The displacement DOFs and the level set share B-spline mesh 0.

**Geometry (GEN).** One level-set geometry on B-spline mesh 0 (`discretization_mesh_index` = 0).
Cases 0 and 1 use a `CIRCLE` field. Cases 2 and 3 use a `USER_DEFINED` field
`Func_Hole_Grid`, the distance to the nearest of 13 hole centres minus the radius, clipped to
the coefficient bounds ±3h. They also set `use_cell_order_intersection`, so the cut follows the
quadratic level set inside each cell rather than its corner values. Cases 1 and 3 set
`adv_filter_radius` to 1.5 h. `IQI_types` is `IQIBulkStrainEnergy`, `IQIBulkVolume`.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization.

**Physics (FEM).** `STRUC_LIN_ISO` with `STRUC_LINEAR_BULK`, a symmetric Nitsche Dirichlet IWG
on side set 4 (`DIRICHLET_NITSCHE` penalty 100), a Neumann traction on side set 2, and
`GHOST_NORMAL_FIELD` with a `GHOST_DISPL` parameter of 0.01. Criteria: a `STRAIN_ENERGY` IQI
and a `VOLUME` IQI.

**Solver (SOL).** Direct `AMESOS_IMPL` solve with `Amesos_Klu`, one Newton iteration.

**Optimization (OPT).** Cases 0 and 1: `SWEEP` with `finite_difference_type = "all"` and
`finite_difference_epsilons = 1E-6` for every design variable (a single finite-difference mode
skips the analytic gradients). Cases 2 and 3: the in-tree `MMA` algorithm with its
conservative inner loop (`step_size` 0.05, `penalty` 5, `max_inner_its` 5, 100 iterations),
objective strain energy divided by its initial value, constraint volume / 1.0 − 1. No
level-set reinitialization.

## Results

### Sensitivities

| | Case 0, no filter | Case 1, filter 1.5 h |
|---|---|---|
| Design variables with a nonzero volume gradient | 36 | 64 |
| Worst analytic vs central-difference error, objective | 5.2 × 10⁻⁷ | 8.2 × 10⁻⁷ |
| Worst analytic vs central-difference error, constraint | 1.7 × 10⁻⁹ | 4.7 × 10⁻⁹ |
| Strain energy at the initial design | 23.773 | 22.717 |
| Volume at the initial design | 1.7243 | 1.7568 |

Errors are relative to the largest analytic entry. The filter couples each variable to its
neighbours, so more variables have a nonzero gradient, and the initial design is the filtered
field, so its criteria differ from case 0.

### Optimization

```{figure} img/LevelSet_ADV_Filter_designs.png
:alt: Cantilever starting design with 13 holes, and the final designs after 100 iterations without and with the filter, with zooms on the cell grid.

Top: the starting design and the final designs after 100 iterations without (case 2) and with
(case 3) the filter. Bottom: the same window of both final designs on the cell grid. Without
the filter the members have kinked, ragged edges and a small stray hole. With the filter the
members are smooth and the layout is simpler.
```

| | Case 2, no filter | Case 3, filter 1.5 h |
|---|---|---|
| Evaluations (100 iterations) | 324 | 340 |
| Final volume (limit 1.0) | 1.000 | 0.9996 |
| Final strain energy | 31.91 | 32.07 |

At the same volume the filtered design is 0.5 % more compliant, a small price for the smoother
boundary and simpler layout. Neither run is converged after 100 iterations; the example shows
the trend, not optimal designs.

## Using the filter in a deck

```cpp
aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
aParameterLists.set( "field_function_name", "Func_My_Seed" );
aParameterLists.set( "discretization_mesh_index", 0 );           // the filter acts on B-spline coefficients
aParameterLists.set( "discretization_lower_bound", -3.0 * h );
aParameterLists.set( "discretization_upper_bound", 3.0 * h );
aParameterLists.set( "adv_filter_radius", 1.5 * h );            // physical units; 0 = off (default)
```

- Use 1.5 cells in 2D and about 1.8 cells in 3D. A radius below the coefficient spacing changes
  nothing (each coefficient is its own only neighbour).
- Keep the radius in physical units when the mesh is refined.
- The filter runs on one processor only.
- The in-place level-set reinitialization (`reinit_in_place`) writes coefficient values directly
  and is not supported with the filter; the workflow stops with an error if both are set.
- Restart files hold the optimizer's variables ψ, so a warm restart continues the filtered run.

## What the test checks

- **Cases 0 and 1:** analytic and central-difference sensitivities agree for every design
  variable, within 10⁻⁵ of the largest analytic entry, for both criteria; the criteria at the
  initial design match reference values.
- **Cases 2 and 3:** the final design meets the volume limit, and its strain energy and volume
  match reference values (relative tolerance 10⁻⁶).

The filter itself is unit-tested by `GEN-test.exe "[adv_filter]"`: identity when off, constants
preserved, linear fields reproduced away from the edges, $F^\mathsf{T}$ the exact adjoint of
$F$, and the design field equal to an unfiltered field given $F\psi$.

## Run it

```bash
ctest -R LevelSet_ADV_Filter --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The four ctests are `LevelSet_ADV_Filter_<0..3>-1-procs`; cases 0 and 1 take about 5 s each,
cases 2 and 3 about 80 s each.

## References

- Kreissl, S. and Maute, K. (2012). Levelset based fluid topology optimization using the
  extended finite element method. *Structural and Multidisciplinary Optimization* 46(3),
  311–326. [doi:10.1007/s00158-012-0782-8](https://doi.org/10.1007/s00158-012-0782-8)
- Barrera, J. L., Geiss, M. J. and Maute, K. (2022). Minimum feature size control in level set
  topology optimization via density fields. *Structural and Multidisciplinary Optimization*
  65(3), 94. [doi:10.1007/s00158-021-03096-7](https://doi.org/10.1007/s00158-021-03096-7).
  The filter is their Eq. (1); they use radius 1.5 h in 2D and 1.8 h in 3D.
