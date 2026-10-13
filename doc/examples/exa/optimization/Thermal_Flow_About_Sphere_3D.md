# Thermal flow about a sphere (3D)

The 3D counterpart of {doc}`Thermal_Flow_About_Sphere_2D`: incompressible laminar flow through a
channel of square cross-section passes a heated spherical solid inclusion, with conjugate heat
transfer between fluid and solid. The sphere radius is the single ADV and the OPT block is a
finite-difference `SWEEP`. The directory is not connected to the build: the parent
`CMakeLists.txt` comments out its `add_subdirectory`, and there is no `example_test_case.cpp`.

```{exa-summary} projects/EXA/optimization/Thermal_Flow_About_Sphere_3D
```

## Problem

Inputs are in SI units and scaled with reference length 1 mm, reference pressure equal to the
inlet pressure, and reference density 0.3289 kg/m³, as in the 2D deck.

| | |
|---|---|
| Channel | 5.0 mm × 1.0 mm × 1.0 mm (fluid region, shifted by 0.02 mm from the background origin) |
| Inclusion | sphere at (0.75, 0.5, 0.5) mm, radius 0.225 mm |
| Phases | 0 fluid, 1 solid, 2 inlet void, 3 outlet void, 4 lateral void (four lateral planes) |
| Fluid | μ = 43.32·10⁻⁶ N s/m², ρ = 0.3289 kg/m³, c_p = 1210 J/(kg K), k = 76.26·10⁻³ W/(m K) |
| Solid | ρ = 7000 kg/m³, c_p = 200 J/(kg K), k = 150 W/(m K) |
| Inlet | imposed pressure 101.325 Pa; temperature 0.0 |
| Walls | zero velocity on the four lateral walls and the sphere surface (Nitsche) |
| Heat source | volumetric load 250·10⁴ (deck units W/m², scaled as a volumetric load) in the solid |
| Analysis | steady, one time step |

The reference values give Re ≈ 133, as in 2D.

## Setup in MORIS

**Mesh (HMR).** A 5.25 × 1.25 × 1.25 mm background domain offset by −0.125 mm in each
direction, with 27 × 7 × 7 elements (`ceil(dim / 0.2)`), no initial refinement, and two
B-spline meshes plus one Lagrange mesh. The interpolation order is the file-local
`gInterpolationOrder = 2`. MSI places VX, VY, VZ and TEMP on B-spline mesh 1.

**Geometry and phases (GEN).** Seven `USER_DEFINED` level sets: `Func_Sphere` and six
`Func_Plane` instances (inlet, outlet, lower, upper, back, front), each offset 0.02 mm off the
grid. The radius is the `Design_Variable` (bounds 0.9r–1.1r, r = 0.225), with
`Func_Sphere_Deriv` as its field sensitivity. `get_phase_index` maps the seven geometry signs to
the five phases, which FEM refers to by name.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization, and Exodus output of the integration mesh.

**Physics (FEM).** The same IWG/IQI/SP set as the 2D deck, extended to VZ (`FLUID_INCOMPRESSIBLE`
with SUPG/PSPG `INCOMPRESSIBLE_FLOW`, `SPATIALDIFF_BULK` and `ADVECTION_BULK` with
`SUPG_ADVECTION`, Nitsche no-slip and inlet conditions, and `GHOST_NORMAL_FIELD` ghost terms).
Differences from 2D:

| | 2D deck | 3D deck |
|---|---|---|
| Fluid–solid thermal interface IWG | `SPATIALDIFF_INTERFACE_SYMMETRIC_NITSCHE` | `SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE` |
| `INCOMPRESSIBLE_FLOW` parameters | `C_inv/1.0` | `C_inv` |
| `SUPG_ADVECTION` parameters | `1.0` | not set |
| Thermal-energy IQIs constitutive model | `CMFluidDiffusion` | `CMFluid` |
| Unused selection properties | `PropSelectUY` | `PropSelectUY`, `PropSelectUZ` |
| Output files | `Thermal_Flow_About_Sphere_2D_<case>.*` | `Channel3D.exo`, `Channel3D_SEN.hdf5`, `Channel3D_GEN.exo` |

**Solver (SOL).** As in 2D: Amesos (`Amesos_Mumps` or `Amesos_Superludist`), a one-way staged
`NLBGS_SOLVER` over the flow block (VX, VY, VZ, P; Newton, drop 10⁻⁸, up to 20 iterations) and
the TEMP block for the forward problem, and a monolithic Newton solver for the sensitivities.

**Optimization (OPT).** `SWEEP` with `finite_difference_type = all`, epsilon 10⁻⁴, one
evaluation per ADV. The objective is criterion 0; nine constraints map one to one to the
criteria.

**Output (VIS).** `Channel3D.exo` with nodal VX, VY, VZ, P and TEMP on the fluid and solid
sets, time offset 10.

## What the test checks

Nothing: the leaf has no `example_test_case.cpp` and is not added to the build. The
`CMakeLists.txt` declares `PROCS 1 4`, which would yield `Thermal_Flow_About_Sphere_3D-1-procs`
and `Thermal_Flow_About_Sphere_3D-4-procs` once a test translation unit is written and the
subdirectory is reconnected.

## Run it

No ctest target exists. The `CMakeLists.txt` comment states that reconnecting requires authoring
the test TU first.

## Notes

- The deck builds its parameter lists with the older `aParameterLists( i ).push_back(
  prm::create_*_parameter_list() )` calls and indexes FEM lists through `tPropIndex`,
  `tIQIIndex` and similar names that are not declared in the file, unlike the 2D deck, which
  uses `add_parameter_list()` with `FEM::`/`SOL::` sub-list names.
- `IQI_types` lacks a comma after `"IQIOutletTotalPressure"`, so the adjacent literals
  concatenate to `"IQIOutletTotalPressureIQIPerimeterItf"` and the list holds eight entries,
  while `tNumConstraints` is 9 and the callbacks index nine criteria
  ([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).
