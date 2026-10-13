# Oscillator

Transient incompressible flow of water through a 2D fluidic oscillator on a body-fitted
mesh. The mesh is read from an Exodus file through STK (workflow `STK_FEM`), so the deck uses
no background mesh, geometry or cut; boundary conditions are applied on the sidesets stored
in the mesh file. The example demonstrates the body-fitted path through FEM and SOL with
SUPG/PSPG-stabilized incompressible Navier–Stokes, Nitsche boundary conditions, backflow
prevention at the outlet, and time continuity.

```{exa-summary} projects/EXA/fluid/laminar/Oscillator
```

## Problem

| | |
|---|---|
| Domain | fluidic-oscillator geometry from `Oscillator.g` (2D, 20 element blocks `block_1`…`block_20`) |
| Material | water: density 999, viscosity 1.12·10⁻³ (all scales 1.0, so dimensional values are used) |
| Inlet | sidesets `surface_1`, `surface_2`: parabolic profile u_x = 0.5 (1 − 4y²/H²), H = 3·10⁻³, u_y = 0, multiplied by a tanh time ramp over 0.25 |
| Walls | sidesets `surface_5`, `surface_6`: no-slip |
| Outlet | sidesets `surface_3`, `surface_4`: backflow prevention (weight 1.0) |
| Time | 2 time steps of 0.005, final time 0.010; initial velocity (0, 0) in the time-continuity term |

The ramp factor is `clamp( (tanh(βη) + tanh(β(t/T_r − η))) / (tanh(βη) + tanh(β(1 − η))), 0, 1 )`
with β = 4, η = 0.5 and T_r = 0.25. A deck comment gives the original final time as 26.0;
the checked-in value 0.010 covers a small fraction of the ramp. The deck and the `.dox` page
both list a domain thickness of 3·10⁻³, which enters only the unused `tDarcy` constant.

## Setup in MORIS

**Mesh (STK).** OPT sets `workflow` to `STK_FEM`. The STK parameter list reads
`<MORIS root>/projects/EXA/fluid/laminar/Oscillator/Oscillator.g`, resolved with
`get_base_moris_dir()`. The `.dox` page describes how that file was produced: COMSOL output
(`sharpGeom_studentMesh.nas`) converted to Exodus by `ProcessMesh.m` and `DataExtractor.m`,
with sidesets added in Cubit (`Oscillator.cub`). The STK callback also sets the logger
severity to 0.

**Physics (FEM).** The block uses `hack_for_legacy_fem()` and explicit `mesh_set_names`.
Constitutive model `FLUID_INCOMPRESSIBLE`. Stabilization parameters `INCOMPRESSIBLE_FLOW`
(`SPSUPGNS`, "36.0 / 1.0") and `VELOCITY_DIRICHLET_NITSCHE` (`SPNitscheU`, "50.0/1.0"). IWGs:

- `INCOMPRESSIBLE_NS_VELOCITY_BULK` and `INCOMPRESSIBLE_NS_PRESSURE_BULK` on all 20 blocks;
- `INCOMPRESSIBLE_NS_VELOCITY_DIRICHLET_UNSYMMETRIC_NITSCHE` and
  `INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_UNSYMMETRIC_NITSCHE` on the wall sidesets
  (`PropNoSlip`) and on the inlet sidesets (`PropInletVel`, callback `Func_InletVel`);
- `INCOMPRESSIBLE_NS_IMPOSED_PRESSURE` on the outlet sidesets with `PropDensity` and
  `PropBackFlowPrevention` (1.0);
- `TIME_CONTINUITY_DOF` on VX, VY with weights ρ · 1.0 = 999 for the current and previous
  step and initial condition `PropInitialVel` (0, 0).

IQIs: `DOF` for VX, VY and P on all blocks.

**Solver (SOL).** `AMESOS_IMPL` with `Solver_Type` `Amesos_Mumps` when MORIS is built with
`MORIS_USE_MUMPS`, otherwise `Amesos_Superludist`. Newton (`NEWTON_SOLVER`) with
`NLA_rel_res_norm_drop` 10⁻³, relaxation 1.0, at most 20 iterations, on `VX,VY;P`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). Time solver:
`TSA_Num_Time_Steps` 2, `TSA_Time_Frame` 0.010, two time levels per DOF type, initial solution
VX = 0.5, VY = 0, P = 0. No preconditioner.

**Output (VIS).** `Oscillator.exo` on all blocks with nodal VX, VY, P, written every step
(`Save_Frequency` 1, `Time_Offset` 10.0).

## What the test checks

The test sets `gInterpolationOrder = 1` and `gPrintReferenceValues = false`, runs the
workflow with `./Oscillator.so`, and on 1 processor reads `Oscillator.exo.e-s.0000` and
compares the first stored time step:

| Quantity | Reference | Tolerance |
|---|---|---|
| Dimensions, nodes, elements | 2, 7127, 11454 | exact |
| Coordinates of node 4403 | (0.0274990890175104, 0.00591698987409472) | relative norm < 10⁻⁵ |
| Time value | 0.005 | relative < 10⁻⁸ |
| Velocity (VX, VY) at node 4403 | (−3.14130625281294·10⁻⁷, −2.69800723469803·10⁻⁷) | relative norm < 10⁻⁵ |
| Pressure at node 4403 | 0.185522355835879 | relative < 10⁻⁵ |

The return code is not checked. The test file also contains a branch for 4 processors
(checked on rank 1) and aborts with an error for other processor counts; the CMake
registration uses 1 processor only. A pressure-IQI check is present but commented out.

## Run it

```bash
ctest -R Oscillator-1-procs --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

The mesh path is built from the MORIS root directory returned by `get_base_moris_dir()`.

## Notes

- `tAlpha` (an effective impermeability for a Darcy term) and `tDarcy` are computed but not
  used by any property or IWG.
- The directory also holds MATLAB scripts (`DataExtractor.m`, `ProcessMesh.m`, `misc/`) for
  regenerating the Exodus mesh; the `.dox` page refers to them under an older path
  `projects/EXA/fluid/BodyFitted/Oscillator/`.
