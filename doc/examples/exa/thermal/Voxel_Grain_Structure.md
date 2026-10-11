# Heat conduction through a voxel grain structure (3D)

A steady 3D heat-conduction problem through a polycrystal-like microstructure: ten grains are
read from a 64 × 64 × 64 voxel file, a phase function maps them to phases 1–10, and phases 1–5
and 6–10 are assigned two different conductivities. A temperature is imposed on one face of the cube and
a heat flux on the opposite face; neighbouring grains are coupled by Nitsche's method across
every grain boundary. The deck demonstrates the GEN voxel geometry and a phase function that maps
voxel grain IDs to XTK phases. A related deck is {doc}`Grain_Structure <Grain_Structure>`.

```{exa-summary} projects/EXA/thermal/diffusion/Voxel_Grain_Structure
```

## Problem

| | |
|---|---|
| Domain | cube [−2, 2]³ |
| Microstructure | `Voxel_file`: 262144 lines (64³) of five integers; the second column takes values 1–10, the last three run from 0 to 63 |
| Phases | 11 (`number_of_phases`): 1–10 from `get_phase_index`, 0 otherwise |
| Conductivity | 1.0 for phases 1–5, 10.0 for phases 6–10 |
| Density, heat capacity | 1.0 |
| Right face (`SideSet_2`) | T = 2.0, symmetric Nitsche |
| Left face (`SideSet_4`) | heat flux 2.0 (Neumann) |
| Time | static, default time-solver settings |

## Setup in MORIS

**Mesh (HMR).** 39 × 39 × 39 linear elements on [−2, 2]³, linear B-splines, refinement and
staircase buffer 1.

**Geometry and phases (GEN).** One voxel geometry (`create_voxel_geometry_parameter_list`)
reading `<MORIS root>/projects/EXA/thermal/diffusion/Voxel_Grain_Structure/Voxel_file` (path from
`moris::get_base_moris_dir()`), with the same domain dimensions and offset as HMR.
`get_phase_index` returns the first index 1 … 10 whose bit is set in the geometry bitset, else 0.

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
and multigrid off (`isGhost`, `isMultigrid`), Exodus output of the integration mesh.

**Physics (FEM).** Legacy set-name input. `Name_func()` builds the set-name lists at run time
for the two grain groups ([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)):
bulk sets `HMR_dummy_{n,c}_p1` … `p5` and `p6` … `p10`, the side sets `SideSet_2_{c,n}_pN` and
`SideSet_4_{c,n}_pN`, and the double-sided grain-boundary sets `dbl_iside_p0_i_p1_j` (i < j)
sorted into group pairs 1–1, 1–2 and 2–2. Each group has a `DIFF_LIN_ISO` model
(`CMDiffusion_1`, conductivity 1.0; `CMDiffusion_3`, conductivity 10.0) and a `SPATIALDIFF_BULK`
IWG. The right face uses `SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE` (`DIRICHLET_NITSCHE`, penalty
100.0) per group; the left face uses `SPATIALDIFF_NEUMANN` on all ten phases. Grain boundaries use
`SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE` with a `NITSCHE_INTERFACE` parameter (penalty 100.0)
per group pair. The only IQI is the `DOF` temperature.

**Solver (SOL).** Direct solve with `Amesos_Mumps` (MUMPS builds) or `Amesos_Superludist`; one
Newton iteration, relative residual drop 10⁻⁴.

**Output (VIS).** `Voxel_Grain_Structure.exo` with nodal `TEMP` on all grain bulk sets.

## What the test checks

The Catch2 case `Voxel_Grain_Structure` sets `gInterpolationOrder = 1`, `gDim = 2` and
`gTestCaseIndex = 0` (the deck reads neither `gDim` nor `gTestCaseIndex`), calls the workflow
without checking the return code, and compares `Voxel_Grain_Structure.exo`:

- 3 dimensions (`REQUIRE`);
- 404689 nodes and 656001 elements (`CHECK`, so a mismatch is reported but the test continues);
- node 127193 at (−1.282051282051282, 1.692307692307692, −0.3589743589743590), relative
  tolerance 10⁻⁸;
- nodal temperature 4.03803 at node 127193, relative tolerance 10⁻⁴.

## Run it

```bash
ctest -R Voxel_Grain_Structure-1-procs --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The deck locates `Voxel_file` (4.8 MB) through the MORIS source directory, so it runs from any
working directory as long as `get_base_moris_dir()` resolves.

## Notes

- The interface constitutive-model pairs in the FEM block reference `CMDiffusion_2` and
  `CMDiffusion_4` as followers; these carry the same conductivities as `CMDiffusion_1` (1.0) and
  `CMDiffusion_3` (10.0).
- The 2–1 interface list (`tInterface_21`) stays empty: its loop visits only pairs with the
  second grain index at least the first, which never places a group-2 grain before a group-1
  grain.
- The `Func_Inlet_U` callback (a parabolic velocity profile) is defined but not referenced.
- The `.dox` documentation file is absent from this directory.
