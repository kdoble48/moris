# Channel with four cylinders, transient

A transient version of {doc}`Channel_with_Four_Cylinders_Static <Channel_with_Four_Cylinders_Static>`:
incompressible flow and heat transport in a 2D channel past four immersed cylinders, solved over
10 time steps with a ramped inlet velocity. Flow and temperature are solved in a staggered
nonlinear block Gauss–Seidel (NLBGS) loop rather than monolithically. The deck is written in the
pre-modern parameter-list API and is currently disconnected from the build (see Notes
below).

```{exa-summary} projects/EXA/thermal/advection/Channel_with_Four_Cylinders_Transient
```

## Problem

| | |
|---|---|
| Background domain | 4.0 × 2.0 rectangle, offset (−1.24, −0.86) |
| Physical domain | channel 0 ≤ x ≤ 2.2, 0 ≤ y ≤ 0.41 |
| Cylinders | radius 0.05, centers x = 0.4, 0.8, 1.2, 1.6 at y = 0.1, 0.3, 0.1, 0.3 |
| Fluid | viscosity 0.001, density 1.0, heat capacity 1.0, conductivity 0.0005 |
| Inlet | parabolic u_x = 4 · m(t) · 1.5 · y (0.41 − y) / 0.41², temperature 0.0 |
| Ramp | m(t) = min(t / Δt / 5, 1): full inlet velocity after 5 time steps |
| Walls and cylinders | no-slip; cylinder heat flux 20.0 |
| Initial conditions | velocity property (10⁻⁵, 0); temperature 0.0 |
| Time | 10 time steps over a time frame of 10.0 (Δt = 1.0) |

Compared with the static deck, the conductivity is 10× larger (0.0005 vs 0.00005) and the
cylinder flux is 10× larger (20.0 vs 2.0).

## Setup in MORIS

**Mesh (HMR).** 22 × 11 elements with four levels of initial refinement
(`pattern_initial_refinement = 4`), refinement buffer 3, staircase buffer 5,
`adaptive_refinement_level = 1`, `use_number_aura = false`; linear Lagrange and B-spline meshes.

**Geometry and phases (GEN).** Eight `USER_DEFINED` level sets (four planes, four cylinders),
each clamped to a minimum magnitude of 10⁻⁶. No phase function is given, so phases are the
default sign-pattern indices; the deck refers to them directly in its set names (fluid `p160`).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost stabilization
on; no integration-mesh Exodus output.

**Physics (FEM).** The same constitutive models, bulk IWGs, Nitsche boundary IWGs, ghost
penalties, and `DOF` IQIs as the static deck, addressed by explicit mesh-set names
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)):
bulk `HMR_dummy_c_p160,HMR_dummy_n_p160`, inlet `iside_b0_160_b1_128`, cylinders
`iside_b0_160_b1_{168,164,162,161}`, walls `iside_b0_160_b1_{32,224}`, ghost `ghost_p160`.
Differences from the static deck: the pressure Dirichlet IWGs use
`INCOMPRESSIBLE_NS_PRESSURE_DIRICHLET_SYMMETRIC_NITSCHE`, the `INCOMPRESSIBLE_FLOW` parameter is
`36.0`, `SUPG_ADVECTION` takes no function parameter, and two `TIME_CONTINUITY_DOF` IWGs (VX, VY
and TEMP; weights 100.0/100.0, `time_continuity = true`) enforce the initial conditions between
time slabs.

**Solver (SOL).** `BELOS_IMPL` with `ifpack_prec_type = ILU`. Nonlinear solver 0 (Newton,
VX, VY, P) and solver 1 (Newton, TEMP) share algorithm 0 (relative drop 10⁻⁴, relaxation 1.0,
at most 10 iterations); solver 2 is `NLBGS_SOLVER` over sub-solvers 0 and 1 with algorithm 1
(one iteration, drop 1.0). The time solver uses solver 2, 10 steps over 10.0, initial vector
`VX,1.5;VY,0.0;P,0.0;TEMP,0.0`
([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)).

**Output (VIS).** `Channel_with_Four_Cylinders_Transient.exo` on the fluid sets with nodal VX,
VY, P, TEMP, written every step.

## What the test checks

The single Catch2 case runs the deck and requires a return code of 0. It compares no field or
mesh values.

## Run it

The leaf registers `Channel_with_Four_Cylinders_Transient` (one process, `NO_PROCS_SUFFIX`), but
`projects/EXA/thermal/advection/CMakeLists.txt` has the `add_subdirectory` line commented out,
so the test does not exist in a build. After the deck is ported and the line is restored:

```bash
ctest -R Channel_with_Four_Cylinders_Transient --output-on-failure   # in a build with MORIS_USE_EXAMPLES=ON
```

## Notes

- The parent CMakeLists states: "these decks use the pre-modern ParameterList API and no longer
  compile", to be reconnected once ported to `Module_Parameter_Lists`. The deck's module
  functions take `Vector< Vector< ParameterList > >&` and append lists with `push_back`.
- `IWGZeroVelocity` and `IWGZeroPressure` are each defined twice with identical names, types,
  and set lists.
- `Func_Initial_Condition` sets the velocity to (10⁻⁵, 0) and does not read its
  `function_parameters` (1.5).
