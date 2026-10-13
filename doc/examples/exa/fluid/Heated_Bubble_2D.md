# Heated compressible bubble (transient)

A transient compressible Navier–Stokes simulation of a square box of air at rest, heated by
a radially symmetric volumetric source centred at the origin. All four sides carry a
prescribed-pressure Nitsche condition with pressure upwinding. The deck exercises the
compressible bulk IWG with a body heat load and GLS stabilization, and time-continuity IWGs
for space–time time stepping. It shares its structure with
{doc}`Heated_Channel_2D <Heated_Channel_2D>`; the differences are tabulated below.

```{exa-summary} projects/EXA/fluid/compressible/Heated_Bubble_2D
```

## Problem

| | |
|---|---|
| Domain | 50 × 50 square centred at the origin, x, y ∈ [−25, 25] |
| Phases | one fluid phase (phase 0) covering the whole domain |
| Material | perfect gas: c_v = 718, R = 287.058, μ = 1.716·10⁻⁵, κ = 24.35·10⁻³ |
| Heat load | Q(r) = 10⁵ · exp(−200 (2r/L)²), r = √(x² + y²), L = 50 |
| Boundary | pressure fixed to the initial pressure on all four sides (side ordinals 1–4), pressure upwind weight 3.774074960608205·10⁻³ |
| Initial state | ρ₀ = 1.0, T₀ = 273, P₀ = R ρ₀ T₀ ≈ 78 367, velocity 0 |
| Time | 30 time steps of Δt = 5/c_s ≈ 0.0151, total ≈ 0.453 |

As in the channel deck, c_s = √(γ R T₀) ≈ 331.2 with γ = 1 + R/c_v, and the step size is
`20 L / 100 / c_s / 2`.

### Differences from Heated_Channel_2D

| | Heated_Bubble_2D | Heated_Channel_2D |
|---|---|---|
| Domain, elements | 50 × 50, offset (−25, −25); 50 × 50 elements | 100 × 1; 100 × 1 elements |
| Heat-load distribution | radial Gaussian about the origin | Gaussian in x about mid-length |
| Side walls (1, 3) zero normal velocity | off (`tHaveTopBottomBCs = false`) | on |
| Pressure IWG sides | "1,2,3,4" | "2,4" |
| Pressure upwind weight | 3.774074960608205·10⁻³ | 1.0 |
| Nitsche penalty VX, VY | 10⁷ | 100 |
| GLS (`DummySP,GLS` on the bulk IWG) | on | off |
| Time penalty → weight | 10 → ≈ 662 | 0.1 → ≈ 3.31 |
| Time steps | 30 | 10 |
| Newton tolerance | 10⁻⁷ | 10⁻⁹ |
| Parameter-list API | `Vector< Vector< ParameterList > >` with `prm::create_*` | `Module_Parameter_Lists` |

## Setup in MORIS

**Mesh (HMR).** 50 × 50 elements, quadratic Lagrange and B-spline orders (`tIpOrder = 2`),
`lagrange_to_bspline = "-1"` (`tUseLagrange = true`).

**Geometry and phases (GEN).** One `USER_DEFINED` level set, `Func_Dummy_Plane`
(φ = y − 10000), which places the whole domain in phase 0; the fluid sets are
`HMR_dummy_c_p0, HMR_dummy_n_p0`
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, no ghost
stabilization; the enriched integration mesh is printed and written to Exodus.

**Physics (FEM).** Material model `PERFECT_GAS` and constitutive model
`FLUID_COMPRESSIBLE_NEWTONIAN`. IWGs: `COMPRESSIBLE_NS_BULK` with `PropHeatLoad` as
`BodyHeatLoad` and GLS through `DummySP` (a `DIRICHLET_NITSCHE` parameter of value 1.0);
`COMPRESSIBLE_NS_DIRICHLET_SYMMETRIC_NITSCHE` on sides "1,2,3,4" with `PropInitialPressure` as
`PrescribedDof1` and `PropUpwind` as `PressureUpwind`, stabilized by
`COMPRESSIBLE_DIRICHLET_NITSCHE` (penalties 0, 10⁷, 10⁷, 0); and three `TIME_CONTINUITY_DOF`
IWGs for P, VX/VY and TEMP with weight 10/Δt. IQIs: `DOF` for P, VX, VY, TEMP and `PROPERTY`
IQIs for Mach number, Reynolds number and heat load.

**Solver (SOL).** `AMESOS_IMPL` direct solve, Newton with `NLA_rel_res_norm_drop` 10⁻⁷,
relaxation 1.0, at most 10 iterations; 30 time steps over ≈ 0.453 with two time levels per DOF
type ([coupling C7](../../../user-guide/interlinks.md#c7-dof-types-fem-msi-sol)). No
preconditioner.

**Output (VIS).** `Heated_Bubble_2D.exo` on the fluid sets with nodal fields
`P, VX, VY, TEMP, Ma, Re, Q`, written every step (`Save_Frequency` 1, `Time_Offset` 10.0).

## What the test checks

The test case calls the workflow with `./Heated_Bubble_2D.so` and requires a return code of 0
(`REQUIRE( tRet == 0 )`). It sets `gInterpolationOrder = 0` and
`gPrintReferenceValues = false`; no field values are compared.

## Run it

The example is not built: `add_subdirectory(Heated_Bubble_2D)` is commented out in
`projects/EXA/fluid/compressible/CMakeLists.txt`, so the ctest below is not registered. After
re-enabling that line, in a build with `MORIS_USE_EXAMPLES=ON`:

```bash
ctest -R Heated_Bubble_2D-1-procs --output-on-failure
```

## Notes

- The IWG on all four sides keeps the name `IWGNitscheOutlets` from the channel deck.
- The Reynolds-number property uses the box height (50) as its length scale.
