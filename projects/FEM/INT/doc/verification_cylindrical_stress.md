# Verification: IQI_Cylindrical_Stress (`IQI_Type::CYLINDRICAL_STRESS`)

Verification record for the cylindrical-stress IQI, following the staged method
(consistency limit → signature/exactness → literature benchmark). Data provenance at the end.

## What the IQI computes

At an evaluation point **x** with polar angle θ about a deck-configured axis origin
(`function_parameters = "cx;cy"`), the IQI rotates the elastic CM's Cauchy stress into
cylindrical components, selected by `vectorial_field_index`:

| index | component |
|---|---|
| 0 | σ_rr = σ_xx c² + σ_yy s² + 2 σ_xy c s |
| 1 | σ_θθ = σ_xx s² + σ_yy c² − 2 σ_xy c s |
| 2 | σ_rθ = (σ_yy − σ_xx) c s + σ_xy (c² − s²) |

with c = cos θ, s = sin θ.

## Tier 1 — unit-level exactness (`UT_FEM_IQI_Cylindrical_Stress`, tag `[IQI_Cylindrical_Stress]`)

Plane-stress `STRUC_LIN_ISO` CM (E = 1, ν = 0.3) on a unit-square QUAD with a displacement
state chosen so the stress at the evaluation point has σ_xy ≠ 0 and σ_xx ≠ σ_yy (both
REQUIREd — a degenerate state would let wrong rotations pass). Checks, all to 1e-10:

1. **Axis-aligned reduction, θ = 0** (origin west of the point): σ_rr = σ_xx, σ_θθ = σ_yy,
   σ_rθ = σ_xy.
2. **Axis-aligned reduction, θ = 90°** (origin south): σ_rr = σ_yy, σ_θθ = σ_xx,
   σ_rθ = −σ_xy.
3. **Oblique angle** (θ = 45°): all three components against the closed-form rotation
   formulas evaluated independently in the test.
4. **Trace invariance**: σ_rr + σ_θθ = σ_xx + σ_yy at the oblique angle (rotation must
   preserve the first invariant).

Status: passes (15 assertions). The full `[IQI]` slice (154 cases, 7609 assertions) passes
with the IQI registered, i.e. no factory/enum regressions.

## Tier 2 — literature benchmark: Kirsch plate-with-hole

Uniaxial tension of an infinite plate with a circular hole (Kirsch 1898): analytic hoop
stress on the hole boundary σ_θθ(θ)/σ₀ = 1 − 2cos 2θ, giving the stress concentration
**σ_θθ/σ₀ = 3 at the crown** and −1 at the sides. Setup: immersed (XTK-cut) hole, R = 0.2,
centre (0.5, 0.5), E = 1, ν = 0.3, exact-solution displacement BCs on the outer boundary;
the IQI evaluates σ_rr/σ_θθ/σ_rθ about the hole centre as VIS nodal fields.

**Crown stress vs analytic 3.0** (p = 1 ladder, n = 24…80): 2.908, 2.965, 3.013, 2.989,
2.932 — within 3.1% of the analytic concentration at every resolution.

**h-refinement of the interface hoop-stress error** |σ_θθ,FE − σ_θθ,exact|/σ₀ on r = R,
and the displacement L2 error (fitted rates over the ladder):

| p | hoop rms error (coarsest → finest) | rate | disp L2 rate |
|---|---|---|---|
| 1 | 2.07e-01 → 8.60e-02 | O(h^0.88) | O(h^2.24) |
| 2 | 4.01e-02 → 7.25e-03 | O(h^1.72) | O(h^2.76) |
| 3 | 2.80e-02 → 4.70e-03 | O(h^2.48) | O(h^3.37) |

The hoop-stress error converges at ≈ O(h^p) — the expected one-order-below-displacement
rate for a stress (gradient) quantity on a cut interface — and improves systematically with
p. This is a benchmark of the *whole chain* (immersed discretization + CM stress + this
IQI's rotation); the rotation itself is exact per Tier 1.

## Provenance

- Unit test: `projects/FEM/INT/test/UT_FEM_IQI_Cylindrical_Stress.cpp` (this PR).
- Benchmark decks (workspace, PyMORIS): `studies/moris_immersed_pathologies/decks/
  kirsch_validation.py` (immersed/enriched) and `kirsch_ersatz.py` (ersatz comparison),
  both consuming `IQI_Type::CYLINDRICAL_STRESS` with `function_parameters = "cx;cy"` and
  `vectorial_field_index` 0/1/2.
- Ladder data: `runs/studies/moris_immersed_pathologies/kirsch_refine/kirsch_refine.json`
  (p = 1, crown values) and `kirsch_order_refine/kirsch_order_refine.json` (p = 1/2/3);
  analysis scripts `studies/moris_immersed_pathologies/analysis/kirsch_refine.py`,
  `kirsch_order_refine.py` (2026-07).
