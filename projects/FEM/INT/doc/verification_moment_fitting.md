# Design & verification: moment-fitted cut-cell quadrature

Record for the feature-flagged moment-fitting rule (`use_moment_fitting`, default OFF).
Two phases: an offline Phase-0 design study that **rejected** the naive design point, and the
in-MORIS Phase-1 implementation verified through five gates, including bitwise flag-off
identity and a GCMMA optimization A/B.

## Phase 0 — why the design is shaped this way (offline study, 2026-07-14)

An offline harness fitted weights at *fixed parent-cell tensor-Gauss points* (27/64 pts,
total degree d ≤ 5) against moments from 3,844 real XTK cut cells (131k material tets,
volume fractions down to 5e-15) and measured element-stiffness consistency against a
degree-11 reference rule. Verdict: **no-go at that design point** —

- unconstrained least squares reaches K_e consistency 2.5e-3 (median) but produces negative
  weights in 92% of cut cells and destroys elementwise semi-definiteness in 79%;
- non-negative fits (NNLS) are exactly SPD-safe but plateau at K_e error ~2e-2: with points
  fixed on the full parent cell, non-negativity cannot match small-cut moments.

The blocker is **point placement**, not moment fitting itself. The Phase-1 design responds:

1. **Candidate points = the cluster's own tessellation Gauss points** (already inside the
   material region), pruned by the NNLS active set;
2. **Lawson–Hanson NNLS** (own implementation: MGS-QR passive solves, stagnation guard) —
   non-negative weights ⇒ machine-PSD element matrices by construction;
3. **Bounding-box-scaled fitting basis** (Legendre products in coordinates scaled to the
   cluster's material bounding box) — thin slivers otherwise make basis columns
   near-collinear; the scaled frame reaches ~1e-15 moment residual at α = 1e-3 where the raw
   frame plateaued at 1e-8;
4. **Per-cluster fallback** to the tessellated rule on residual > tolerance, inverted tet,
   or non-TRI/TET integration cells.

## Seam

A cut BULK cluster evaluates ONE host element (largest-|detJ| tet) under a temporary
per-cluster integration-rule override on `fem::Set`; fitted points are mapped to the host
tet's parametric frame and weights divided by its parametric detJ, so the standard
`w · detJ` product in the element loop reproduces the fitted physical rule exactly (affine
IP maps = HMR lattices). Fitted paths: residual, Jacobian, QI, dQIdu. **Deliberately
tessellated:** dRdp/dQIdp (weights are not re-fitted under FD geometry perturbations), VIS
output, cluster measures, ghost/side/double-side sets.

## Phase 1 — verification gates (all measured, single rank, sliver deck n=6 / p=2)

| gate | check | result |
|---|---|---|
| identity | flag-off Jacobian export vs pre-change build | **bitwise identical** (run-to-run determinism established first) |
| 0 | in-MORIS rule dump vs the offline harness math, α ∈ {0.5, 0.05, 1e-3, +ghost} | moment residual ~1e-15, volume ~1e-15, machine-PSD; 1120 candidates → ~35 retained (d=4) |
| 1 | VOLUME criterion fitted vs tessellated | identical to all printed digits at every α |
| 2a | assembled K export | fitted K Cholesky-SPD at all α; ‖K_f−K_t‖/‖K_t‖ = 4.7e-5 → 9.1e-8 (α 0.5 → 1e-3); ghost-ON κ₂ ratio 1.000 |
| 2b | loaded probe (interface traction) | strain energy identical to 6 digits; SE-density L2 diff 5.6e-3 @ d=4 → 4.9e-5 @ d=5 → 1.1e-5 @ d=6 |
| 3 | boxbeam GCMMA A/B, 7 aligned evaluations | objective trajectory rel. diff ≤ 4.05e-9 per iteration; FD sensitivities + remeshing + .so reload ran through with fitting active |

Overhead at the probe size: fitted end-to-end 0.83 s vs 1.71 s tessellated (**2.1×
faster**); weights re-fit automatically per design iteration when FEM rebuilds.

## Full-length production A/B (3D box-beam, 151 GCMMA iterations per arm)

Beyond the 7-evaluation gate, two complete optimizations of the 3D frameless box-beam were
run to 151 iterations each (~2h45m tessellated, ~2h53m fitted, both exit 0; 2026-07-14,
`phase1_3d_production/`):

- objective trajectories: median per-iteration rel. difference **8.6e-7**, p90 6.8e-4,
  max 3.0e-3 (a brief transient at iteration 110); final objectives 0.267846 vs 0.267822
  (rel. diff **9.1e-5**);
- final designs: section-for-section visually identical (figure below), with identical
  measured topology — 1 material component, 3 void components, 0 enclosed cavities;
  volume fraction 0.202238 (tess) vs 0.202249 (fit).

Design evolution per arm (isosurface renders every 10 evaluations, with per-eval objective /
strain energy and the objective trajectory):

![Tessellated-rule arm: design evolution over 151 GCMMA
evaluations](verification_moment_fitting_evolution_tess.png)

![Moment-fitted arm: design evolution over 151 GCMMA
evaluations](verification_moment_fitting_evolution_fit.png)

Pixel-level identity check of the final designs (section slices):

![Final designs after 151 GCMMA iterations: tessellated rule (top) vs moment-fitted rule
(bottom), four cross-sections and four elevations each](verification_moment_fitting_ab_designs.png)

## Usage guidance

- Deck keys (FEM computation list): `use_moment_fitting` (default false),
  `moment_fitting_degree` (default 4), `moment_fitting_tolerance` (default 1e-10).
- **Degree lever:** d=4 gives K_e consistency ~5e-3 (the stiffness integrand lies outside
  the d=4 space); accuracy-critical work should use **d ≥ 5** (K_e 4.3e-4 @ 56 pts;
  1.7e-4 @ d=6/84 pts) — Gate 2b shows d ≥ 5 sits below the ~1e-3 discretization scale.
- Applies to non-trivial cut BULK clusters with TRI/TET tessellation only; everything else
  keeps the tessellated rule automatically.
- Flag-off isolation is structural, not incidental: production paths never call the new
  `solve_least_squares`; it is exercised only inside the gated fitting path.

## Unit tests (`UT_FEM_Moment_Fitting.cpp`, tag `[MomentFitting]`)

Basis evaluation, NNLS solver, TET-rule degree map, the uncut-hex rigor gate (fitted
weights recover tensor-Gauss to ~1e-15), slab cells at α ∈ {0.5, 0.05, 1e-3}, and the
host-tet frame conversion.

## Provenance

Workspace study `studies/moment_fitting/` (harness + mflib) and
`runs/studies/moment_fitting/`: `FINDINGS.md` (Phase-0 tables + figures),
`PHASE1_LOG.md` (gate protocols, dirs `phase1/identity_*`, `gate3_tess/`, `gate3_fit*/`),
`results.json`, `figures/fig_*.png`. Rule-dump escape hatch:
`MORIS_MOMENT_FITTING_DUMP=<path>` appends per-cluster fitted points/weights.
