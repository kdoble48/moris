# Deck API v2

The fork's single-entry-point deck format. It removes the boilerplate of legacy decks — twelve
`extern "C"` module functions, trivial callbacks, seven OPT callbacks — and replaces most
string coupling with typed handles. The design record is {doc}`../internal/DECK_API_RFC`.

Legacy decks keep working unchanged: the loader looks for the new entry symbol first and falls
back to the legacy symbol scan. One shared object must use one style or the other.

## Entry point

```cpp
#include "cl_Input_Deck.hpp"
#include "cl_Input_Deck_Vocabulary.hpp"
#include "fn_FEM_Presets.hpp"

MORIS_DECK( aDeck )
{
    Module_Parameter_Lists& tHmr = aDeck.hmr();   // touching a module activates it
    tHmr.set( "number_of_elements_per_dimension", 45, 15 );
    // aDeck.xtk(), .gen(), .fem(), .sol(), .opt(), .msi(), .vis(), .morisgeneral()
}
```

Each accessor returns the same `Module_Parameter_Lists` a legacy module function receives, so
every legacy `set()` / `add_parameter_list()` call works as before.

## Callbacks without `extern "C"`

Property and field functions are registered in-process and looked up by name, before the
dynamic-symbol fallback:

```cpp
aDeck.register_function( "Func_Neumann_U", &Func_Neumann_U );
```

The registry is type-checked, and the common trivial callbacks are built in: `Func_Const`
(a constant property) and `Output_Criterion` (always write output) need no definition in the
deck.

## Typed set names

`deck::Phase` generates the XTK set-name grammar ([coupling C5](interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)),
so a phase index or side number is written once:

```cpp
const deck::Phase tVoid( 0 );
const deck::Phase tFrame( 2 );

tFrame.bulk();                        // "HMR_dummy_n_p2,HMR_dummy_c_p2"
tFrame.side( deck::Side::Left );      // "SideSet_4_n_p2,SideSet_4_c_p2"
tFrame.ghost();                       // "ghost_p2"
deck::interface( tFrame, tVoid );     // single-sided interface sets
deck::between( tInterior, tFrame );   // double-sided interface set
deck::join( a, b );                   // comma-join set lists
```

`deck::Dofs::Displacement2D` and its siblings supply the DOF strings shared by FEM, MSI, and SOL
([coupling C7](interlinks.md#c7-dof-types-fem-msi-sol)).

## Physics presets

`fem::presets::linear_elastic()` emits a complete linear-elastic block — properties,
constitutive model, bulk IWG, Nitsche Dirichlet, Neumann, ghost, bedding, strain-energy and
volume IQIs — and returns the generated names for further wiring:

```cpp
fem::presets::Linear_Elastic_Config tFrame;
tFrame.mPrefix        = "Frame_";
tFrame.mBulkSets      = tFramePhase.bulk();
tFrame.mDirichletSets = tFramePhase.side( deck::Side::Left );
tFrame.mNeumannSets   = tFramePhase.side( deck::Side::Right );
tFrame.mTraction      = "1.0";
tFrame.mGhostSets     = tFramePhase.ghost();

fem::presets::Linear_Elastic_Names tNames = fem::presets::linear_elastic( tFem, tFrame );
// tNames.mCM, tNames.mYoungsProp, tNames.mStrainEnergyIQI, tNames.mVolumeIQI, …
```

`fem::presets::diffusion()` is the thermal counterpart. Anything a preset does not cover is
added with raw `set()` calls next to it.

## Objective and constraints as expressions

Instead of `compute_objectives`, `compute_constraints`, and their five derivative callbacks,
write the functions over named criteria:

```cpp
aDeck.objective(
        deck::criterion( tNames.mStrainEnergyIQI ) / tInitialStrainEnergy
        + 0.2 * deck::criterion( "IQIPerimeter" ) / tInitialPerimeter );

aDeck.constraint( deck::criterion( tNames.mVolumeIQI ) / tMaxMass - 1.0 <= 0.0 );
```

The values are computed by evaluating the expression tree; the derivatives with respect to the
criteria by reverse-mode differentiation. GEN's `IQI_types` is written automatically in
first-appearance order, and each constraint's type follows from its comparison operator. The
derivatives with respect to the ADVs default to zero, as in nearly every deck, since the ADV
sensitivities flow through the FEM adjoint.

## Worked examples

- `projects/EXA/optimization/LevelSet_Boxbeam_NewIO/` — a two-phase level-set topology
  optimization with GCMMA (616 lines in the legacy twin, 332 here, no OPT callbacks).
- `projects/EXA/optimization/Shape_Sensitivity_Sweep_NewIO/` — a shape-sensitivity sweep that
  passes the same analytic-versus-finite-difference gates as its legacy twin.
- `python3 share/scripts/moris_init.py MyDeck` — scaffolds a new deck from a run-verified
  template.
