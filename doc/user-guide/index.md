# User guide

This guide is for people who write and run MORIS input decks.

- {doc}`deck-anatomy` — what a deck is, how MORIS loads it, and the two authoring styles.
- {doc}`modules` — the parameter lists of each module (HMR, XTK, GEN, FEM, SOL, OPT, MSI,
  VIS) and the types and keys each one exposes.
- {doc}`interlinks` — the eight string couplings between modules. Most deck bugs are broken
  couplings, and most of them fail silently.
- {doc}`deck-api-v2` — the single-entry-point deck style: typed set names, physics presets,
  and objective/constraint expressions with automatic criteria gradients.
- {doc}`workflow` — the recommended path from a new deck to a trustworthy optimization run:
  forward check, short run, sensitivity sweep, production run.

## Sources of truth

This guide summarizes. When it disagrees with the code, the code wins:

1. `projects/PRM/src/fn_PRM_<MODULE>_Parameters.hpp` — every parameter key and its default.
2. `projects/PRM/ENM/src/cl_<MODULE>_Enums.hpp` and `GEN_Data_Types.hpp` — the valid values
   of every enumerated parameter.
3. The factories (`projects/FEM/INT/src/<SEC>/cl_FEM_<SEC>_Factory.cpp`,
   `projects/OPT/src/fn_OPT_create_algorithm.cpp`) — proof that a type is actually wired.
4. The example decks in `projects/EXA/` — working configurations to adapt.
