# RFC: Input-Deck Format Redesign ("Deck API v2")

Status: **accepted** (2026-07-17) — staged implementation in progress on `feature/deck-api-v2`.
Scope: hand-written C++ input decks. The 78 EXA decks are regression tests and keep
compiling/running unchanged throughout (dual-path support); the PyMORIS generator keeps emitting
the legacy format.

## 1. Problem — measured

Survey across the 78 decks under `projects/EXA/` (2026-07):

| Boilerplate item | Occurrence | Notes |
|---|---|---|
| `Func_Const` (4-line constant-property callback) | 67/78, byte-identical | already redundant: `fem::Property` defaults to a constant value function (`cl_FEM_Property.hpp:84`) |
| `Output_Criterion` (`return true;`) | 78/78 | pure ceremony |
| 7-function OPT callback set (`compute_objectives`, `compute_constraints`, `compute_d*_d*`, `get_constraint_types`) | 40/78 | includes pure-analysis decks (e.g. `HeatConduction.cpp:312-393`); the two `*_dadv` are ~always zero-returns |
| 9 `<MODULE>ParameterList` wrappers | 77/78 | MSI/MORISGENERAL usually empty stubs; XTK bodies ~78% identical across decks |
| FEM property→CM→SP→IWG→IQI block | every deck | ~150 lines recurring per physics family (`Mach_Leading_Edge.cpp`: 33 properties → one `Func_Const`, ~130 lines of constant-binding) |

Typical deck: 295–2632 lines, mean ~742; roughly 30–35% trivial-or-verbatim. Every cross-module
link is an unchecked string (`value_function = "tYoungsFunc"`, GEN `IQI_types` order = OPT
criteria index, XTK set-name grammar, FEM⇄MSI⇄SOL DOF agreement): typos compile and silently
misbehave. There is no scaffolding — new decks start as a copy of an existing example.
`projects/EXA/include/EXA_Common.hpp` (built-in `Func_Const`/`Output_Criterion`/
`get_constraint_types`) exists but has zero adopters.

## 2. Loading contract today (verified)

- One `.so` per deck; `Library_IO` dlopens it (`cl_Library_IO.cpp:207`, `RTLD_NOW`) and dlsyms up
  to 12 symbols `<MODULE>ParameterList( Module_Parameter_Lists& )`.
- **Missing symbol ⇒ module cleared; present-but-empty function ⇒ fn_PRM defaults kept**
  (`cl_Library_IO.cpp:321-324`). Load-bearing (omitting a symbol disables a module) but a trap.
- Every user callback resolves through one choke point, `Library_IO::load_function<T>`
  (`cl_Library_IO.hpp:231-264`). FEM/GEN/HMR/SOL callbacks are string-named parameters; the 7 OPT
  callbacks are fixed names, and OPT **re-dlopens the same .so** via its `"library"` parameter
  (`cl_OPT_Problem_User_Defined.cpp:26-39`).
- The XML input path is module-granular, not key-granular: a module absent from the XML file is
  **reset to defaults during the XML pass**, discarding `.so`-provided values
  (`cl_Library_IO.cpp:363-374`), and an absent submodule likewise (`:401-408`). Additionally, a
  module that IS in the XML gets its parsed lists **appended behind a fresh ctor-seeded default**
  (the XML value lands at index 1; consumers reading index 0 never see it). Both pinned in
  `UT_IOS_Library_IO_Deck_Semantics.cpp`.
- Decks compile standalone (`g++ -shared -fPIC`, nothing linked); symbols resolve at dlopen
  against the `-rdynamic` host binary — so **core-provided free functions are callable from
  decks**. Caveat: moris libs are static archives, so a TU referenced only by decks is dropped
  from the executable at link time (presets need a force-link anchor).
- `fem::PropertyFunc` is already a `std::function`; `Problem_User_Defined` already has a
  function-pointer constructor (`cl_OPT_Problem_User_Defined.cpp:44-63`).

## 3. Design

### Principles

1. **Dual path.** The loader dlsyms a new single entry symbol first (`MORISInputDeck`); if absent,
   the legacy 12-symbol scan runs byte-for-byte. Both styles in one `.so` is a hard error.
2. **Lowering layer, not a parallel schema.** The new `Deck` builder writes into today's
   `Vector<Module_Parameter_Lists>` plus a `Function_Registry` (name → type-erased
   `std::function`). Receipt/XML/validation machinery unchanged; registered lambdas get stable
   generated names so the receipt stays meaningful.
3. **Registry-first resolution** at the `load_function` choke point: registry → dlsym → builtin
   defaults. Legacy decks hit dlsym unchanged; builtins make the trivial callbacks optional even
   for legacy decks.
4. **Presets emit through public builder calls** and return handle bundles. Escape hatches at
   every tier: preset → builder → per-entry raw `.set()` → whole-module raw access.
5. **Typed handles make string coupling impossible by construction**: `Phase`/`Region`/`Side`
   generate the XTK set-name grammar; criteria order is assigned by the objective-expression
   collector; DOF lists thread FEM→MSI/SOL; VIS fields come from IQI handles.

### Target authoring experience

One entry point, zero free functions: `MORIS_DECK( Deck& d )` +
`d.mesh()/.geometry()/.fem()/.optimization()/.solver()/.output()` sections, real C++ lambdas for
spatial properties, `presets::IsotropicElasticity{...}` for the standard elastic block,
`opt.objective( tSE/tSE0 + 0.2*tPerim/tP0 )` replacing the 7 OPT callbacks (small-tree
reverse-mode for `dObjective/dCriteria`). Projected: `Levelset_Boxbeam.cpp` 616 → ~105 lines;
`Shape_Sensitivity_Sweep` 416 → ~180; a minimal analysis deck < 120.

## 4. Stages

| Stage | Content | Key files |
|---|---|---|
| 0 | Safety-net UTs pinning current loader semantics; this RFC | `projects/MRS/IOS/test/UT_IOS_Library_IO_Deck_Semantics.cpp` (+ fixture) |
| 1 | Builtin fallbacks (`Func_Const`, `Output_Criterion`); OPT `*_dadv` optional (zero default, loud log); `get_constraint_types` → `constraint_types` parameter; fence the missing-symbol trap | `fn_Library_Builtin_Functions.*`, `cl_OPT_Problem_User_Defined.cpp`, `cl_Library_IO.cpp` |
| 2 | `MORISInputDeck` entry point; `Function_Registry`; `Input_Deck` core (touch-to-activate modules); OPT double-dlopen fix (library plumb) | `cl_Function_Registry.*`, `cl_Input_Deck.*`, `cl_Library_IO.*`, `cl_OPT_Manager.*` |
| 3 | Typed layer: handles, vocabulary (`Phase`/`Region`/`Side`/`Dofs`), typed registry entries, OPT expression API | new `projects/DCK/` |
| 4 | C++ physics presets (independent; usable from legacy decks) | `projects/FEM/INT/src/CORE/fn_FEM_Presets.*` |
| 5 | Load-time interlink validation (C1-C8) + `moris --validate` dry-run with did-you-mean; warning-first rollout (`MORIS_STRICT_INPUT=1`) | `fn_DCK_Validate.*`, `fn_WRK_Workflow_Main_Interface.cpp` |
| 6 | `moris_runner.py init/validate/new-from`; `_NewIO` EXA twin ports (same reference values); knowledge-base updates | `share/scripts/moris_runner.py`, `projects/EXA/` |

Sequencing: 0 → 1 → 2 → 3; 4 independent (after 0); 5 after 2; 6 last. Gate per stage: full
build + `ctest` including all EXA decks, no new failures vs baseline.

## 5. Risks

- **Silent zero gradients** if `*_dadv` defaults are wrong for a deck that truly needs them:
  mitigated by loud log + only defaulting the pair that is zero in every surveyed deck.
- **Static-archive symbol dropout** for core-hosted presets: force-link anchor in WRK
  (dependency-direction safe); fallback header-inline.
- **C++ types across dlopen**: already the shipping contract (`Module_Parameter_Lists&`);
  same-toolchain requirement documented.
- **C++ standard**: preset config structs use designated initializers (C++20 in the deck TU only;
  core stays C++17 until a measured decision).
- Unverified at design time (checked during implementation): `get_constraint_types` call timing;
  SOL behavior on empty `TSA_Output_Criteria`; GEN field-function `std::function` compatibility.
