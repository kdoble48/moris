# Coding style

Formatting is defined by `.clang-format` at the repository root; run `clang-format -i` on the
files you change. The conventions below cover what the formatter cannot.

## Variables

Prefix every variable with its kind:

| Prefix | Kind | Example |
|---|---|---|
| `m` | member | `mNumRows` |
| `a` | function argument | `aElementInd` |
| `t` | local (temporary) | `tTotalSum` |
| `s` | static constant in an input deck | `sUseGhost` |

Use camelCase after the prefix, with qualifiers on the left (`tUserId`, not `tIdUser`).
Approved abbreviations: `Glb` (global), `Loc` (local), `Id` (1-based identifier),
`Ind` (0-based index), `Num` (number of).

- Declare and initialize one variable per statement.
- Attach `*` and `&` to the type: `uint& aCount`.
- Write real literals with a decimal point (`5.0`, not `5`).
- Spell out types; use `auto` only where the type is not expressible, such as expression
  templates or third-party return types.

## Classes and functions

- Class names are capitalized words joined by underscores: `Cut_Integration_Mesh`.
- Members are ordered `public`, `protected`, `private`.
- Function names are lower case with underscores: `compute_jacobian`.
- Separate function definitions with a blank line, a line of dashes, and another blank line.

## Files

One class per file, named after it with a `cl_` prefix: `cl_FEM_IWG_Struc_Linear_Bulk.hpp`.
Free-function files use `fn_`, as in `fn_PRM_FEM_Parameters.hpp`. Every file starts with the
University of Colorado MIT-license header.

## Comments

Comment the intent of each block, with references for any formulation taken from the
literature:

```cpp
// Compute first- and second-order ghost-penalty contributions using a
// displacement-gradient formulation. Cf. Burman and Hansbo, "Fictitious domain
// methods using cut elements: III. A stabilized Nitsche method for Stokes' problem"
this->compute_ghost_pen_disp_grad();
```

Document classes and public functions with Doxygen comments, so they appear in the
{doc}`API reference <../api/index>`; see {doc}`documentation`.

## Errors and assertions

- `MORIS_ERROR( condition, format, … )` checks a condition in every build; use it for invalid
  user input and states the program cannot recover from.
- `MORIS_ASSERT( condition, format, … )` is compiled only with `MORIS_HAVE_DEBUG`; use it for
  internal invariants in hot code.

Write messages that name the function and the offending value, so the log line alone locates the
problem.
