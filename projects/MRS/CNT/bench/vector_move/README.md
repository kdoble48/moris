# moris::Vector move-semantics benchmark

`moris::Vector` (`MRS/CNT/src/cl_Vector.hpp`) declares a copy constructor and no move constructor, so the
compiler does not generate a move. When an outer `Vector< Vector< T > >` reallocates (growth by
`push_back( value, increment )`, or the outer `shrink_to_fit` in `shrink_to_fit_all`), `std::vector` therefore
copies every inner list: one new allocation plus a copy per row, then frees the old rows.

`bench.cpp` replays the container pattern of `xtk::Integration_Mesh_Generator::create_edges_from_element_to_node()`
(capacity policy, inner reserves, `edge_exists()` lookup, final `shrink_to_fit_all`) on a structured tet mesh,
using the real `cl_Vector.hpp`:

| variant | header | what it shows |
|---|---|---|
| `copy` | as in MORIS | today's behaviour |
| `move` | `patch_header.py` adds a `noexcept` move ctor + move assignment | the proposed change |
| `move_throwing` | same, without `noexcept` | `std::vector` only moves `noexcept` elements on reallocation; otherwise it still copies |
| `csr` | — | the same connectivity as flat arrays (fixed width + CSR), for comparison |

## Run

```bash
./run.sh $MORIS_RUNS_DIR/benchmarks/xtk_vector_move/<date> 16 32 48 64
```

Compile defines and include paths are copied from the production `build.ninja` (`BUILD_NINJA`), flags match the
production release build (`-O3 -std=c++17`), and every variant is built against glibc malloc and tcmalloc. Output:
`results.csv` (one line per repeat: phase times, allocation count, peak live bytes, max RSS, `noexcept`-move flag,
an order-independent checksum that must agree across variants).
