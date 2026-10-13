# Running an input deck

A MORIS input deck is a C++ source file. It is compiled into a shared object (`.so`) that the
`moris` executable loads at run time and queries for its parameter lists and user callbacks.
(The executable also accepts an `.xml` parameter file for decks without callbacks.)

## Compile and run

From the directory containing `MyDeck.cpp`:

```bash
$MORISROOT/share/scripts/create_shared_object.sh . build_opt MyDeck
$MORISROOT/build_opt/projects/mains/moris ./MyDeck.so
```

`create_shared_object.sh` compiles the deck against an existing build tree (which must have been
configured with `BUILD_ALL=ON`). For a debug run, pass `build_dbg` and use that tree's
executable. The wrapper `share/scripts/moris.sh` does both steps:

```bash
moris.sh opt 1 MyDeck      # optimized build, 1 process
moris.sh dbg 4 MyDeck      # debug build, 4 MPI ranks
```

It reads the build locations from `MORISBUILDOPT` and `MORISBUILDDBG`.

## Parallel runs

```bash
mpirun -np 4 $MORISROOT/build_opt/projects/mains/moris ./MyDeck.so
```

Use a parallel linear solver for parallel runs. The default Amesos direct solver is serial; set
`Solver_Type` to `Amesos_Superludist` or use an iterative solver (see {doc}`../user-guide/modules`).

## Command-line options

| Option | Meaning |
|---|---|
| `--outputlog FILE`, `-ol FILE` | Write the log to `FILE` (one file per rank in parallel runs) |
| `--severity N`, `-sl N` | Log severity threshold |
| `--directoutput N`, `-do N` | Structured console output level |
| `--meshgen`, `-mg` | Only generate and write the foreground and background meshes |
| `--perf-level N` | Performance-report granularity (0 = run total … 3 = full tree, < 0 off) |
| `--perf-report FILE` | Path of the JSON performance report |
| `--pause`, `-p` | Print process IDs and wait for input (attach a debugger) |
| `--version`, `-v` | Print the version and exit |
| `--help`, `-h` | List the options |

The deck file (`.so` or `.xml`) is the last argument.

## What a run writes

- **Exodus files** (`.exo`) from the VIS module: the solution fields on the cut mesh, one time
  step per optimization iteration. Open them in ParaView.
- **`xtk_temp.exo`** and the GEN output mesh, when requested: the integration mesh and the
  level-set fields.
- **`ADV_Alg_0_Iter_N.hdf5`**: the design variables of every optimization iteration. Any of
  these files restarts a run (`restart_file` in the OPT parameter list).
- **`Parameter_Receipt.xml`**: every parameter the run actually used, including defaults. Diff
  two receipts to find why two runs differ.

## Starting from an existing deck

The quickest way to a working deck is to adapt an example from `projects/EXA/` (optimization,
structure, thermal, fluid), or to scaffold one in the {doc}`Deck API v2 <../user-guide/deck-api-v2>`
style:

```bash
python3 share/scripts/moris_init.py MyDeck --dir my_case
```

This emits a deck that compiles and runs as written: an immersed linear-elastic plate with a
finite-difference sensitivity check. Adapt it from there.
