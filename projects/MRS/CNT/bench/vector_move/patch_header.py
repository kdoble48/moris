#!/usr/bin/env python3
"""Derive a move-enabled copy of MORIS's cl_Vector.hpp for the vector-of-vectors benchmark.

Usage: patch_header.py <path/to/cl_Vector.hpp> <out_dir>

Inserts, right after the existing copy constructor and copy assignment:
    Vector( Vector< T >&& aVector ) BENCH_MOVE_NOEXCEPT;
    Vector< T >& operator=( Vector< T >&& aVector ) BENCH_MOVE_NOEXCEPT;
BENCH_MOVE_NOEXCEPT is `noexcept` unless the benchmark is compiled with -DBENCH_MOVE_NOT_NOEXCEPT,
which shows that std::vector falls back to copying elements whose move constructor may throw.
Everything else in the header is left byte-identical, so the benchmark tracks the real container.
"""
import pathlib
import sys

COPY_CTOR = """        Vector( const Vector< T >& aVector )
        {
            mVector = aVector.data();
        }
"""

MOVE_CTOR = """
        //------------------------------------------------------------------
        // BENCHMARK PATCH: move constructor

        Vector( Vector< T >&& aVector ) BENCH_MOVE_NOEXCEPT
                : mVector( std::move( aVector.mVector ) )
        {
        }
"""

COPY_ASSIGN = """        operator=( Vector< T > const & aVector )
        {
            mVector = aVector.data();

            return *this;
        }
"""

MOVE_ASSIGN = """
        //------------------------------------------------------------------
        // BENCHMARK PATCH: move assignment

        Vector< T >&
        operator=( Vector< T >&& aVector ) BENCH_MOVE_NOEXCEPT
        {
            mVector = std::move( aVector.mVector );

            return *this;
        }
"""

PRELUDE = """
#ifdef BENCH_MOVE_NOT_NOEXCEPT
#define BENCH_MOVE_NOEXCEPT
#else
#define BENCH_MOVE_NOEXCEPT noexcept
#endif
"""


def main() -> None:
    src = pathlib.Path(sys.argv[1]).read_text()
    out_dir = pathlib.Path(sys.argv[2])
    for anchor in (COPY_CTOR, COPY_ASSIGN):
        if src.count(anchor) != 1:
            sys.exit(f"anchor not found exactly once in header:\n{anchor}")
    src = src.replace(COPY_CTOR, COPY_CTOR + MOVE_CTOR)
    src = src.replace(COPY_ASSIGN, COPY_ASSIGN + MOVE_ASSIGN)
    # the macro must exist before the class; put it after the include guard's #define
    guard = "#define MORIS_CONTAINERS_CL_VECTOR_HPP_\n"
    if src.count(guard) != 1:
        sys.exit("include guard not found")
    src = src.replace(guard, guard + PRELUDE)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "cl_Vector.hpp").write_text(src)


if __name__ == "__main__":
    main()
