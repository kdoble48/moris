#!/usr/bin/env bash
# Build and run the vector-of-vectors benchmark (copy / move / move_throwing / csr) under glibc malloc
# and tcmalloc (production MORIS links gperftools tcmalloc).
#
# Usage: run.sh <out_dir> [sizes...]     e.g. run.sh $MORIS_RUNS_DIR/benchmarks/xtk_vector_move 16 32 48
# Env:   BUILD_NINJA  production build.ninja to copy compile DEFINES/INCLUDES from
#        REPEATS      repeats per (variant, size, allocator), default 3
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECTS="$(cd "$HERE/../../../.." && pwd)"
OUT="${1:?usage: run.sh <out_dir> [sizes...]}"
shift
SIZES=(16 32 48)
[ $# -gt 0 ] && SIZES=("$@")
REPEATS="${REPEATS:-3}"
BUILD_NINJA="${BUILD_NINJA:-/home/doble/moris_workspace/moris-develop-local/build_opt/build.ninja}"

mkdir -p "$OUT/bin" "$OUT/patched"

# production defines + includes for the XTK sources, re-pointed at this checkout's projects/
DEFINES=$(grep -m1 -A8 "cl_XTK_Model.cpp.o:" "$BUILD_NINJA" | grep "DEFINES =" | sed 's/^ *DEFINES = //')
INCLUDES=$(grep -m1 -A8 "cl_XTK_Model.cpp.o:" "$BUILD_NINJA" | grep "INCLUDES =" | sed 's/^ *INCLUDES = //' \
    | sed "s#/[^ ]*/projects/#$PROJECTS/#g")
CXXFLAGS="-O3 -DNDEBUG -std=c++17 -m64 -DMORIS_HAVE_PARALLEL"
TCMALLOC_FLAGS="-fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free -ltcmalloc"

python3 "$HERE/patch_header.py" "$PROJECTS/MRS/CNT/src/cl_Vector.hpp" "$OUT/patched"

build() {    # name, extra flags...
    local name=$1
    shift
    # shellcheck disable=SC2086
    g++ $CXXFLAGS $DEFINES "$@" $INCLUDES "$HERE/bench.cpp" -o "$OUT/bin/$name"
}

for alloc in glibc tcmalloc; do
    extra=""
    [ "$alloc" = tcmalloc ] && extra="$TCMALLOC_FLAGS"
    # shellcheck disable=SC2086
    build "copy_$alloc" $extra
    # shellcheck disable=SC2086
    build "move_$alloc" -I"$OUT/patched" -DBENCH_PATCHED $extra
    # shellcheck disable=SC2086
    build "move_throwing_$alloc" -I"$OUT/patched" -DBENCH_PATCHED -DBENCH_MOVE_NOT_NOEXCEPT $extra
    # shellcheck disable=SC2086
    build "csr_$alloc" -DBENCH_CSR $extra
done

CSV="$OUT/results.csv"
echo "allocator,variant,n,tets,edges,build_s,shrink_s,total_s,allocs,peak_live_mb,maxrss_mb,nothrow_move,checksum" > "$CSV"
for n in "${SIZES[@]}"; do
    for alloc in glibc tcmalloc; do
        for variant in copy move move_throwing csr; do
            "$OUT/bin/${variant}_$alloc" "$n" "$REPEATS" | sed "s/^/$alloc,/" >> "$CSV"
        done
    done
done
echo "wrote $CSV"
