#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
MODE=${1:-checked}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-manual-$MODE"

case "$MODE" in
  checked)
    FLAGS="-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"
    ;;
  optimized)
    FLAGS="-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -O3"
    ;;
  *)
    echo "usage: $0 [checked|optimized]" >&2
    exit 2
    ;;
esac

rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"

compile_module() {
  src=$1
  obj=$2
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/$src" -o "$BUILD/obj/$obj"
}

compile_module src/ppcor_kinds.f90 ppcor_kinds.o
compile_module src/ppcor_special.f90 ppcor_special.o
compile_module src/ppcor_linalg.f90 ppcor_linalg.o
compile_module src/ppcor_stats.f90 ppcor_stats.o
compile_module src/ppcor.f90 ppcor.o

OBJS="$BUILD/obj/ppcor_kinds.o $BUILD/obj/ppcor_special.o $BUILD/obj/ppcor_linalg.o $BUILD/obj/ppcor_stats.o $BUILD/obj/ppcor.o"

for src in "$ROOT"/test/*.f90; do
  name=$(basename "$src" .f90)
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" $OBJS "$src" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done

for src in "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
  name=$(basename "$src" .f90)
  "$FC" $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" $OBJS "$src" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name" >/dev/null
done

echo "All $MODE tests and examples passed."
