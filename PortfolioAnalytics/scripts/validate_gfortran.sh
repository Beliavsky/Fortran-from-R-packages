#!/usr/bin/env sh
set -eu
MODE=${1:-checked}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build-$MODE"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
if [ "$MODE" = optimized ]; then
  FLAGS="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Wimplicit-interface"
else
  FLAGS="-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wimplicit-interface -fcheck=all -fbacktrace"
fi
SOURCES="pa_kinds pa_types pa_linalg pa_statistics pa_constraints pa_objectives pa_optimizers pa_portfolios pa_views pa_factor_models pa_robust portfolio_analytics"
for name in $SOURCES; do
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/src/$name.f90" -o "$BUILD/obj/$name.o"
done
gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" -c "$ROOT/test/test_support.f90" -o "$BUILD/obj/test_support.o"
for source in "$ROOT"/test/test_*.f90; do
  [ "$(basename "$source")" = test_support.f90 ] && continue
  name=$(basename "$source" .f90)
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$BUILD"/obj/*.o "$source" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
for source in "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $FLAGS -J "$BUILD/mod" -I "$BUILD/mod" "$BUILD"/obj/pa_*.o "$BUILD/obj/portfolio_analytics.o" "$source" -o "$BUILD/bin/$name"
  "$BUILD/bin/$name"
done
