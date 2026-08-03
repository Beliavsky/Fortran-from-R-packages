#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
BUILD="$ROOT/build/checked"
if [ -f "$ROOT/backend/lib/libhighs_fortran_bridge.so" ]; then
  export HIGHS_FORTRAN_BRIDGE="$ROOT/backend/lib/libhighs_fortran_bridge.so"
elif [ -f "$ROOT/backend/lib/libhighs_fortran_bridge.dylib" ]; then
  export HIGHS_FORTRAN_BRIDGE="$ROOT/backend/lib/libhighs_fortran_bridge.dylib"
fi
rm -rf "$BUILD"
mkdir -p "$BUILD"
cc -std=c11 -O0 -g -Wall -Wextra -Wpedantic -c "$ROOT/src/highs_dynamic_loader.c" -o "$BUILD/highs_dynamic_loader.o"
cd "$BUILD"
gfortran -std=f2018 -O0 -g -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace -J . -I . -c \
  "$ROOT/src/highs_kinds.f90" "$ROOT/src/highs_constants.f90" "$ROOT/src/highs_sparse.f90" \
  "$ROOT/src/highs_types.f90" "$ROOT/src/highs_c_bindings.f90" "$ROOT/src/highs_solver_api.f90" \
  "$ROOT/src/highs_model_api.f90" "$ROOT/src/highs.f90"
for src in "$ROOT"/test/*.f90 "$ROOT"/example/*.f90 "$ROOT"/app/*.f90; do
  exe="$BUILD/$(basename "${src%.f90}")"
  gfortran -std=f2018 -O0 -g -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace -J . -I . "$src" ./*.o -o "$exe"
  "$exe"
done
