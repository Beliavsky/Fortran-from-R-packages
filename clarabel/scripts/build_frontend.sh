#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
mode=${1:-checked}
case "$mode" in
  checked)
    fflags="-std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace -O0 -g"
    cflags="-std=c11 -Wall -Wextra -Wpedantic -O0 -g"
    ;;
  optimized)
    fflags="-std=f2018 -Wall -Wextra -Wpedantic -O3"
    cflags="-std=c11 -Wall -Wextra -Wpedantic -O3"
    ;;
  *) echo "usage: $0 [checked|optimized]" >&2; exit 2 ;;
esac
build="$root/build/frontend-$mode"
rm -rf "$build"
mkdir -p "$build/mock"
cd "$build"

gcc $cflags -fPIC -shared -I"$root/include" "$root/test/mock_backend.c" \
  -o "$build/mock/libclarabel_fortran_bridge.so"
gcc $cflags -c "$root/src/clarabel_dynamic_loader.c"
for source in \
  clarabel_kinds.f90 clarabel_sparse.f90 clarabel_types.f90 clarabel_psd.f90 \
  clarabel_c_api.f90 clarabel_solver.f90 clarabel.f90; do
  gfortran $fflags -J. -I. -c "$root/src/$source"
done
objects=(clarabel_dynamic_loader.o clarabel_kinds.o clarabel_sparse.o clarabel_types.o \
         clarabel_psd.o clarabel_c_api.o clarabel_solver.o clarabel.o)
export CLARABEL_FORTRAN_BRIDGE="$build/mock/libclarabel_fortran_bridge.so"
for source in "$root"/test/*.f90; do
  exe=$(basename "$source" .f90)
  gfortran $fflags -J. -I. "$source" "${objects[@]}" -o "$exe"
  "./$exe"
done
for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  exe=$(basename "$source" .f90)
  gfortran $fflags -J. -I. "$source" "${objects[@]}" -o "$exe"
done
./unconstrained_qp
./equality_qp
./persistent_update
./demo_clarabel

echo "Frontend/runtime-loader $mode validation passed."
