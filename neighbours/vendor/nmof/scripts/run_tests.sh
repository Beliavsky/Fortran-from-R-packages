#!/usr/bin/env bash
set -euo pipefail

mode=${1:-strict}
root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
build="$root/build-$mode"

case "$mode" in
  strict)
    flags=(-O0 -g -std=f2018 -ffree-line-length-none -Wall -Wextra -Werror \
           -fcheck=all -ffpe-trap=invalid,zero,overflow)
    ;;
  optimized)
    flags=(-O3 -std=f2018 -ffree-line-length-none -Wall -Wextra)
    ;;
  *)
    echo "usage: $0 [strict|optimized]" >&2
    exit 2
    ;;
esac

rm -rf "$build"
mkdir -p "$build/obj" "$build/mod" "$build/bin"

sources=(
  nmof_kinds
  nmof_types
  nmof_rng
  nmof_linalg
  nmof_math
  nmof_qp
  nmof_utilities
  nmof_simulation
  nmof_optimization
  nmof_finance
  nmof_portfolio
  nmof
)

objects=()
for name in "${sources[@]}"; do
  obj="$build/obj/$name.o"
  gfortran "${flags[@]}" -J"$build/mod" -I"$build/mod" \
    -c "$root/src/$name.f90" -o "$obj"
  objects+=("$obj")
done

build_and_run() {
  local src=$1
  local exe="$build/bin/$(basename "${src%.f90}")"
  gfortran "${flags[@]}" -Wno-unused-dummy-argument -I"$build/mod" \
    "$root/$src" "${objects[@]}" -llapack -lblas -o "$exe"
  "$exe"
}

for src in test/test_core.f90 test/test_finance.f90 \
           test/test_optimization.f90 test/test_portfolio.f90; do
  build_and_run "$src"
done

for src in app/demo_nmof.f90 example/finance_example.f90 \
           example/optimization_example.f90; do
  build_and_run "$src"
done

echo "NMOF Fortran $mode build: PASS"
