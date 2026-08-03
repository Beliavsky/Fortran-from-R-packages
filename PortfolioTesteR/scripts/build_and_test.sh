#!/usr/bin/env sh
set -eu

mode=${1:-checked}
fc=${FC:-gfortran}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build/obj" "$build/bin"

common="-std=f2018 -Wall -Wextra -Werror -ffree-line-length-none -J$build/obj -I$build/obj"
case "$mode" in
  checked) flags="$common -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow" ;;
  release) flags="$common -O3" ;;
  *) echo "usage: $0 [checked|release]" >&2; exit 2 ;;
esac

sources="
ptr_kinds
ptr_types
ptr_utils
ptr_data
ptr_indicators
ptr_filters
ptr_weighting
ptr_performance
ptr_backtest
ptr_cross_sectional
ptr_ml
ptr_strategies
ptr_optimization
ptr_walk_forward
portfolio_tester
"
objects=""
for name in $sources; do
  "$fc" $flags -c "$root/src/$name.f90" -o "$build/obj/$name.o"
  objects="$objects $build/obj/$name.o"
done

run_group() {
  directory=$1
  for source in "$root/$directory"/*.f90; do
    [ -f "$source" ] || continue
    name=$(basename "$source" .f90)
    exe="$build/bin/$name"
    "$fc" $flags "$source" $objects -o "$exe"
    "$exe"
  done
}

run_group test
run_group example
run_group app
printf '%s\n' "PortfolioTesteR $mode build: PASS"
