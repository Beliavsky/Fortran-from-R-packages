#!/usr/bin/env bash
set -euo pipefail

mode="${1:-check}"
fc="${FC:-gfortran}"
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build/$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

common=(-std=f2018 -Wall -Wextra -Wpedantic)
case "$mode" in
  check) flags=("${common[@]}" -O0 -g -fcheck=all -fbacktrace) ;;
  release) flags=("${common[@]}" -O3) ;;
  *) echo "usage: $0 [check|release]" >&2; exit 2 ;;
esac

sources=(
  pmwr_kinds.f90
  pmwr_types.f90
  pmwr_utils.f90
  pmwr_returns.f90
  pmwr_portfolio.f90
  pmwr_trades.f90
  pmwr_analysis.f90
  pmwr_backtest.f90
  pmwr_identifiers.f90
  pmwr.f90
)

objects=()
for src in "${sources[@]}"; do
  obj="$build/${src%.f90}.o"
  "$fc" "${flags[@]}" -J"$build/mod" -I"$build/mod" -c "$root/src/$src" -o "$obj"
  objects+=("$obj")
done

for source in "$root"/test/*.f90; do
  name="$(basename "$source" .f90)"
  "$fc" "${flags[@]}" -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name="$(basename "$source" .f90)"
  "$fc" "${flags[@]}" -I"$build/mod" "$source" "${objects[@]}" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done

echo "PMwR-fortran $mode build: PASS"
