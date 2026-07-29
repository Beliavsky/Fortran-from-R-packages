#!/usr/bin/env bash
set -euo pipefail

mode="${1:-strict}"
fc="${FC:-gfortran}"
root="$(cd "$(dirname "$0")" && pwd)"
build="$root/build-${mode}"
moddir="$build/mod"
objdir="$build/obj"
bindir="$build/bin"

case "$mode" in
  strict)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all \
           -ffpe-trap=invalid,zero,overflow -fbacktrace)
    ;;
  optimized)
    flags=(-std=f2018 -O3 -Wall -Wextra -Werror)
    ;;
  *)
    echo "Usage: $0 [strict|optimized]" >&2
    exit 2
    ;;
esac

rm -rf "$build"
mkdir -p "$moddir" "$objdir" "$bindir"

sources=(
  frapo_kinds.f90
  frapo_types.f90
  frapo_linalg.f90
  frapo_statistics.f90
  frapo_series.f90
  frapo_risk.f90
  frapo_optimization.f90
  frapo_portfolios.f90
  frapo.f90
)

objects=()
for source in "${sources[@]}"; do
  object="$objdir/${source%.f90}.o"
  "$fc" "${flags[@]}" -J "$moddir" -I "$moddir" \
    -c "$root/src/$source" -o "$object"
  objects+=("$object")
done

for source in "$root"/test/*.f90; do
  name="$(basename "${source%.f90}")"
  "$fc" "${flags[@]}" -I "$moddir" "$source" "${objects[@]}" \
    -llapack -lblas -o "$bindir/$name"
  "$bindir/$name"
done

for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  name="$(basename "${source%.f90}")"
  "$fc" "${flags[@]}" -I "$moddir" "$source" "${objects[@]}" \
    -llapack -lblas -o "$bindir/$name"
done

echo "All $mode builds and tests passed."
