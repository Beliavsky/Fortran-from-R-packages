#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build/validation
mkdir -p build/validation/mod build/validation/obj build/validation/bin
FC=${FC:-gfortran}
MODE=${1:-debug}
if [[ "$MODE" == "optimized" ]]; then
  FLAGS=(-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror -ffree-line-length-132)
else
  FLAGS=(-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace -ffree-line-length-132)
fi
SOURCES=(
  src/yc_kinds.f90 src/yc_types.f90 src/yc_linalg.f90 src/yc_utils.f90
  src/yc_splines.f90 src/yc_models.f90 src/yc_curve_ops.f90
  src/yc_analysis.f90 src/yc_pca_mod.f90 src/yieldcurves.f90
)
for src in "${SOURCES[@]}"; do
  obj="build/validation/obj/$(basename "${src%.f90}").o"
  "$FC" "${FLAGS[@]}" -J build/validation/mod -I build/validation/mod -c "$src" -o "$obj"
done
for src in test/*.f90; do
  exe="build/validation/bin/$(basename "${src%.f90}")"
  "$FC" "${FLAGS[@]}" -J build/validation/mod -I build/validation/mod "$src" build/validation/obj/*.o -o "$exe"
  "$exe"
done
for src in app/*.f90 example/*.f90; do
  exe="build/validation/bin/$(basename "${src%.f90}")"
  "$FC" "${FLAGS[@]}" -J build/validation/mod -I build/validation/mod "$src" build/validation/obj/*.o -o "$exe"
  "$exe" >/dev/null
done
printf 'validation (%s): PASS\n' "$MODE"
