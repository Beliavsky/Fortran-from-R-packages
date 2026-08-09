#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
rm -rf build_strict
mkdir -p build_strict/mod build_strict/obj build_strict/bin
flags=(-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0)
src=(src/coneproj_kinds.f90 src/coneproj_types.f90 src/coneproj_linalg.f90 src/coneproj_stats.f90 \
     src/coneproj_core.f90 src/coneproj_shape.f90 src/coneproj_regression.f90 src/coneproj.f90)
gfortran "${flags[@]}" -Jbuild_strict/mod -Ibuild_strict/mod -c "${src[@]}"
mv ./*.o build_strict/obj/
for f in test/*.f90; do
  exe="build_strict/bin/$(basename "$f" .f90)"
  gfortran "${flags[@]}" -Ibuild_strict/mod "$f" build_strict/obj/*.o -o "$exe"
  "$exe"
done
for f in example/*.f90; do
  exe="build_strict/bin/$(basename "$f" .f90)"
  gfortran "${flags[@]}" -Ibuild_strict/mod "$f" build_strict/obj/*.o -o "$exe"
  "$exe"
done
