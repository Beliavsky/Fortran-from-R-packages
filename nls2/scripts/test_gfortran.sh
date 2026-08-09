#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build_strict
mkdir -p build_strict/mod build_strict/bin
flags=(-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all)
src=(
  src/nls2_kinds.f90
  src/nls2_types.f90
  src/nls2_linalg.f90
  src/nls2_random.f90
  src/nls2_core.f90
  src/nls2_search.f90
  src/nls2_stats.f90
  src/nls2.f90
)
gfortran "${flags[@]}" -c -J build_strict/mod -I build_strict/mod "${src[@]}"
mv ./*.o build_strict/
for file in test/test_*.f90; do
  name=$(basename "$file" .f90)
  gfortran "${flags[@]}" -J build_strict/mod -I build_strict/mod build_strict/*.o "$file" -o "build_strict/bin/$name"
  "build_strict/bin/$name"
done
for file in example/*.f90; do
  name=$(basename "$file" .f90)
  gfortran "${flags[@]}" -J build_strict/mod -I build_strict/mod build_strict/*.o "$file" -o "build_strict/bin/$name"
  "build_strict/bin/$name"
done
