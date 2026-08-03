#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build/optimized
mkdir -p build/optimized/{mod,obj,bin}
flags=(-std=f2018 -O3 -march=native -Wall -Wextra -Werror)
sources=(
  ragtop_kinds ragtop_constants ragtop_types ragtop_math
  ragtop_term_structures ragtop_cashflows ragtop_black_scholes
  ragtop_instruments ragtop_pde ragtop_calibration ragtop_greeks
  ragtop_validation ragtop
)
for name in "${sources[@]}"; do
  gfortran "${flags[@]}" -J build/optimized/mod -I build/optimized/mod \
    -c "src/${name}.f90" -o "build/optimized/obj/${name}.o"
done
objects=(build/optimized/obj/*.o)
for source in test/*.f90 example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" -J build/optimized/mod -I build/optimized/mod \
    "$source" "${objects[@]}" -o "build/optimized/bin/$name"
  "build/optimized/bin/$name"
done
