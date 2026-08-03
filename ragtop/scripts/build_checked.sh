#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
rm -rf build/checked
mkdir -p build/checked/{mod,obj,bin}
flags=(-std=f2018 -Wall -Wextra -Werror -fcheck=all \
       -ffpe-trap=invalid,zero,overflow -fbacktrace)
sources=(
  ragtop_kinds ragtop_constants ragtop_types ragtop_math
  ragtop_term_structures ragtop_cashflows ragtop_black_scholes
  ragtop_instruments ragtop_pde ragtop_calibration ragtop_greeks
  ragtop_validation ragtop
)
for name in "${sources[@]}"; do
  gfortran "${flags[@]}" -J build/checked/mod -I build/checked/mod \
    -c "src/${name}.f90" -o "build/checked/obj/${name}.o"
done
objects=(build/checked/obj/*.o)
for source in test/*.f90 example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" -J build/checked/mod -I build/checked/mod \
    "$source" "${objects[@]}" -o "build/checked/bin/$name"
  "build/checked/bin/$name"
done
