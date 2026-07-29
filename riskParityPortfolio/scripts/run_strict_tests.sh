#!/usr/bin/env sh
set -eu

fc=${FC:-gfortran}
flags="-std=f2018 -O0 -g -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
mkdir -p build
rm -f ./*.o build/*.mod

$fc $flags -J build -I build -c \
  src/rpp_kinds.f90 \
  src/rpp_types.f90 \
  src/rpp_linalg.f90 \
  src/rpp_core.f90 \
  src/rpp_formulations.f90 \
  src/rpp_qp.f90 \
  src/rpp_solvers.f90 \
  src/rpp_api.f90 \
  src/risk_parity_portfolio.f90

for source in test/*.f90; do
  case "$source" in
    *tmp_*|*debug_*) continue ;;
  esac
  exe="build/$(basename "$source" .f90)"
  $fc $flags -I build "$source" ./*.o -llapack -lblas -o "$exe"
  "$exe"
done
