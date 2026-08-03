#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build_gfortran
mkdir -p build_gfortran/mod build_gfortran/bin build_gfortran/obj
FLAGS="${FFLAGS:--std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow}"
for src in src/indgen_kinds.f90 src/indgen_types.f90 src/indgen_special.f90 \
  src/indgen_cvm_tables.f90 src/indgen_moebius.f90 src/indgen_core.f90 \
  src/indgenerrors.f90; do
  obj="build_gfortran/obj/$(basename "${src%.f90}").o"
  gfortran $FLAGS -J build_gfortran/mod -I build_gfortran/mod -c "$src" -o "$obj"
done
for src in test/*.f90; do
  exe="build_gfortran/bin/$(basename "${src%.f90}")"
  gfortran $FLAGS -I build_gfortran/mod build_gfortran/obj/*.o "$src" -o "$exe"
  "$exe"
done
for src in example/*.f90 app/*.f90; do
  exe="build_gfortran/bin/$(basename "${src%.f90}")"
  gfortran $FLAGS -I build_gfortran/mod build_gfortran/obj/*.o "$src" -o "$exe"
  "$exe"
done
