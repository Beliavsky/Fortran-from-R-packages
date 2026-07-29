#!/bin/sh
# SPDX-License-Identifier: GPL-2.0-or-later
set -eu
mode=${1:-strict}
case "$mode" in
  strict)
    flags='-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow -O0 -g'
    ;;
  optimized)
    flags='-std=f2018 -Wall -Wextra -Werror -pedantic -O3'
    ;;
  *)
    echo "usage: $0 [strict|optimized]" >&2
    exit 2
    ;;
esac
build="build-$mode"
rm -rf "$build"
mkdir "$build"
cd "$build"
gfortran $flags -c ../src/vares_kinds.f90
gfortran $flags -c ../src/vares_special.f90
gfortran $flags -c ../src/vares_quadrature.f90
for source in ../src/vares_distributions_*.f90; do
  gfortran $flags -c "$source"
done
gfortran $flags -c ../src/vares.f90
for test_source in ../test/*.f90; do
  name=$(basename "$test_source" .f90)
  gfortran $flags "$test_source" ./*.o -o "$name"
  "./$name"
done
gfortran $flags ../app/main.f90 ./*.o -o vares_demo
./vares_demo
gfortran $flags ../example/vector_example.f90 ./*.o -o vector_example
./vector_example
