#!/usr/bin/env sh
set -eu
rm -rf build-optimized
mkdir build-optimized
gfortran -c -std=f2018 -O3 src/lbfgsb3_core.f90 \
  -J build-optimized -o build-optimized/lbfgsb3_core.o
gfortran -c -std=f2018 -O3 src/lbfgsb3.f90 \
  -I build-optimized -J build-optimized -o build-optimized/lbfgsb3.o
gfortran -std=f2018 -O3 test/test_lbfgsb3.f90 \
  build-optimized/lbfgsb3.o build-optimized/lbfgsb3_core.o \
  -I build-optimized -J build-optimized -o build-optimized/test_lbfgsb3
./build-optimized/test_lbfgsb3
