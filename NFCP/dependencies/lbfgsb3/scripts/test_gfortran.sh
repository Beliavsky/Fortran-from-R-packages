#!/usr/bin/env sh
set -eu
rm -rf build-gfortran
mkdir build-gfortran
gfortran -c -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow src/lbfgsb3_core.f90 \
  -J build-gfortran -o build-gfortran/lbfgsb3_core.o
gfortran -c -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow src/lbfgsb3.f90 \
  -I build-gfortran -J build-gfortran -o build-gfortran/lbfgsb3.o
gfortran -std=f2018 -O0 -g -fcheck=all \
  -ffpe-trap=invalid,zero,overflow test/test_lbfgsb3.f90 \
  build-gfortran/lbfgsb3.o build-gfortran/lbfgsb3_core.o \
  -I build-gfortran -J build-gfortran -o build-gfortran/test_lbfgsb3
./build-gfortran/test_lbfgsb3
