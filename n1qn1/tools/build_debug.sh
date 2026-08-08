#!/usr/bin/env sh
set -eu
rm -rf build-debug
mkdir build-debug
gfortran -std=f2018 -Wall -Wextra -Wpedantic -Wno-compare-reals \
  -Wimplicit-interface -Werror=implicit-interface -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -fbacktrace \
  -J build-debug -I build-debug src/n1qn1.f90 test/test_n1qn1.f90 \
  -o build-debug/test_n1qn1
./build-debug/test_n1qn1
gfortran -std=f2018 -Wall -Wextra -Wpedantic -Wno-compare-reals \
  -Wimplicit-interface -Werror=implicit-interface -fcheck=all \
  -ffpe-trap=invalid,zero,overflow -fbacktrace \
  -J build-debug -I build-debug src/n1qn1.f90 example/banana.f90 \
  -o build-debug/banana
./build-debug/banana
