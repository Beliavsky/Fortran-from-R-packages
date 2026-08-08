#!/usr/bin/env sh
set -eu
rm -rf build-gfortran
mkdir -p build-gfortran/mod
flags="-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
gfortran $flags -J build-gfortran/mod src/mize.f90 test/test_mize.f90 -o build-gfortran/test_mize
./build-gfortran/test_mize
