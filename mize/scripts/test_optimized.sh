#!/usr/bin/env sh
set -eu
rm -rf build-optimized
mkdir -p build-optimized/mod
flags="-std=f2018 -O3 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface"
gfortran $flags -J build-optimized/mod src/mize.f90 test/test_mize.f90 -o build-optimized/test_mize
./build-optimized/test_mize
gfortran $flags -J build-optimized/mod src/mize.f90 example/rosenbrock.f90 -o build-optimized/rosenbrock
./build-optimized/rosenbrock
gfortran $flags -J build-optimized/mod src/mize.f90 example/stateful.f90 -o build-optimized/stateful
./build-optimized/stateful
