#!/usr/bin/env sh
set -eu
fc=${FC:-gfortran}
flags="-std=f2018 -O2 -Wall -Wextra -Wpedantic -Wimplicit-interface -Werror=implicit-interface"
mkdir -p build/manual
$fc $flags -J build/manual -I build/manual -c src/nonneg_cg.f90 -o build/manual/nonneg_cg.o
$fc $flags -J build/manual -I build/manual build/manual/nonneg_cg.o test/test_nonneg_cg.f90 -o build/manual/test_nonneg_cg
$fc $flags -J build/manual -I build/manual build/manual/nonneg_cg.o example/rosenbrock.f90 -o build/manual/rosenbrock
build/manual/test_nonneg_cg
build/manual/rosenbrock
