#!/usr/bin/env sh
set -eu
FFLAGS="-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow" \
  sh scripts/build_gfortran.sh
./build/bin/test_ycevo
./build/bin/ycevo_example
