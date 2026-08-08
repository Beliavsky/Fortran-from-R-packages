#!/usr/bin/env sh
set -eu
mkdir -p build/debug
gfortran -std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface \
  -fcheck=all -ffpe-trap=invalid,zero,overflow -J build/debug -I build/debug \
  src/stochqn_kinds.f90 src/stochqn_core.f90 src/stochqn_guided.f90 \
  src/stochqn_logistic.f90 test/test_stochqn.f90 -o build/debug/test_stochqn
./build/debug/test_stochqn
