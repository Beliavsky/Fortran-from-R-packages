#!/usr/bin/env sh
set -eu
mkdir -p build/release
gfortran -std=f2018 -O3 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface \
  -J build/release -I build/release src/stochqn_kinds.f90 src/stochqn_core.f90 \
  src/stochqn_guided.f90 src/stochqn_logistic.f90 test/test_stochqn.f90 \
  -o build/release/test_stochqn
./build/release/test_stochqn
