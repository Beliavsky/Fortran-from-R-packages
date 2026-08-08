#!/usr/bin/env sh
set -eu
mkdir -p build/reference_c
gcc -std=c99 -O0 \
  -I validation/c_reference \
  -I original/stochQN-master/inst/include \
  original/stochQN-master/src/stochqn.c \
  validation/c_reference/reference_olbfgs.c \
  -lm -o build/reference_c/reference_olbfgs
./build/reference_c/reference_olbfgs
