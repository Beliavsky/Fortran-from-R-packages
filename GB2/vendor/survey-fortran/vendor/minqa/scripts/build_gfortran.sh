#!/usr/bin/env sh
set -eu

rm -rf build-standalone
mkdir -p build-standalone/mod

gfortran -std=f2018 -Wall -Wextra -Wpedantic -fcheck=all -fbacktrace \
  -J build-standalone/mod \
  src/minqa_module.f90 test/test_minqa.f90 \
  -o build-standalone/test_minqa

./build-standalone/test_minqa

gfortran -std=f2018 -O3 \
  -J build-standalone/mod \
  src/minqa_module.f90 example/minqa_example.f90 \
  -o build-standalone/minqa_example

./build-standalone/minqa_example
