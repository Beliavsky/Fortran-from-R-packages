#!/usr/bin/env bash
set -euo pipefail
FC=${FC:-gfortran}
FLAGS="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0"
rm -rf build-strict
mkdir -p build-strict/mod build-strict/bin
$FC $FLAGS -c src/cec2005benchmark.f90 -J build-strict/mod -o build-strict/cec2005benchmark.o
for src in test/*.f90; do
  name=$(basename "$src" .f90)
  $FC $FLAGS -I build-strict/mod "$src" build-strict/cec2005benchmark.o -o "build-strict/bin/$name"
  "build-strict/bin/$name"
done
for src in example/*.f90; do
  name=$(basename "$src" .f90)
  $FC $FLAGS -I build-strict/mod "$src" build-strict/cec2005benchmark.o -o "build-strict/bin/example_$name"
  "build-strict/bin/example_$name"
done
