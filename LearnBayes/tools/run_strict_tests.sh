#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
build="$root/build/strict"
rm -rf "$build"
mkdir -p "$build"
cd "$build"

python "$root/tools/check_source_rules.py"

flags=(
  -std=f2018
  -Wall
  -Wextra
  -Werror
  -Wimplicit-interface
  -fimplicit-none
  -fcheck=all
  -ffpe-trap=invalid,zero,overflow
  -I.
  -J.
)

sources=(
  learnbayes_kinds
  learnbayes_special
  learnbayes_linalg
  learnbayes_rng
  learnbayes_distributions
  learnbayes_types
  learnbayes_discrete
  learnbayes_models
  learnbayes_sampling
  learnbayes_regression
  learnbayes_hierarchical
  learnbayes_utilities
  learnbayes
)

for name in "${sources[@]}"; do
  gfortran "${flags[@]}" -c "$root/src/$name.f90"
done

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" "$source" ./*.o -o "$name"
  "./$name"
done

for source in "$root"/example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran "${flags[@]}" "$source" ./*.o -o "$name"
  "./$name"
done
