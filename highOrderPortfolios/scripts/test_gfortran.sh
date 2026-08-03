#!/bin/sh
set -eu

cd "$(dirname "$0")/.."

sources="
src/fitheavytail_kinds.f90
src/fitheavytail_status.f90
src/fitheavytail_types.f90
src/fitheavytail_linalg.f90
src/fitheavytail_rng.f90
src/fitheavytail_special.f90
src/fitheavytail_tail.f90
src/fitheavytail_mvst.f90
src/highorder_types.f90
src/highorder_linalg.f90
src/highorder_moments.f90
src/highorder_optimization.f90
src/highorderportfolios.f90
"

build_one() {
  name=$1
  flags=$2
  dir="build/$name"
  rm -rf "$dir"
  mkdir -p "$dir/mod"
  objects=""
  for source in $sources; do
    object="$dir/$(basename "$source" .f90).o"
    gfortran $flags -J "$dir/mod" -I "$dir/mod" -c "$source" -o "$object"
    objects="$objects $object"
  done
  for source in test/*.f90; do
    exe="$dir/$(basename "$source" .f90)"
    gfortran $flags -J "$dir/mod" -I "$dir/mod" "$source" $objects -o "$exe"
    "$exe"
  done
  for source in example/*.f90 app/*.f90; do
    exe="$dir/$(basename "$source" .f90)"
    gfortran $flags -J "$dir/mod" -I "$dir/mod" "$source" $objects -o "$exe"
    "$exe"
  done
}

common="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fbacktrace"
build_one debug "$common -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"
build_one release "$common -O3"
