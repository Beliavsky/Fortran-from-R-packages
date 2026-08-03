#!/usr/bin/env sh
set -eu

mode=${1:-checked}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/$mode"

case "$mode" in
  checked)
    flags='-O0 -g -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace'
    ;;
  optimized)
    flags='-O3 -std=f2018 -Wall -Wextra -Werror'
    ;;
  *)
    echo "usage: $0 [checked|optimized]" >&2
    exit 2
    ;;
esac

rm -rf "$build"
mkdir -p "$build"
cd "$build"

sources='multiatsm_kinds multiatsm_linalg multiatsm_types multiatsm_random multiatsm_pca multiatsm_var multiatsm_jll multiatsm_affine multiatsm_likelihood multiatsm_outputs multiatsm_optimization multiatsm_bootstrap multiatsm_bias multiatsm'
for source in $sources; do
  gfortran $flags -J. -I. -c "$root/src/$source.f90"
done

for source in "$root"/test/*.f90; do
  exe=$(basename "$source" .f90)
  gfortran $flags -I. "$source" ./*.o -llapack -lblas -o "$exe"
  "./$exe"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  exe=$(basename "$source" .f90)
  gfortran $flags -I. "$source" ./*.o -llapack -lblas -o "$exe"
  "./$exe"
done
