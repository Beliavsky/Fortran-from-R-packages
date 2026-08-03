#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-debug"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"
sources="
$root/src/kernlab_kinds.f90
$root/src/kernlab_types.f90
$root/src/kernlab_linalg.f90
$root/src/kernlab_kernels.f90
$root/src/kernlab_core.f90
$root/src/kernlab_unsupervised.f90
$root/src/kernlab_mmd.f90
$root/src/kernlab_supervised.f90
$root/src/kernlab.f90"
# shellcheck disable=SC2086
gfortran $flags -J "$build/mod" -I "$build/mod" -c $sources
mv ./*.o "$build/"
for source in "$root"/test/*.f90 "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  # shellcheck disable=SC2086
  gfortran $flags -J "$build/mod" -I "$build/mod" "$source" "$build"/*.o -o "$build/bin/$name"
  "$build/bin/$name"
done
