#!/usr/bin/env sh
# SPDX-License-Identifier: Artistic-2.0
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-optimized"
rm -rf "$build"
mkdir -p "$build"

flags="-std=f2018 -O3 -march=native -Wall -Wextra -Wimplicit-interface"
sources="
ldhmm_kinds
ldhmm_status
ldhmm_math
ldhmm_types
ldhmm_distribution
ldhmm_parameters
ldhmm_modeling
ldhmm_simulation
ldhmm_optimization
ldhmm_series
ldhmm
"
objects=""
for name in $sources; do
    gfortran $flags -J "$build" -I "$build" -c "$root/src/$name.f90" \
        -o "$build/$name.o"
    objects="$objects $build/$name.o"
done

for source in "$root"/test/*.f90; do
    name=$(basename "$source" .f90)
    echo "Running $name"
    gfortran $flags -I "$build" "$source" $objects -o "$build/$name"
    "$build/$name"
done

echo "All optimized GNU Fortran tests passed."
