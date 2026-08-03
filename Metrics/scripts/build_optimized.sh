#!/usr/bin/env sh
set -eu

fc=${FC:-gfortran}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/optimized"
rm -rf "$build"
mkdir -p "$build"

flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals -Wno-intrinsic-shadow"
modules="metrics_kinds metrics_utils metrics_regression metrics_binary_classification metrics_classification metrics_information_retrieval metrics_time_series metrics metrics_test_support"

for name in $modules; do
    "$fc" $flags -J "$build" -I "$build" -c "$root/src/$name.f90" -o "$build/$name.o"
done

objects=""
for name in $modules; do objects="$objects $build/$name.o"; done

for source in "$root"/test/*.f90; do
    name=$(basename "$source" .f90)
    "$fc" $flags -J "$build" -I "$build" "$source" $objects -o "$build/$name"
    "$build/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
    name=$(basename "$source" .f90)
    "$fc" $flags -J "$build" -I "$build" "$source" $objects -o "$build/$name"
    "$build/$name"
done
