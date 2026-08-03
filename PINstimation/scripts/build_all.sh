#!/usr/bin/env sh
set -eu
mode=${1:-check}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-$mode"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
common="-std=f2018 -Wall -Wextra -Wpedantic -Wimplicit-interface -ffree-line-length-none"
if [ "$mode" = "release" ]; then
  flags="$common -O3 -DNDEBUG"
else
  flags="$common -O0 -g -fcheck=all -fbacktrace -finit-real=snan -finit-integer=-999999"
fi
order="pinstimation_kinds pinstimation_types pinstimation_math pinstimation_optimization pinstimation_pin pinstimation_mpin pinstimation_adjpin pinstimation_data pinstimation_vpin pinstimation"
for unit in $order; do
  gfortran $flags -J"$build/mod" -I"$build/mod" -c "$root/src/$unit.f90" -o "$build/obj/$unit.o"
done
objects=$(find "$build/obj" -name '*.o' -print | sort)
for source in "$root"/test/*.f90 "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I"$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name"
done
printf '%s\n' "All $mode builds and runs passed."
