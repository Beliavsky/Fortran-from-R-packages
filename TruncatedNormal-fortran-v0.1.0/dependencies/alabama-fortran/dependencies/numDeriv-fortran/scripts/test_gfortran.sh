#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-strict"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags='-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace'
for source in numderiv_kinds numderiv_types numderiv_callbacks numderiv_core numderiv; do
  gfortran $flags -J "$build/mod" -I "$build/mod" -c "$root/src/$source.f90" -o "$build/obj/$source.o"
done
objects="$build/obj/numderiv_kinds.o $build/obj/numderiv_types.o $build/obj/numderiv_callbacks.o $build/obj/numderiv_core.o $build/obj/numderiv.o"
for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J "$build/mod" -I "$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name"
done
