#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran-optimized"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-compare-reals"

compile_module() {
  src=$1
  obj=$2
  gfortran $flags -J"$build/mod" -I"$build/mod" -c "$src" -o "$obj"
}

compile_module "$root/src/quadprog_kinds.f90" "$build/obj/quadprog_kinds.o"
compile_module "$root/src/quadprog_core.f90" "$build/obj/quadprog_core.o"
compile_module "$root/src/quadprog.f90" "$build/obj/quadprog.o"
compile_module "$root/src/test_support.f90" "$build/obj/test_support.o"

objects="$build/obj/quadprog_kinds.o $build/obj/quadprog_core.o \
$build/obj/quadprog.o $build/obj/test_support.o"

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I"$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I"$build/mod" "$source" $objects -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done

echo "GNU Fortran optimized validation: PASS"
