#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-optimized"
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"

flags="-std=f2018 -Wall -Wextra -Wimplicit-interface -O3"
sources="magic_kinds magic_status magic_tensor magic_square magic_hypercube magic_combinatorics magic"
objects=""
for name in $sources; do
  gfortran $flags -J"$build/mod" -I"$build/mod" -c "$root/src/$name.f90" -o "$build/$name.o"
  objects="$objects $build/$name.o"
done

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J"$build/mod" -I"$build/mod" $objects "$source" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J"$build/mod" -I"$build/mod" $objects "$source" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done
