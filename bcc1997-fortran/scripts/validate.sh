#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/.validation-build"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"

flags="-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
sources="bcc1997_kinds bcc1997_types bcc1997_quadrature bcc1997_model bcc1997"
objects=""

for name in $sources; do
  gfortran $flags -J "$build/mod" -I "$build/mod" -c \
    "$root/src/$name.f90" -o "$build/obj/$name.o"
  objects="$objects $build/obj/$name.o"
done

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I "$build/mod" $objects "$source" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I "$build/mod" $objects "$source" -o "$build/bin/$name"
  "$build/bin/$name" >/dev/null
done

rm -rf "$build"
printf '%s\n' 'validation: PASS'
