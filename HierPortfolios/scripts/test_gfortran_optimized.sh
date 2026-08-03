#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran-optimized"
rm -rf "$build"
mkdir -p "$build"

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3"
sources="
$root/src/hierportfolios_kinds.f90
$root/src/hierportfolios_types.f90
$root/src/hierportfolios_hierarchy.f90
$root/src/hierportfolios_core.f90
$root/src/hierportfolios.f90
"

cd "$build"
gfortran $flags -J . -I . -c $sources
for source in "$root"/test/*.f90 "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -J . -I . ./*.o "$source" -o "$name"
  "./$name"
done
