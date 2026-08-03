#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build_optimized
mkdir -p build_optimized/mod build_optimized/obj build_optimized/bin

flags="-std=f2018 -O3 -Wall -Wextra -Werror -pedantic -ffree-line-length-none"
modules="r4gpf_kinds r4gpf_status r4gpf_linalg r4gpf_random r4gpf_optimization r4gpf_mortality r4gpf_finance r4gpf_portfolio r4gpf_household r4gpf_simulation r4good_personal_finances"
objects=""
for module in $modules; do
  gfortran $flags -J build_optimized/mod -I build_optimized/mod -c "src/$module.f90" -o "build_optimized/obj/$module.o"
  objects="$objects build_optimized/obj/$module.o"
done

for source in test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I build_optimized/mod $objects "$source" -o "build_optimized/bin/$name"
  "build_optimized/bin/$name"
done
for source in example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I build_optimized/mod $objects "$source" -o "build_optimized/bin/$name"
  "build_optimized/bin/$name"
done
