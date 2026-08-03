#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build_checked
mkdir -p build_checked/mod build_checked/obj build_checked/bin

flags="-std=f2018 -Wall -Wextra -Werror -pedantic -ffree-line-length-none -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
modules="r4gpf_kinds r4gpf_status r4gpf_linalg r4gpf_random r4gpf_optimization r4gpf_mortality r4gpf_finance r4gpf_portfolio r4gpf_household r4gpf_simulation r4good_personal_finances"
objects=""
for module in $modules; do
  gfortran $flags -J build_checked/mod -I build_checked/mod -c "src/$module.f90" -o "build_checked/obj/$module.o"
  objects="$objects build_checked/obj/$module.o"
done

for source in test/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I build_checked/mod $objects "$source" -o "build_checked/bin/$name"
  "build_checked/bin/$name"
done
for source in example/*.f90 app/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I build_checked/mod $objects "$source" -o "build_checked/bin/$name"
  "build_checked/bin/$name"
done
