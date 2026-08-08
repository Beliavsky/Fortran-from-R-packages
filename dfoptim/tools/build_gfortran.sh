#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-gfortran"
rm -rf "$build"
mkdir -p "$build"

sources="
$root/src/dfoptim_kinds.f90
$root/src/dfoptim_interfaces.f90
$root/src/dfoptim_rng.f90
$root/src/dfoptim_utils.f90
$root/src/dfoptim_hooke_jeeves.f90
$root/src/dfoptim_nelder_mead.f90
$root/src/dfoptim_mads.f90
$root/src/dfoptim.f90
"

flags="-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
# shellcheck disable=SC2086
gfortran $flags -J "$build" -I "$build" $sources "$root/test/test_dfoptim.f90" -o "$build/test_dfoptim"
"$build/test_dfoptim"
# shellcheck disable=SC2086
gfortran $flags -J "$build" -I "$build" $sources "$root/example/dfoptim_example.f90" -o "$build/dfoptim_example"
"$build/dfoptim_example"
