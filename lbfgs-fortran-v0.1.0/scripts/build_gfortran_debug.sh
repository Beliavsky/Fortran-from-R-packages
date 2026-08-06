#!/usr/bin/env sh
set -eu

build_dir=${1:-build-gfortran-debug}
rm -rf "$build_dir"
mkdir -p "$build_dir"

flags="-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -fcheck=all"
flags="$flags -ffpe-trap=invalid,zero,overflow -fbacktrace"
mods="-J $build_dir -I $build_dir"
sources="src/lbfgs_kinds.f90 src/lbfgs_status.f90 src/lbfgs_solver.f90 src/lbfgs.f90"

gfortran $flags $mods $sources test/test_lbfgs.f90 -o "$build_dir/test_lbfgs"
"$build_dir/test_lbfgs"
