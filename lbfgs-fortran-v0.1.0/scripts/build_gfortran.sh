#!/usr/bin/env sh
set -eu

build_dir=${1:-build-gfortran}
rm -rf "$build_dir"
mkdir -p "$build_dir"

flags="-std=f2018 -O2 -Wall -Wextra -Wpedantic"
mods="-J $build_dir -I $build_dir"
sources="src/lbfgs_kinds.f90 src/lbfgs_status.f90 src/lbfgs_solver.f90 src/lbfgs.f90"

gfortran $flags $mods $sources test/test_lbfgs.f90 -o "$build_dir/test_lbfgs"
"$build_dir/test_lbfgs"

gfortran $flags $mods $sources example/rosenbrock.f90 -o "$build_dir/rosenbrock"
"$build_dir/rosenbrock"

gfortran $flags $mods $sources example/owlqn_soft_threshold.f90 \
    -o "$build_dir/owlqn_soft_threshold"
"$build_dir/owlqn_soft_threshold"
