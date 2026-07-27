#!/usr/bin/env sh
set -eu

project_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

build_and_run() {
  build_dir=$1
  flags=$2
  rm -rf "$build_dir"
  mkdir -p "$build_dir"
  cd "$build_dir"

  # shellcheck disable=SC2086
  gfortran -c $flags \
    "$project_root/src/svdnf_kinds.f90" \
    "$project_root/src/svdnf_types.f90" \
    "$project_root/src/svdnf_stats.f90" \
    "$project_root/src/svdnf_models.f90" \
    "$project_root/src/svdnf_grids.f90" \
    "$project_root/src/svdnf_filter.f90" \
    "$project_root/src/svdnf_simulation.f90" \
    "$project_root/src/svdnf_optimization.f90" \
    "$project_root/src/svdnf.f90" \
    "$project_root/test/test_support.f90" \
    "$project_root/test/test_custom_callbacks.f90"

  objects="svdnf_kinds.o svdnf_types.o svdnf_stats.o svdnf_models.o svdnf_grids.o"
  objects="$objects svdnf_filter.o svdnf_simulation.o svdnf_optimization.o svdnf.o"
  test_objects="$objects test_support.o test_custom_callbacks.o"

  for name in probabilities models filter jumps forecast optimization; do
    # shellcheck disable=SC2086
    gfortran $flags "$project_root/test/test_${name}.f90" $test_objects -o "test_${name}"
    "./test_${name}"
  done

  for source in "$project_root/app/svdnf_demo.f90" \
    "$project_root/example/built_in_models.f90" \
    "$project_root/example/filter_and_forecast.f90" \
    "$project_root/example/estimate_taylor.f90"; do
    executable=$(basename "$source" .f90)
    # shellcheck disable=SC2086
    gfortran $flags "$source" $objects -o "$executable"
    "./$executable" >/dev/null
  done
}

checked_flags="-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
optimized_flags="-O2 -std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror"

build_and_run "$project_root/build/checked" "$checked_flags"
build_and_run "$project_root/build/optimized" "$optimized_flags"

printf '%s\n' 'validation: PASS'
