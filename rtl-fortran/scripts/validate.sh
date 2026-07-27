#!/usr/bin/env sh
# SPDX-License-Identifier: MIT
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

compile_and_run() {
  build_dir=$1
  optimization=$2
  rm -rf "$build_dir"
  mkdir -p "$build_dir/mod" "$build_dir/obj"

  flags="-std=f2018 $optimization -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -fbacktrace"
  if [ "$optimization" = "-O0" ]; then
    flags="$flags -fcheck=all"
  fi

  for unit in rtl_kinds rtl_types rtl_stats rtl_calendar rtl_options rtl_processes \
      rtl_fixed_income rtl_portfolio rtl_market rtl; do
    gfortran $flags -J"$build_dir/mod" -I"$build_dir/mod" \
      -c "src/$unit.f90" -o "$build_dir/obj/$unit.o"
  done

  objects=""
  for unit in rtl_kinds rtl_types rtl_stats rtl_calendar rtl_options rtl_processes \
      rtl_fixed_income rtl_portfolio rtl_market rtl; do
    objects="$objects $build_dir/obj/$unit.o"
  done

  for source in test/test_*.f90; do
    name=$(basename "$source" .f90)
    gfortran $flags -I"$build_dir/mod" "$source" $objects -o "$build_dir/$name"
    "$build_dir/$name"
  done

  for source in app/*.f90 example/*.f90; do
    name=$(basename "$source" .f90)
    gfortran $flags -I"$build_dir/mod" "$source" $objects -o "$build_dir/$name"
    "$build_dir/$name" >/dev/null
  done
}

compile_and_run build_validate_debug -O0
compile_and_run build_validate_opt -O2
printf '%s\n' 'validation: PASS'
