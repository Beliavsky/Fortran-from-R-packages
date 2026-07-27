#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build"

sources="
$root/src/bidask_kinds.f90
$root/src/bidask_types.f90
$root/src/bidask_statistics.f90
$root/src/bidask_rng.f90
$root/src/bidask_estimators.f90
$root/src/bidask_windows.f90
$root/src/bidask_simulation.f90
$root/src/bidask.f90
"

run_build() {
  name=$1
  flags=$2
  dir="$build/$name"
  mkdir -p "$dir"
  cd "$dir"
  # shellcheck disable=SC2086
  gfortran $flags -J . -I . -c $sources
  objects="bidask_kinds.o bidask_types.o bidask_statistics.o bidask_rng.o bidask_estimators.o bidask_windows.o bidask_simulation.o bidask.o"
  for source in "$root"/test/*.f90; do
    exe=$(basename "$source" .f90)
    # shellcheck disable=SC2086
    gfortran $flags -J . -I . "$source" $objects -o "$exe"
    ./"$exe"
  done
  for source in "$root"/app/*.f90 "$root"/example/*.f90; do
    exe=$(basename "$source" .f90)
    # shellcheck disable=SC2086
    gfortran $flags -J . -I . "$source" $objects -o "$exe"
    ./"$exe" >/dev/null
  done
}

run_build checked "-std=f2018 -O0 -g -Wall -Wextra -Werror -fcheck=all -fbacktrace"
run_build optimized "-std=f2018 -O2 -Wall -Wextra -Werror"
printf '%s\n' 'validation: PASS'
