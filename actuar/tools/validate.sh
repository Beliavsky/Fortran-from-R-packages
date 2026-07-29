#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
FC=${FC:-gfortran}

SOURCES=(
  "$ROOT/src/actuar_kinds.f90"
  "$ROOT/src/actuar_special.f90"
  "$ROOT/src/actuar_rng.f90"
  "$ROOT/src/actuar_types.f90"
  "$ROOT/src/actuar_continuous.f90"
  "$ROOT/src/actuar_supplements.f90"
  "$ROOT/src/actuar_discrete.f90"
  "$ROOT/src/actuar_aggregate.f90"
  "$ROOT/src/actuar_phase_type.f90"
  "$ROOT/src/actuar_credibility.f90"
  "$ROOT/src/actuar_grouped.f90"
  "$ROOT/src/actuar_risk.f90"
  "$ROOT/src/actuar.f90"
)

OBJECTS=(
  actuar_kinds.o
  actuar_special.o
  actuar_rng.o
  actuar_types.o
  actuar_continuous.o
  actuar_supplements.o
  actuar_discrete.o
  actuar_aggregate.o
  actuar_phase_type.o
  actuar_credibility.o
  actuar_grouped.o
  actuar_risk.o
  actuar.o
)

run_build() {
  local mode=$1
  local flags=$2
  local build_dir="$ROOT/build/validation_$mode"

  rm -rf "$build_dir"
  mkdir -p "$build_dir"
  cd "$build_dir"

  "$FC" $flags -c "${SOURCES[@]}"

  for source in "$ROOT"/test/*.f90; do
    local name
    name=$(basename "$source" .f90)
    "$FC" $flags -I. "$source" "${OBJECTS[@]}" -o "$name"
    "./$name"
  done

  "$FC" $flags -I. "$ROOT/app/actuar_demo.f90" "${OBJECTS[@]}" -o actuar_demo
  ./actuar_demo >/dev/null

  for source in "$ROOT"/example/*.f90; do
    local name
    name=$(basename "$source" .f90)
    "$FC" $flags -I. "$source" "${OBJECTS[@]}" -o "$name"
    "./$name" >/dev/null
  done

  echo "validation ($mode): PASS"
}

DEBUG_FLAGS="-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -fcheck=all -fbacktrace -Werror"
OPT_FLAGS="-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -Werror"

run_build debug "$DEBUG_FLAGS"
run_build optimized "$OPT_FLAGS"
