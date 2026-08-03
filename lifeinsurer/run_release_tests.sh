#!/usr/bin/env bash
set -euo pipefail
FC=${FC:-gfortran}
ROOT=$(cd "$(dirname "$0")" && pwd)
BUILD="$ROOT/build-release"
rm -rf "$BUILD"
mkdir -p "$BUILD/mod" "$BUILD/obj" "$BUILD/bin"
FLAGS=(-std=f2018 -Wall -Wextra -Wpedantic -Werror -O3 -J"$BUILD/mod" -I"$BUILD/mod")
SOURCES=(lifeinsurer_kinds lifeinsurer_types lifeinsurer_helpers lifeinsurer_pv lifeinsurer_cashflows lifeinsurer_actuarial lifeinsurer_profit lifeinsurer_contract lifeinsurer)
for source in "${SOURCES[@]}"; do "$FC" "${FLAGS[@]}" -c "$ROOT/src/$source.f90" -o "$BUILD/obj/$source.o"; done
OBJECTS=(); for source in "${SOURCES[@]}"; do OBJECTS+=("$BUILD/obj/$source.o"); done
for file in "$ROOT"/test/*.f90; do name=$(basename "${file%.f90}"); "$FC" "${FLAGS[@]}" "$file" "${OBJECTS[@]}" -o "$BUILD/bin/$name"; "$BUILD/bin/$name"; done
