#!/usr/bin/env sh
set -eu

FC=${FC:-gfortran}
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$ROOT"

build_one() {
  mode=$1
  flags=$2
  bdir="build-$mode"
  rm -rf "$bdir"
  mkdir -p "$bdir/mod" "$bdir/obj" "$bdir/bin"

  for src in \
    src/apt_kinds.f90 \
    src/apt_special.f90 \
    src/apt_regression.f90 \
    src/apt_cointegration.f90 \
    src/apt_ecm.f90 \
    src/apt.f90; do
    obj="$bdir/obj/$(basename "${src%.f90}").o"
    $FC $flags -J"$bdir/mod" -I"$bdir/mod" -c "$src" -o "$obj"
  done

  objs=$(find "$bdir/obj" -name '*.o' -print | sort | tr '\n' ' ')
  $FC $flags -J"$bdir/mod" -I"$bdir/mod" -c test/test_support.f90 -o "$bdir/obj/test_support.o"

  for name in test_statistics test_reference test_search test_ecm_tests; do
    $FC $flags -J"$bdir/mod" -I"$bdir/mod" $objs "$bdir/obj/test_support.o" \
      "test/$name.f90" -llapack -lblas -o "$bdir/bin/$name"
    "$bdir/bin/$name"
  done

  $FC $flags -J"$bdir/mod" -I"$bdir/mod" $objs app/apt_demo.f90 \
    -llapack -lblas -o "$bdir/bin/apt_demo"
  "$bdir/bin/apt_demo" >/dev/null

  for name in threshold_search asymmetric_ecm; do
    $FC $flags -J"$bdir/mod" -I"$bdir/mod" $objs "example/$name.f90" \
      -llapack -lblas -o "$bdir/bin/$name"
    "$bdir/bin/$name" >/dev/null
  done
}

COMMON="-std=f2018 -Wall -Wextra -Werror"
build_one strict "$COMMON -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow"
build_one optimized "$COMMON -O3"

echo "All apt-fortran builds and tests passed."
