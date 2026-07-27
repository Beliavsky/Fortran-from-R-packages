#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build-validation
mkdir build-validation
cd build-validation
flags="-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Wno-compare-reals -Werror -fcheck=all -fbacktrace -O0"
sources="../src/ob_kinds.f90 ../src/ob_types.f90 ../src/ob_utils.f90 ../src/ob_alignment.f90 ../src/ob_io.f90 ../src/ob_events.f90 ../src/ob_trades.f90 ../src/ob_depth.f90 ../src/ob_book.f90 ../src/ob_processing.f90 ../src/ob_analytics.f90"
gfortran $flags -c -J . -I . $sources
objects="ob_kinds.o ob_types.o ob_utils.o ob_alignment.o ob_io.o ob_events.o ob_trades.o ob_depth.o ob_book.o ob_processing.o ob_analytics.o"
for source in ../test/test_*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I . "$source" $objects -o "$name"
done
for source in ../app/*.f90 ../example/*.f90; do
  name=$(basename "$source" .f90)
  gfortran $flags -I . "$source" $objects -o "$name"
done
cd ..
for executable in build-validation/test_*; do "$executable"; done
build-validation/ob_analytics_demo >/dev/null
build-validation/event_matching >/dev/null
build-validation/depth_and_book >/dev/null
echo "validation: PASS"
