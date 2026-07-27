#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build-validation"
rm -rf "$build"
mkdir -p "$build"
cd "$build"
flags="-std=f2018 -O0 -g -Wall -Wextra -Werror -Wconversion-extra -Wimplicit-interface -fcheck=all -fbacktrace -J . -I ."
for f in garchsk_kinds garchsk_types garchsk_stats garchsk_linalg garchsk_models garchsk_estimation garchsk; do
  gfortran $flags -c "$root/src/$f.f90"
done
gfortran $flags -c "$root/test/test_support.f90"
lib="garchsk_kinds.o garchsk_types.o garchsk_stats.o garchsk_linalg.o garchsk_models.o garchsk_estimation.o garchsk.o"
for t in test_constraints test_estimation test_likelihood_forecast test_moments; do
  gfortran $flags "$root/test/$t.f90" $lib test_support.o -o "$t"
  "./$t"
done
for source in "$root/app/garchsk_demo.f90" "$root/example/basic_garchsk.f90" "$root/example/estimate_and_forecast.f90"; do
  exe=$(basename "$source" .f90)
  gfortran $flags "$source" $lib -o "$exe"
  "./$exe" >/dev/null
done
echo "validation: PASS"
