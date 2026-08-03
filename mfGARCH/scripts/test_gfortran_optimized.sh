#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
build="$root/build/gfortran-release"
rm -rf "$build"
mkdir -p "$build/obj" "$build/bin"
cd "$build/obj"

flags='-std=f2018 -O3 -march=native -Wall -Wextra -Werror -Wimplicit-interface'

# shellcheck disable=SC2086
gfortran $flags -c \
  "$root/src/mfgarch_kinds.f90" \
  "$root/src/mfgarch_status.f90" \
  "$root/src/mfgarch_types.f90" \
  "$root/src/mfgarch_math.f90" \
  "$root/src/mfgarch_components.f90" \
  "$root/src/mfgarch_low_level.f90" \
  "$root/src/mfgarch_optimization.f90" \
  "$root/src/mfgarch_random.f90" \
  "$root/src/mfgarch_simulation.f90" \
  "$root/src/mfgarch_prediction.f90" \
  "$root/src/mfgarch_fit.f90" \
  "$root/src/mfgarch_io.f90" \
  "$root/src/mfgarch.f90"

for source in "$root"/test/*.f90; do
  name=$(basename "$source" .f90)
  # shellcheck disable=SC2086
  gfortran $flags ./*.o "$source" -o "$build/bin/$name"
  "$build/bin/$name"
done

for source in "$root"/example/*.f90 "$root"/app/*.f90; do
  name=$(basename "$source" .f90)
  # shellcheck disable=SC2086
  gfortran $flags ./*.o "$source" -o "$build/bin/$name"
  (cd "$build" && "$build/bin/$name")
done
