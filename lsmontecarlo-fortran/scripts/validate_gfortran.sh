#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."
rm -rf build
mkdir -p build/mod build/obj build/bin

flags=(
  -std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra
  -Wimplicit-interface -Werror -fcheck=all -fbacktrace
  -ffree-line-length-132
)

sources=(
  src/lsmc_kinds.f90
  src/lsmc_math.f90
  src/lsmc_random.f90
  src/lsmc_linear_algebra.f90
  src/lsmc_types.f90
  src/lsmc_european.f90
  src/lsmc_simulation.f90
  src/lsmc_pricing.f90
  src/lsmc_surface.f90
  src/lsmontecarlo.f90
)

objects=()
for source in "${sources[@]}"; do
  object="build/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -Jbuild/mod -Ibuild/mod -c "$source" -o "$object"
  objects+=("$object")
done

for source in test/*.f90 app/*.f90 example/*.f90; do
  executable="build/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" -Jbuild/mod -Ibuild/mod "$source" "${objects[@]}" -o "$executable"
  "$executable"
done
