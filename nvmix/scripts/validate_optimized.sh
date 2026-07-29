#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
build=build/optimized
rm -rf "$build"
mkdir -p "$build/mod" "$build/bin"
flags=(-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror)
sources=(
  src/nvmix_kinds.f90 src/nvmix_types.f90 src/nvmix_special.f90
  src/nvmix_random.f90 src/nvmix_linalg.f90 src/nvmix_mixing.f90
  src/nvmix_core.f90 src/nvmix_gamma_mix.f90 src/nvmix_distributions.f90
  src/nvmix_risk_dependence.f90 src/nvmix_skewt.f90 src/nvmix_fitting.f90
  src/nvmix_compat.f90 src/nvmix.f90
)
gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" -c "${sources[@]}"
mv ./*.o "$build/"
for source in test/*.f90; do
  name=$(basename "${source%.f90}")
  gfortran "${flags[@]}" -J "$build/mod" -I "$build/mod" "$build"/*.o "$source" -o "$build/bin/$name"
  "$build/bin/$name"
done
echo "optimized validation: PASS"
