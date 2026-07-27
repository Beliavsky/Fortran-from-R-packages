#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
mode=${1:-debug}
rm -rf build
mkdir -p build/mod build/obj build/bin
common=(-std=f2018 -Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror -Jbuild/mod -Ibuild/mod)
if [[ "$mode" == "optimized" ]]; then
  flags=("${common[@]}" -O2)
else
  flags=("${common[@]}" -O0 -fcheck=all -fbacktrace)
fi
sources=(
  src/hdshop_kinds.f90
  src/hdshop_linalg.f90
  src/hdshop_stats.f90
  src/hdshop_shrinkage.f90
  src/hdshop_portfolio.f90
  src/hdshop_inference.f90
  src/hdshop_random.f90
  src/hdshop_formulas.f90
  src/hdshop_compat.f90
  src/hdshop.f90
)
objects=()
for source in "${sources[@]}"; do
  object="build/obj/$(basename "${source%.f90}").o"
  gfortran "${flags[@]}" -c "$source" -o "$object"
  objects+=("$object")
done
for source in test/*.f90; do
  exe="build/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$exe"
  "$exe"
done
for source in app/*.f90 example/*.f90; do
  exe="build/bin/$(basename "${source%.f90}")"
  gfortran "${flags[@]}" "$source" "${objects[@]}" -o "$exe"
  "$exe" >/dev/null
done
echo "validation ($mode): PASS"
