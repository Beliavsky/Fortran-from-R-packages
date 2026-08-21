#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-validate"
rm -rf "$build"
mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface \
       -fcheck=all -ffpe-trap=invalid,zero,overflow -J"$build/mod" -I"$build/mod")
src=(
  compositions_kinds.f90
  bayesm_kinds.f90
  bayesm_math.f90
  bayesm_linalg.f90
  bayesm_rng.f90
  robustbase_kinds.f90
  robustbase_sort.f90
  robustbase_scale.f90
  robustbase_probability.f90
  robustbase_linalg.f90
  robustbase_covariance.f90
  robustbase_detmcd.f90
  tensora_kinds.f90
  tensora_types.f90
  tensora_index.f90
  tensora_core.f90
  tensora_linalg.f90
  tensora_stats.f90
  tensora.f90
  compositions_linalg.f90
  compositions_geometry.f90
  compositions_distributions.f90
  compositions_zero.f90
  compositions_imputation.f90
  compositions_imputation_cache.f90
  compositions_stats.f90
  compositions_geostat.f90
  compositions_counts.f90
  compositions_gof.f90
  compositions_energy_gof.f90
  compositions_outliers.f90
  compositions_tensor.f90
  compositions.f90
)
objs=()
for s in "${src[@]}"; do
  o="$build/obj/${s%.f90}.o"
  gfortran "${flags[@]}" -c "$root/src/$s" -o "$o"
  objs+=("$o")
done
for t in "$root"/test/*.f90; do
  exe="$build/bin/$(basename "${t%.f90}")"
  gfortran "${flags[@]}" "$t" "${objs[@]}" -llapack -lblas -o "$exe"
  "$exe"
done
gfortran "${flags[@]}" "$root/example/demo_compositions.f90" "${objs[@]}" \
  -llapack -lblas -o "$build/bin/demo_compositions"
"$build/bin/demo_compositions"
