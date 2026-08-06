#!/usr/bin/env bash
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
build_dir="${1:?build directory required}"
flags="${2:?compiler flags required}"
fc="${FC:-gfortran}"

rm -rf "$build_dir"
mkdir -p "$build_dir/mod" "$build_dir/obj" "$build_dir/bin"

compile_source() {
  local source="$1"
  local object="$build_dir/obj/$(basename "${source%.f90}").o"
  "$fc" $flags -J"$build_dir/mod" -I"$build_dir/mod" -c "$source" -o "$object"
}

robustbase="$root/vendor/robustbase-modern-fortran/src"
for unit in \
  robustbase_kinds robustbase_sort robustbase_linalg robustbase_probability \
  robustbase_pca robustbase_scale robustbase_psi robustbase_medcouple \
  robustbase_covariance robustbase_detmcd robustbase_regression robustbase_fastlts \
  robustbase_fast_algorithms robustbase_bylogreg robustbase_lmrob; do
  compile_source "$robustbase/$unit.f90"
done

rrcov="$root/vendor/rrcov/src"
for unit in \
  rrcov_kinds rrcov_types rrcov_random rrcov_sort rrcov_linalg rrcov_stats \
  rrcov_robust rrcov_pca; do
  compile_source "$rrcov/$unit.f90"
done

for unit in \
  robstattm_kinds robstattm_types robstattm_utils robstattm_psi \
  robstattm_regression robstattm_logistic robstattm_multivariate robstattm_pca \
  robstattm_compat robstattm; do
  compile_source "$root/src/$unit.f90"
done

objects=("$build_dir"/obj/*.o)
for source in "$root"/test/test_*.f90; do
  executable="$build_dir/bin/$(basename "${source%.f90}")"
  "$fc" $flags -I"$build_dir/mod" "$source" "${objects[@]}" -llapack -lblas -o "$executable"
done

if compgen -G "$root/example/*.f90" > /dev/null; then
  for source in "$root"/example/*.f90; do
    executable="$build_dir/bin/$(basename "${source%.f90}")"
    "$fc" $flags -I"$build_dir/mod" "$source" "${objects[@]}" -llapack -lblas -o "$executable"
  done
fi
