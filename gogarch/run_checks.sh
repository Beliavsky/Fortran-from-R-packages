#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")" && pwd)
fc=${FC:-gfortran}
libs='-llapack -lblas'
sources=(
  gogarch_kinds.f90
  gogarch_linalg.f90
  gogarch_optimizer.f90
  gogarch_rng.f90
  gogarch_distributions.f90
  gogarch_types.f90
  gogarch_univariate.f90
  gogarch_orthogonal.f90
  gogarch_core.f90
  gogarch_ica.f90
  gogarch_model.f90
  gogarch_estimators.f90
  gogarch.f90
)
objects=(
  gogarch.o gogarch_estimators.o gogarch_model.o gogarch_ica.o
  gogarch_core.o gogarch_orthogonal.o gogarch_univariate.o
  gogarch_types.o gogarch_distributions.o gogarch_rng.o
  gogarch_optimizer.o gogarch_linalg.o gogarch_kinds.o
)

"$root/test/check_license.sh"

build_and_run() {
  local name=$1
  shift
  local flags=("$@")
  local dir="$root/build/$name"
  rm -rf "$dir"
  mkdir -p "$dir"
  cd "$dir"
  for source in "${sources[@]}"; do
    "$fc" "${flags[@]}" -J. -I. -c "$root/src/$source"
  done
  "$fc" "${flags[@]}" -J. -I. -c "$root/test/test_helpers.f90"
  "$fc" "${flags[@]}" -J. -I. "$root/test/test_core.f90" test_helpers.o "${objects[@]}" $libs -o test_core
  "$fc" "${flags[@]}" -J. -I. "$root/test/test_univariate_extended.f90" test_helpers.o \
    "${objects[@]}" $libs -o test_univariate_extended
  "$fc" "${flags[@]}" -J. -I. "$root/test/test_estimators.f90" test_helpers.o \
    "${objects[@]}" $libs -o test_estimators
  "$fc" "${flags[@]}" -J. -I. "$root/test/test_gogarch_extensions.f90" test_helpers.o \
    "${objects[@]}" $libs -o test_gogarch_extensions
  "$fc" "${flags[@]}" -J. -I. "$root/app/demo_gogarch.f90" "${objects[@]}" $libs -o demo_gogarch
  "$fc" "${flags[@]}" -J. -I. "$root/example/fit_csv.f90" "${objects[@]}" $libs -o fit_csv
  ./test_core
  ./test_univariate_extended
  ./test_estimators
  ./test_gogarch_extensions
  ./demo_gogarch > demo_output.txt
  ./fit_csv "$root/data/sample_returns.csv" ica > fit_csv_ica_output.txt
  ./fit_csv "$root/data/sample_returns.csv" mm > fit_csv_mm_output.txt
  ./fit_csv "$root/data/sample_returns.csv" nls > fit_csv_nls_output.txt
  ./fit_csv "$root/data/sample_returns.csv" ml > fit_csv_ml_output.txt
  ./fit_csv "$root/data/sample_returns.csv" ica garch std 2 0 1 > fit_csv_garch21_std_output.txt
  ./fit_csv "$root/data/sample_returns.csv" ica aparch sstd 1 1 1 > fit_csv_aparch_sstd_output.txt
  echo "$name build and execution checks passed."
}

common=(-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
build_and_run debug "${common[@]}" -O0 -g -fcheck=all
build_and_run optimized "${common[@]}" -O2
