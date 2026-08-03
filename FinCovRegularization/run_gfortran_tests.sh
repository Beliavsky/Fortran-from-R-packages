#!/usr/bin/env sh
set -eu

fc=${FC:-gfortran}
flags=${FFLAGS:--std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g}
build_dir=build-gfortran

rm -rf "$build_dir"
mkdir -p "$build_dir"

sources="
src/fincov_kinds.f90
src/fincov_status.f90
src/fincov_types.f90
src/fincov_utils.f90
src/fincov_rng.f90
src/fincov_linalg.f90
src/fincov_regularization.f90
src/fincov_norms.f90
src/fincov_factor_models.f90
src/fincov_portfolio.f90
src/fincov_cv.f90
src/fincovregularization.f90
"

objects=""
for source in $sources; do
  name=$(basename "$source" .f90)
  "$fc" $flags -J "$build_dir" -I "$build_dir" -c "$source" -o "$build_dir/$name.o"
  objects="$objects $build_dir/$name.o"
done

for source in test/*.f90; do
  name=$(basename "$source" .f90)
  "$fc" $flags -J "$build_dir" -I "$build_dir" "$source" $objects -o "$build_dir/$name"
  "$build_dir/$name"
done

for source in app/*.f90 example/*.f90; do
  name=$(basename "$source" .f90)
  "$fc" $flags -J "$build_dir" -I "$build_dir" "$source" $objects -o "$build_dir/$name"
  "$build_dir/$name"
done
