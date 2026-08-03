#!/usr/bin/env sh
set -eu

fc=${FC:-gfortran}
flags=${FFLAGS:--std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -O0 -g}
build_dir=build-gfortran

rm -rf "$build_dir"
mkdir -p "$build_dir"

$fc $flags -J "$build_dir" -I "$build_dir" -c \
  src/fincal_kinds.f90 \
  src/fincal_status.f90 \
  src/fincal_types.f90 \
  src/fincal_rates.f90 \
  src/fincal_tvm.f90 \
  src/fincal_ratios.f90 \
  src/fincal_statistics.f90 \
  src/fincal_accounting.f90 \
  src/fincal.f90

objects="$build_dir/fincal_kinds.o $build_dir/fincal_status.o $build_dir/fincal_types.o $build_dir/fincal_rates.o $build_dir/fincal_tvm.o $build_dir/fincal_ratios.o $build_dir/fincal_statistics.o $build_dir/fincal_accounting.o $build_dir/fincal.o"

# Some compilers place objects in the current directory when -c is used.
for object in fincal_kinds fincal_status fincal_types fincal_rates fincal_tvm fincal_ratios fincal_statistics fincal_accounting fincal; do
  if [ -f "$object.o" ]; then mv "$object.o" "$build_dir/$object.o"; fi
done

for source in test/*.f90; do
  name=$(basename "$source" .f90)
  $fc $flags -J "$build_dir" -I "$build_dir" "$source" $objects -o "$build_dir/$name"
  "$build_dir/$name"
done

for source in app/*.f90 example/*.f90; do
  name=$(basename "$source" .f90)
  $fc $flags -J "$build_dir" -I "$build_dir" "$source" $objects -o "$build_dir/$name"
  "$build_dir/$name"
done
