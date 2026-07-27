#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

fc=${FC:-gfortran}
build_dir=${BUILD_DIR:-build-validation}
flags="-std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0 -J$build_dir -I$build_dir"

rm -rf "$build_dir"
mkdir -p "$build_dir"

for source in \
  src/sharper_kinds.f90 \
  src/sharper_types.f90 \
  src/sharper_math.f90 \
  src/sharper_linalg.f90 \
  src/sharper_distributions.f90 \
  src/sharper_estimation.f90 \
  src/sharper_inference.f90 \
  src/sharper_tests.f90 \
  src/sharper.f90
do
  object="$build_dir/$(basename "${source%.f90}").o"
  "$fc" $flags -c "$source" -o "$object"
done

objects=""
for name in sharper_kinds sharper_types sharper_math sharper_linalg \
  sharper_distributions sharper_estimation sharper_inference sharper_tests sharper
do
  objects="$objects $build_dir/$name.o"
done

for source in test/*.f90 app/*.f90 example/*.f90
do
  target="$build_dir/$(basename "${source%.f90}")"
  "$fc" $flags "$source" $objects -o "$target"
  "$target"
done
