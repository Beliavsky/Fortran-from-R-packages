#!/usr/bin/env sh
set -eu

root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"

rm -rf build-validation
mkdir -p build-validation/mod build-validation/obj build-validation/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0 -g -Jbuild-validation/mod -Ibuild-validation/mod"
objects=""

for source in \
  src/blmodel_kinds.f90 \
  src/blmodel_types.f90 \
  src/blmodel_linalg.f90 \
  src/blmodel_utils.f90 \
  src/blmodel_distributions.f90 \
  src/blmodel_equilibrium.f90 \
  src/blmodel_posterior.f90 \
  src/blmodel.f90
do
  object="build-validation/obj/$(basename "${source%.f90}").o"
  gfortran $flags -c "$source" -o "$object"
  objects="$objects $object"
done

gfortran $flags -c test/test_support.f90 -o build-validation/obj/test_support.o

for source in test/test_distributions.f90 test/test_equilibrium.f90 test/test_posterior.f90 test/test_blmodel.f90
do
  name=$(basename "${source%.f90}")
  gfortran $flags "$source" build-validation/obj/test_support.o $objects -o "build-validation/bin/$name"
  "build-validation/bin/$name"
done

for source in app/blmodel_demo.f90 example/basic_black_litterman.f90 example/view_distributions.f90
do
  name=$(basename "${source%.f90}")
  gfortran $flags "$source" $objects -o "build-validation/bin/$name"
  "build-validation/bin/$name" >/dev/null
done

printf '%s\n' 'validation: PASS'
