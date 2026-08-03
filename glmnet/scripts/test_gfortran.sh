#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran
mkdir -p build-gfortran/mod build-gfortran/obj build-gfortran/bin

flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace"
sources="src/glmnet_kinds.f90 src/glmnet_status.f90 src/glmnet_types.f90 src/glmnet_utils.f90 src/glmnet_gaussian.f90 src/glmnet_glm.f90 src/glmnet_multinomial.f90 src/glmnet_cox.f90 src/glmnet_predict.f90 src/glmnet_assess.f90 src/glmnet_data.f90 src/glmnet_cv.f90 src/glmnet_relax.f90 src/glmnet_control.f90 src/glmnet_api.f90 src/glmnet_families.f90 src/glmnet.f90"

for source in $sources; do
  object="build-gfortran/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran/mod -Ibuild-gfortran/mod -c "$source" -o "$object"
done
ar rcs build-gfortran/libglmnet.a build-gfortran/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libglmnet.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran/mod "$source" build-gfortran/libglmnet.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'strict GNU Fortran validation: PASS'
