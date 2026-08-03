#!/bin/sh
set -eu

cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin

flags="-std=f2018 -O3 -Wall -Wextra -Wpedantic -Werror -fbacktrace"
sources="src/glmnet_kinds.f90 src/glmnet_status.f90 src/glmnet_types.f90 src/glmnet_utils.f90 src/glmnet_gaussian.f90 src/glmnet_glm.f90 src/glmnet_multinomial.f90 src/glmnet_cox.f90 src/glmnet_predict.f90 src/glmnet_assess.f90 src/glmnet_data.f90 src/glmnet_cv.f90 src/glmnet_relax.f90 src/glmnet_control.f90 src/glmnet_api.f90 src/glmnet_families.f90 src/glmnet.f90"

for source in $sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/libglmnet.a build-gfortran-optimized/obj/*.o

for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libglmnet.a -o "$executable"
  "$executable"
done

for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libglmnet.a -o "$executable"
  "$executable" >/dev/null
done

printf '%s\n' 'optimized GNU Fortran validation: PASS'
