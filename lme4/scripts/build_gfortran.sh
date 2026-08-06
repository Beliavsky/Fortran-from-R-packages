#!/usr/bin/env sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
cd "$root"
rm -rf build_gfortran
mkdir -p build_gfortran/mod build_gfortran/obj
flags="-std=f2018 -O2 -Wall -Wextra -Wpedantic"
gfortran $flags -J build_gfortran/mod -I build_gfortran/mod -c \
  dependencies/minqa/src/minqa_module.f90 -o build_gfortran/obj/minqa_module.o
for source in \
  src/lme4_kinds.f90 src/lme4_types.f90 src/lme4_linalg.f90 \
  src/lme4_covariance.f90 src/lme4_quadrature.f90 src/lme4_family.f90 \
  src/lme4_lmm.f90 src/lme4_lmm_pls.f90 src/lme4_glmm.f90 \
  src/lme4_custom_glmm.f90 src/lme4_aghq.f90 src/lme4_aghq_nd.f90 \
  src/lme4_nlmm.f90 src/lme4_simulation.f90 src/lme4_diagnostics.f90 \
  src/lme4_inference.f90 src/lme4_grouped.f90 src/lme4.f90
do
  object=build_gfortran/obj/$(basename "${source%.f90}").o
  gfortran $flags -J build_gfortran/mod -I build_gfortran/mod -c "$source" -o "$object"
done
objects=$(find build_gfortran/obj -name '*.o' -print)
gfortran $flags -I build_gfortran/mod test/test_lme4.f90 $objects -o build_gfortran/test_lme4
./build_gfortran/test_lme4
gfortran $flags -I build_gfortran/mod example/lme4_example.f90 $objects -o build_gfortran/lme4_example
./build_gfortran/lme4_example

gfortran $flags -I build_gfortran/mod example/glmm_extensions_example.f90 $objects -o build_gfortran/glmm_extensions_example
./build_gfortran/glmm_extensions_example

gfortran $flags -I build_gfortran/mod example/advanced_algorithms_example.f90 $objects -o build_gfortran/advanced_algorithms_example
./build_gfortran/advanced_algorithms_example
