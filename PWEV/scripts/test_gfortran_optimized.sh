#!/bin/sh
set -eu
cd "$(dirname "$0")/.."
rm -rf build-gfortran-optimized
mkdir -p build-gfortran-optimized/mod build-gfortran-optimized/obj build-gfortran-optimized/bin
flags="-std=f2018 -Wall -Wextra -Wpedantic -Werror -Wimplicit-interface -O3"
maxlik_sources="vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_kinds.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_status.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_types.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_linalg.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_random.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_evaluation.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_solvers.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_inference.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_utilities.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik_api.f90 vendor/rumidas-fortran/vendor/maxLik-fortran/src/maxlik.f90"
rumidas_sources="vendor/rumidas-fortran/src/rumidas_kinds.f90 vendor/rumidas-fortran/src/rumidas_status.f90 vendor/rumidas-fortran/src/rumidas_types.f90 vendor/rumidas-fortran/src/rumidas_weights.f90 vendor/rumidas-fortran/src/rumidas_statistics.f90 vendor/rumidas-fortran/src/rumidas_garch_midas.f90 vendor/rumidas-fortran/src/rumidas_mem.f90 vendor/rumidas-fortran/src/rumidas_fit.f90 vendor/rumidas-fortran/src/rumidas_forecast.f90 vendor/rumidas-fortran/src/rumidas.f90"
rugarch_sources="vendor/rugarch-modern-fortran/src/rugarch_kinds.f90 vendor/rugarch-modern-fortran/src/rugarch_math.f90 vendor/rugarch-modern-fortran/src/rugarch_rng.f90 vendor/rugarch-modern-fortran/src/rugarch_gh.f90 vendor/rugarch-modern-fortran/src/rugarch_distributions.f90 vendor/rugarch-modern-fortran/src/rugarch_optimizer.f90 vendor/rugarch-modern-fortran/src/rugarch_linalg.f90 vendor/rugarch-modern-fortran/src/rugarch_resampling.f90 vendor/rugarch-modern-fortran/src/rugarch_types.f90 vendor/rugarch-modern-fortran/src/rugarch_models.f90 vendor/rugarch-modern-fortran/src/rugarch_fit.f90 vendor/rugarch-modern-fortran/src/rugarch_risk.f90 vendor/rugarch-modern-fortran/src/rugarch_arfima.f90 vendor/rugarch-modern-fortran/src/rugarch_backtests.f90 vendor/rugarch-modern-fortran/src/rugarch_inference.f90 vendor/rugarch-modern-fortran/src/rugarch_evaluation.f90 vendor/rugarch-modern-fortran/src/rugarch_workflows.f90 vendor/rugarch-modern-fortran/src/rugarch_complete.f90 vendor/rugarch-modern-fortran/src/rugarch.f90"
pwev_sources="src/pwev_kinds.f90 src/pwev_status.f90 src/pwev_types.f90 src/pwev_metrics.f90 src/pwev_pso.f90 src/pwev_models.f90 src/pwev_core.f90 src/pwev.f90"
for source in $maxlik_sources $rumidas_sources $rugarch_sources $pwev_sources; do
  object="build-gfortran-optimized/obj/$(basename "$source" .f90).o"
  gfortran $flags -Jbuild-gfortran-optimized/mod -Ibuild-gfortran-optimized/mod -c "$source" -o "$object"
done
ar rcs build-gfortran-optimized/libpwev.a build-gfortran-optimized/obj/*.o
for source in test/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libpwev.a -o "$executable"
  "$executable"
done
for source in example/*.f90 app/*.f90; do
  executable="build-gfortran-optimized/bin/$(basename "$source" .f90)"
  gfortran $flags -Ibuild-gfortran-optimized/mod "$source" build-gfortran-optimized/libpwev.a -o "$executable"
  "$executable" >/dev/null
done
printf '%s\n' 'optimized GNU Fortran validation: PASS'
