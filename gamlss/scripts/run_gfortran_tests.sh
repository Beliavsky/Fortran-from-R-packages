#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
build="$root/build-gfortran"
rm -rf "$build"; mkdir -p "$build/mod" "$build/obj" "$build/bin"
flags=(-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0 -J"$build/mod" -I"$build/mod")
compile() { gfortran "${flags[@]}" -c "$1" -o "$build/obj/$2.o"; }
# splines
for s in splines_kinds splines_linalg splines_core splines_basis splines; do
  compile "$root/vendor/splines-fortran-v0.1.0/src/$s.f90" "$s"
done
# nlme
for s in nlme_kinds nlme_status nlme_types nlme_linalg nlme_correlation nlme_variance \
  nlme_pdmat nlme_optimize nlme_covariance nlme_gls nlme_lme nlme_nonlinear \
  nlme_diagnostics nlme_grouped nlme_models nlme_test_support nlme; do
  compile "$root/vendor/nlme-fortran/src/$s.f90" "$s"
done
# survival
for s in survival_kinds survival_types survival_linalg survival_utils survival_stats \
  survival_nonparametric survival_cox survival_aft survival_pspline survival; do
  compile "$root/vendor/survival-fortran-v0.1.0/src/$s.f90" "$s"
done
# gamlss.dist
for s in gamlss_kinds gamlss_special gamlss_random gamlss_linalg gamlss_optim gamlss_links \
  gamlss_base gamlss_student_t gamlss_continuous gamlss_discrete gamlss_boxcox \
  gamlss_v02_numerics gamlss_continuous_v02 gamlss_discrete_v02 gamlss_continuous_v03 \
  gamlss_discrete_v03 gamlss_flexible_v03 gamlss_fit gamlss_fit_v03 gamlss_dist; do
  compile "$root/vendor/gamlss-dist-fortran-v0.3.0/src/$s.f90" "$s"
done
# gamlss core
for s in gamlss_types gamlss_smoothers gamlss_smoothers_v02 gamlss_additive_v03 gamlss_core gamlss_family_support \
  gamlss_censoring gamlss_random_effects gamlss_random_effects_v03 gamlss_multi_random_v04 gamlss_multi_random_v05 gamlss_correlation_v04 gamlss_correlated_rs_v05 gamlss_copula_v06 gamlss_mvn_v07 gamlss_copula_mixed_v07 gamlss_joint_random_v06 gamlss_joint_random_ghq_v07 gamlss_joint_random_ais_v08 gamlss_marginal_v09 gamlss_pcat gamlss_diagnostics gamlss_diagnostics_v04 gamlss_validation_v05 gamlss_lms \
  gamlss_selection gamlss_model_selection_v02 gamlss_model_selection_v03 gamlss_model_selection_v04 gamlss_bootstrap_v03 gamlss; do
  compile "$root/src/$s.f90" "$s"
done
objs=("$build"/obj/*.o)
for t in test_core test_lms_selection test_v02 test_v03 test_v04 test_v05 test_v06 test_v07 test_v08 test_v09; do
  gfortran "${flags[@]}" "$root/test/$t.f90" "${objs[@]}" -o "$build/bin/$t"
  "$build/bin/$t"
done
for e in basic v02_extended v03_extended v04_extended v05_extended v06_extended v07_extended v08_extended v09_extended; do
  gfortran "${flags[@]}" "$root/example/$e.f90" "${objs[@]}" -o "$build/bin/$e"
  "$build/bin/$e"
done
echo "gamlss-fortran GNU Fortran validation: PASS"
