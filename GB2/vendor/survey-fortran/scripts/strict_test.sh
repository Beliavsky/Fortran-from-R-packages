#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD="${TMPDIR:-/tmp}/survey-fortran-strict"
FC="${FC:-gfortran}"
STRICT=(-std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all -ffree-line-length-none)

rm -rf "$BUILD"
mkdir -p "$BUILD"
cd "$BUILD"

# Vendored Powell/minqa code uses exact floating-point sentinels and is kept
# close to the supplied translation; do not promote those warnings to errors.
"$FC" -std=f2018 -O2 -fcheck=all -ffree-line-length-none -w -c \
    "$ROOT/vendor/minqa/src/minqa_module.f90"

for f in numderiv_kinds.f90 numderiv_types.f90 numderiv_callbacks.f90 \
         numderiv_core.f90 numderiv.f90; do
    "$FC" "${STRICT[@]}" -c "$ROOT/vendor/numDeriv-fortran/src/$f"
done

for f in survival_kinds.f90 survival_types.f90 survival_linalg.f90 \
         splines_kinds.f90 splines_linalg.f90 splines_core.f90 \
         splines_basis.f90 splines.f90 survival_nonparametric.f90 \
         survival_cox.f90 survival_aft.f90 survival_stats.f90 \
         survival_utils.f90 survival_pspline.f90 survival.f90; do
    "$FC" "${STRICT[@]}" -c "$ROOT/vendor/survival-fortran/src/$f"
done

for f in survey_kinds.f90 survey_types.f90 survey_linalg.f90 survey_design.f90 \
         survey_taylor.f90 survey_estimators.f90 survey_replicates.f90 \
         survey_calibration.f90 survey_quantiles.f90 survey_glm.f90 \
         survey_survival.f90 survey_multivariate.f90 survey_inference.f90 \
         survey_pps.f90 survey_special.f90 survey_chisq.f90 survey_ivreg.f90 \
         survey_mle.f90 survey_nls.f90 survey_ordinal.f90 survey_loglinear.f90 survey_factor.f90 survey_phase.f90 survey_multiframe.f90 survey_mrb.f90 survey_score.f90 survey.f90; do
    "$FC" "${STRICT[@]}" -c "$ROOT/src/$f"
done

for t in "$ROOT"/test/*.f90; do
    exe="$(basename "${t%.f90}")"
    "$FC" "${STRICT[@]}" "$t" ./*.o -o "$exe"
    "./$exe"
done

"$FC" "${STRICT[@]}" "$ROOT/example/basic_survey.f90" ./*.o -o basic_survey
./basic_survey
