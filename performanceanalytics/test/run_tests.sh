#!/usr/bin/env bash
set -euo pipefail
mode="${1:-debug}"
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
fc="${FC:-gfortran}"
base_flags=(-std=f2018 -ffree-line-length-none -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
if [[ "$mode" == "release" ]]; then
  flags=("${base_flags[@]}" -O2)
else
  flags=("${base_flags[@]}" -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow)
fi
build="build/$mode"
rm -rf "$build"
mkdir -p "$build"
modules=(kinds_mod probability_mod statistics_mod returns_mod drawdown_mod risk_mod capm_mod performance_ratios_mod comoments_mod portfolio_risk_mod rolling_mod rng_mod advanced_moments_mod tail_models_mod analytics_extensions_mod finite_sample_moments_mod csv_mod performanceanalytics_mod)
for mod in "${modules[@]}"; do
  "$fc" "${flags[@]}" -c "src/$mod.f90" -J"$build" -I"$build" -o "$build/$mod.o"
done
"$fc" "${flags[@]}" -c test/test_support_mod.f90 -J"$build" -I"$build" -o "$build/test_support_mod.o"
libobjs=()
for mod in "${modules[@]}"; do libobjs+=("$build/$mod.o"); done
for test_name in test_returns_drawdowns test_risk_ratios test_capm_moments test_rolling_cleaning test_advanced_algorithms test_finite_sample_corrections; do
  "$fc" "${flags[@]}" "test/$test_name.f90" "$build/test_support_mod.o" "${libobjs[@]}" -I"$build" -J"$build" -o "$build/$test_name"
  "$build/$test_name"
done
"$fc" "${flags[@]}" app/demo_performanceanalytics.f90 "${libobjs[@]}" -I"$build" -J"$build" -o "$build/demo_performanceanalytics"
"$fc" "${flags[@]}" app/analyze_csv.f90 "${libobjs[@]}" -I"$build" -J"$build" -o "$build/analyze_csv"
"$fc" "${flags[@]}" example/portfolio_contributions.f90 "${libobjs[@]}" -I"$build" -J"$build" -o "$build/portfolio_contributions"
"$fc" "${flags[@]}" example/advanced_estimators.f90 "${libobjs[@]}" -I"$build" -J"$build" -o "$build/advanced_estimators"
"$fc" "${flags[@]}" example/exact_moment_shrinkage.f90 "${libobjs[@]}" -I"$build" -J"$build" -o "$build/exact_moment_shrinkage"
"$build/demo_performanceanalytics" >/dev/null
"$build/analyze_csv" data/example_returns.csv 1 2 12 >/dev/null
"$build/portfolio_contributions" >/dev/null
"$build/advanced_estimators" >/dev/null
"$build/exact_moment_shrinkage" >/dev/null
printf '%s build, tests, and applications passed.\n' "$mode"
