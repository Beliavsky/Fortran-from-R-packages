#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
FC=${FC:-gfortran}
case "$mode" in
  debug) flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace) ;;
  release) flags=(-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace) ;;
  *) echo "Usage: $0 {debug|release}" >&2; exit 2 ;;
esac
root=$(cd "$(dirname "$0")" && pwd)
bdir="$root/build/$mode"
rm -rf "$bdir"
mkdir -p "$bdir"
cd "$bdir"
sources=(
  "$root/src/fextremes_kinds.f90"
  "$root/src/fextremes_rng.f90"
  "$root/src/fextremes_stats.f90"
  "$root/src/fextremes_optimize.f90"
  "$root/src/fextremes_distributions.f90"
  "$root/src/fextremes_preprocess.f90"
  "$root/src/fextremes_fit.f90"
  "$root/src/fextremes_extremal_index.f90"
  "$root/src/fextremes_diagnostics.f90"
  "$root/src/fextremes_risk.f90"
  "$root/src/fextremes_metrics.f90"
  "$root/src/fextremes_csv.f90"
)
"$FC" "${flags[@]}" -c "${sources[@]}"
for testsrc in "$root"/test/test_*.f90; do
  name=$(basename "$testsrc" .f90)
  "$FC" "${flags[@]}" "$testsrc" ./*.o -o "$name"
  "./$name"
done
for appsrc in "$root"/app/*.f90 "$root"/example/*.f90; do
  name=$(basename "$appsrc" .f90)
  "$FC" "${flags[@]}" "$appsrc" ./*.o -o "$name"
done
./demo_fextremes >/dev/null
./fit_csv "$root/data/danishClaims.csv" gpd 10 >/dev/null
(cd "$root" && "$bdir/tail_analysis" >/dev/null)
"$root/test/check_license.sh"
echo "$mode build, tests, and applications passed."
