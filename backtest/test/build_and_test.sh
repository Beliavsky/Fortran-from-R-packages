#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
mode=${1:-debug}
fc=${FC:-gfortran}
build="$root/build/$mode"

case "$mode" in
  debug)
    flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace)
    ;;
  release)
    flags=(-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace)
    ;;
  *)
    echo "unknown mode: $mode" >&2
    exit 2
    ;;
esac

"$root/test/check_license.sh"
rm -rf "$build"
mkdir -p "$build"
flags+=("-J$build" "-I$build")

sources=(
  src/backtest_kinds.f90
  src/backtest_math.f90
  src/backtest_bucket.f90
  src/backtest_portfolio.f90
  src/backtest_engine.f90
  src/backtest_summary.f90
)
objects=()
for source in "${sources[@]}"; do
  object="$build/$(basename "${source%.f90}").o"
  "$fc" "${flags[@]}" -c "$root/$source" -o "$object"
  objects+=("$object")
done

for source in test/test_bucket.f90 test/test_portfolio.f90 test/test_engine.f90; do
  exe="$build/$(basename "${source%.f90}")"
  "$fc" "${flags[@]}" "$root/$source" "${objects[@]}" -o "$exe"
  "$exe"
done

for source in app/demo_backtest.f90 app/backtest_csv.f90 example/overlap_example.f90; do
  exe="$build/$(basename "${source%.f90}")"
  "$fc" "${flags[@]}" "$root/$source" "${objects[@]}" -o "$exe"
done

"$build/demo_backtest" >/dev/null
"$build/backtest_csv" "$root/data/example_panel.csv" 2 1 >/dev/null
"$build/backtest_csv" "$root/data/example_panel.csv" 2 2 >/dev/null
"$build/overlap_example" >/dev/null

echo "$mode build, tests, and applications passed."
