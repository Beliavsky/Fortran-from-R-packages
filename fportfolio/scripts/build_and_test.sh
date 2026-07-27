#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
case "$mode" in
  debug) flags=(-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -ffree-line-length-none) ;;
  release) flags=(-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace -ffree-line-length-none) ;;
  *) echo "usage: $0 debug|release" >&2; exit 2 ;;
esac
b="build/$mode"
rm -rf "$b"; mkdir -p "$b/mod" "$b/obj" "$b/bin"
sources=(
  fportfolio_kinds fportfolio_types fportfolio_probability fportfolio_linalg
  fportfolio_statistics fportfolio_risk fportfolio_optimization
  fportfolio_backtest fportfolio_monitor fportfolio_io
)
objects=()
for name in "${sources[@]}"; do
  obj="$b/obj/$name.o"
  gfortran "${flags[@]}" -J"$b/mod" -I"$b/mod" -c "src/$name.f90" -o "$obj"
  objects+=("$obj")
done
link_program() {
  local src=$1 out=$2
  gfortran "${flags[@]}" -J"$b/mod" -I"$b/mod" "$src" "${objects[@]}" -llapack -lblas -o "$b/bin/$out"
}
link_program test/test_statistics_risk.f90 test_statistics_risk
link_program test/test_optimization.f90 test_optimization
link_program test/test_backtest_monitor.f90 test_backtest_monitor
link_program example/portfolio_demo.f90 portfolio_demo
link_program app/fit_csv.f90 fit_csv
link_program app/backtest_csv.f90 backtest_csv
"$b/bin/test_statistics_risk"
"$b/bin/test_optimization"
"$b/bin/test_backtest_monitor"
"$b/bin/portfolio_demo" >/dev/null
"$b/bin/fit_csv" data/example_returns.csv minvariance >/dev/null
"$b/bin/fit_csv" data/example_returns.csv tangency 0.0 >/dev/null
"$b/bin/fit_csv" data/example_returns.csv riskparity >/dev/null
"$b/bin/fit_csv" data/example_returns.csv maxdiv >/dev/null
"$b/bin/fit_csv" data/example_returns.csv mad >/dev/null
"$b/bin/fit_csv" data/example_returns.csv cvar 0.05 >/dev/null
"$b/bin/fit_csv" data/example_returns.csv efficient 0.0005 >/dev/null
"$b/bin/fit_csv" data/example_returns.csv maxreturn >/dev/null
"$b/bin/backtest_csv" data/example_returns.csv minvariance 60 20 >/dev/null
"$b/bin/backtest_csv" data/example_returns.csv tangency 60 20 >/dev/null
"$b/bin/backtest_csv" data/example_returns.csv riskparity 60 20 >/dev/null
"$root/scripts/check_license.sh"
echo "$mode build, tests, and applications passed."
