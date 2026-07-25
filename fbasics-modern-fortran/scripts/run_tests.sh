#!/usr/bin/env sh
set -eu
mode=${1:-debug}
b="build/$mode"
"$b/test_matrix_stats"
"$b/test_distributions"
"$b/test_inference_interp"
"$b/test_extended_algorithms"
"$b/test_gmm_prewhitening"
"$b/demo_fbasics" >/dev/null
"$b/distribution_example" >/dev/null
"$b/analyze_csv" data/example_returns.csv normal >/dev/null
"$b/analyze_csv" data/example_returns.csv student >/dev/null
"$b/analyze_csv" data/example_returns.csv nig >/dev/null
"$b/analyze_csv" data/example_returns.csv gld >/dev/null
"$b/analyze_csv" data/example_returns.csv fmkl >/dev/null
"$b/analyze_csv" data/example_returns.csv fm5 >/dev/null
"$b/analyze_csv" data/example_returns.csv stable >/dev/null
"$b/analyze_csv" data/example_returns.csv ssd >/dev/null
"$b/analyze_csv" data/example_returns.csv hyp >/dev/null
"$b/analyze_csv" data/example_returns.csv snig >/dev/null
