#!/usr/bin/env bash
set -euo pipefail
mode=${1:-debug}
root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"
./scripts/build.sh "$mode"
./test/check_license.sh
for t in test_scales test_covariance test_regression test_next_batch test_remaining; do
  "build/$mode/bin/$t"
done
"build/$mode/bin/demo_robustbase" > "build/$mode/demo.out"
"build/$mode/bin/fit_csv" data/example_xy.csv mm > "build/$mode/fit_mm.out"
"build/$mode/bin/fit_csv" data/example_xy.csv lts > "build/$mode/fit_lts.out"
"build/$mode/bin/fit_csv" data/example_xy.csv fastlts > "build/$mode/fit_fastlts.out"
"build/$mode/bin/fit_csv" data/example_xy.csv partlts > "build/$mode/fit_partlts.out"
"build/$mode/bin/fit_csv" data/example_xy.csv lar > "build/$mode/fit_lar.out"
"build/$mode/bin/fit_csv" data/example_xy.csv lmrob > "build/$mode/fit_lmrob.out"
"build/$mode/bin/fit_csv" data/example_xy.csv smdm > "build/$mode/fit_smdm.out"
"build/$mode/bin/fit_csv" data/example_binary.csv by > "build/$mode/fit_by.out"
"build/$mode/bin/fit_csv" data/example_binary.csv mqle-binomial > "build/$mode/fit_mqle_binomial.out"
"build/$mode/bin/fit_csv" data/example_binary.csv mt > "build/$mode/fit_mt.out"
"build/$mode/bin/nonlinear_example" > "build/$mode/nonlinear.out"
"build/$mode/bin/next_batch_example" > "build/$mode/next_batch.out"
echo "$mode build, tests, and applications passed."
