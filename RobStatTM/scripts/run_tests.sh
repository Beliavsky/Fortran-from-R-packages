#!/usr/bin/env bash
set -euo pipefail
build_dir="${1:?build directory required}"
for executable in "$build_dir"/bin/test_*; do
  echo "running $(basename "$executable")"
  "$executable"
done
