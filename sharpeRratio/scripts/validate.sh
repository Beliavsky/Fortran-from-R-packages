#!/usr/bin/env bash
set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
cd "$root"

for mode in checked optimized; do
  scripts/build.sh "$mode"
  for executable in "build-$mode"/bin/test_*; do
    "$executable"
  done
  "build-$mode/bin/sharpe_rratio_demo" >/dev/null
  "build-$mode/bin/known_tail_exponent" >/dev/null
done

python - <<'PY'
from pathlib import Path
import tomllib
with Path('fpm.toml').open('rb') as handle:
    manifest = tomllib.load(handle)
assert manifest['name'] == 'sharperratio'
assert manifest['license'] == 'GPL-3.0-only'
print('fpm.toml: PASS')
PY

echo "validation: PASS"
