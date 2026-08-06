#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
./scripts/build_checked.sh
./scripts/build_optimized.sh
