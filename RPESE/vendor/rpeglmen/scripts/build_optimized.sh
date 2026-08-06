#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
if command -v fpm >/dev/null 2>&1; then
  fpm test --profile release --flag "-std=f2018 -O3 -Wall -Wextra -Wpedantic -ffree-line-length-none"
else
  make MODE=optimized clean test
fi
