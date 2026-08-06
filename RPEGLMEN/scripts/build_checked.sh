#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
if command -v fpm >/dev/null 2>&1; then
  fpm test --flag "-std=f2018 -O0 -g -Wall -Wextra -Wpedantic -Wconversion-extra -fcheck=all -fbacktrace -ffree-line-length-none"
else
  make MODE=checked clean test
fi
