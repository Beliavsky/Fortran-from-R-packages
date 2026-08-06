#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
make clean
make manifest
make check
make optimized
make demo
