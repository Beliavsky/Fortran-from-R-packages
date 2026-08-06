#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
make clean test
make MODE=optimized clean test
make MODE=optimized example
