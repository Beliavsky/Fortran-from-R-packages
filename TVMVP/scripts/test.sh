#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
make clean
make MODE=checked test example
make MODE=optimized test example
