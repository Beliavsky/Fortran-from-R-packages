#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
make clean
make check
