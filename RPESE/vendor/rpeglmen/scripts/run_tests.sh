#!/usr/bin/env sh
set -eu
cd "$(dirname "$0")/.."
make MODE=checked test
