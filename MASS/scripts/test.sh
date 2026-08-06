#!/usr/bin/env sh
set -eu
make check
make optimized
make BUILD=build/check MODE=check example
