#!/usr/bin/env sh
set -eu
make check
make optimized
make example
