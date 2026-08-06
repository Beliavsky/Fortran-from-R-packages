#!/usr/bin/env sh
set -eu
make clean
make MODE=optimized test examples
