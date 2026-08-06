#!/usr/bin/env sh
set -eu
make clean
make MODE=debug test
make clean
make MODE=release test
