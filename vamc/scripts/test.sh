#!/usr/bin/env sh
set -eu
make clean
make check
make clean
make optimized
make MODE=optimized FLAGS='-std=f2018 -Wall -Wextra -Werror -pedantic -O3' example
