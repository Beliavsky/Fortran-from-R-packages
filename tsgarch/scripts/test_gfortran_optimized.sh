#!/usr/bin/env sh
FLAGS="-std=f2018 -Wall -Wextra -Werror -pedantic -O3 -march=native"
export FLAGS
. "$(dirname "$0")/common.sh"
cd "$ROOT"
compile_library
run_tests
