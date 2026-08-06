#!/usr/bin/env sh
FLAGS="-std=f2018 -Wall -Wextra -Werror -pedantic -fcheck=all -fbacktrace -O0"
export FLAGS
. "$(dirname "$0")/common.sh"
cd "$ROOT"
compile_library
run_tests
