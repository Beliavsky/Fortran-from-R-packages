#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
flags="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -Wno-error=maybe-uninitialized -Wno-error=uninitialized -fimplicit-none -ffree-line-length-none -O3 -march=native"
"$root/scripts/build_common.sh" "$root/build/optimized" "$flags"
"$root/scripts/run_tests.sh" "$root/build/optimized"
