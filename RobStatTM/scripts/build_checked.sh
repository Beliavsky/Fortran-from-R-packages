#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
flags="-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace -fimplicit-none -ffree-line-length-none -O0 -g"
"$root/scripts/build_common.sh" "$root/build/checked" "$flags"
"$root/scripts/run_tests.sh" "$root/build/checked"
