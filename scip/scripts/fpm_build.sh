#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd "$(dirname "$0")/.." && pwd)
if [ ! -f "$ROOT/vendor/lib/libscipfortran_backend.a" ]; then
    "$ROOT/scripts/build_vendor.sh"
fi
export LIBRARY_PATH="$ROOT/vendor/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
cd "$ROOT"
exec fpm build "$@"
