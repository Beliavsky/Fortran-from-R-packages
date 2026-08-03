#!/usr/bin/env sh
set -eu
ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
"$ROOT/scripts/build_backend.sh"
cd "$ROOT"
if [ "$#" -eq 0 ]; then
  exec fpm run
else
  exec fpm "$@"
fi
