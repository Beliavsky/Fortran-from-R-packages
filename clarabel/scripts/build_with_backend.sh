#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
"$root/scripts/build_backend.sh"
case "$(uname -s)" in
  Darwin) bridge="$root/rust_bridge/bin/libclarabel_fortran_bridge.dylib" ;;
  *) bridge="$root/rust_bridge/bin/libclarabel_fortran_bridge.so" ;;
esac
export CLARABEL_FORTRAN_BRIDGE="$bridge"
cd "$root"
if [ "$#" -eq 0 ]; then
  exec fpm build
else
  exec fpm "$@"
fi
