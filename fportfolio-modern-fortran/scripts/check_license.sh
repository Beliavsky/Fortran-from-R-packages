#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0
while IFS= read -r -d '' f; do
  if ! grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$f"; then
    echo "missing SPDX header: $f" >&2; fail=1
  fi
  if ! grep -qi 'GPL version 2 or later' "$f"; then
    echo "missing GPL notice: $f" >&2; fail=1
  fi
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -print0)
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml" || { echo 'wrong fpm license' >&2; fail=1; }
test -s "$root/LICENSE" || { echo 'missing LICENSE' >&2; fail=1; }
if (( fail )); then exit 1; fi
echo 'GPL-2.0-or-later source license checks passed.'
