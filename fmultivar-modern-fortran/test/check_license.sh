#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
fail=0
while IFS= read -r file; do
  if ! grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"; then
    echo "Missing GPL-2.0-or-later SPDX header: $file" >&2
    fail=1
  fi
  if ! grep -q 'GNU General Public License, version 2 or later' "$file"; then
    echo "Missing GPL version 2-or-later notice: $file" >&2
    fail=1
  fi
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)
grep -q '^license = "GPL-2.0-or-later"' "$root/fpm.toml" || { echo 'Incorrect fpm license' >&2; fail=1; }
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE" || { echo 'LICENSE does not contain GPL text' >&2; fail=1; }
if (( fail != 0 )); then exit 1; fi
echo 'GPL-2.0-or-later source license checks passed.'
