#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
missing=0
while IFS= read -r file; do
  if ! grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"; then
    echo "missing GPL-2.0-or-later SPDX identifier: $file" >&2
    missing=1
  fi
  if ! grep -q 'Distributed under the GNU General Public License, version 2 or later.' "$file"; then
    echo "missing GPL version 2-or-later notice: $file" >&2
    missing=1
  fi
  if ! grep -q 'Bernhard Pfaff' "$file"; then
    echo "missing original-author attribution: $file" >&2
    missing=1
  fi
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)

grep -q '^license = "GPL-2.0-or-later"' "$root/fpm.toml" || {
  echo 'fpm.toml license is not GPL-2.0-or-later' >&2
  missing=1
}
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE" || {
  echo 'LICENSE does not contain the GNU GPL text' >&2
  missing=1
}
grep -q 'Version 2, June 1991' "$root/LICENSE" || {
  echo 'LICENSE is not the GPL version 2 text' >&2
  missing=1
}
if [[ $missing -ne 0 ]]; then
  exit 1
fi
echo 'GPL-2.0-or-later source license checks passed.'
