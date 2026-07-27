#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

status=0
while IFS= read -r file; do
  if ! grep -Fq 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"; then
    echo "missing GPL-2.0-or-later SPDX header: $file" >&2
    status=1
  fi
  if ! grep -Fq 'Distributed under the GNU General Public License, version 2 or later.' "$file"; then
    echo "missing GPL version 2-or-later notice: $file" >&2
    status=1
  fi
done < <(find src app example test -type f -name '*.f90' | sort)

grep -Fq 'license = "GPL-2.0-or-later"' fpm.toml || { echo 'fpm.toml license mismatch' >&2; status=1; }
grep -Fq 'GNU GENERAL PUBLIC LICENSE' LICENSE || { echo 'LICENSE does not contain GNU GPL text' >&2; status=1; }
grep -Fq 'License: GPL (>= 2)' reference/DESCRIPTION.original || { echo 'original DESCRIPTION license not preserved in reference metadata' >&2; status=1; }

if [ "$status" -ne 0 ]; then
  exit "$status"
fi
echo 'GPL-2.0-or-later source license checks passed.'
