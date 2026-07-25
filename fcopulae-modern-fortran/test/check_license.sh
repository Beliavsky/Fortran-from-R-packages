#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"

grep -q 'GNU GENERAL PUBLIC LICENSE' LICENSE
grep -q 'Version 2, June 1991' LICENSE
grep -q 'license = "GPL-2.0-or-later"' fpm.toml
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file" || {
    echo "Missing GPL-2.0-or-later SPDX header: $file" >&2
    exit 1
  }
  grep -q 'GNU General Public License, version 2 or later' "$file" || {
    echo "Missing GPL notice: $file" >&2
    exit 1
  }
done < <(find src app example test -type f -name '*.f90' | sort)
echo 'GPL-2.0-or-later source license checks passed.'
