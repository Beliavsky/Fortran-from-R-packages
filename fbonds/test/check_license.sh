#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-or-later
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"
  grep -q 'GNU General Public License, version 2 or later' "$file"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -name '*.f90' -type f | sort)
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"
test -s "$root/LICENSE"
echo 'GPL-2.0-or-later source license checks passed.'
