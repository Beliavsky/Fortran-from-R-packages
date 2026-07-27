#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
[[ -s "$root/LICENSE" ]]
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'license = "GPL-2.0-or-later"' "$root/fpm.toml"
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"
  grep -q 'GNU GPL version 2 or any later version' "$file"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)
echo 'GPL-2.0-or-later source license checks passed.'
