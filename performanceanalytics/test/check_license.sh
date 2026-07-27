#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
test -s LICENSE
grep -q 'GNU GENERAL PUBLIC LICENSE' LICENSE
grep -q 'license = "GPL-2.0-or-later"' fpm.toml
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file"
  grep -q 'GNU GPL version 2 or, at your option, any later version' "$file"
done < <(find src app example test -type f -name '*.f90' | sort)
echo 'GPL-2.0-or-later source license checks passed.'
