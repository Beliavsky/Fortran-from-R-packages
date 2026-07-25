#!/usr/bin/env bash
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
if ! grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"; then
  echo "LICENSE does not contain the GPL text" >&2
  exit 1
fi
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-or-later' "$file" || {
    echo "missing GPL-2.0-or-later SPDX header: $file" >&2
    exit 1
  }
  grep -qi 'GNU General Public License version 2 or later' "$file" || {
    echo "missing GPL notice: $file" >&2
    exit 1
  }
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f -name '*.f90' | sort)
echo "GPL-2.0-or-later source license checks passed."
