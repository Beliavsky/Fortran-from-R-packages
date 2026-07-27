#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-2.0-only
# Copyright (C) 2010-2023 Yang Lu and David Kane
# Copyright (C) 2026 Modern Fortran translation contributors
# This program is free software under GNU GPL version 2 only.
set -euo pipefail
root=$(cd "$(dirname "$0")/.." && pwd)
while IFS= read -r file; do
  grep -q 'SPDX-License-Identifier: GPL-2.0-only' "$file"
  grep -q 'GNU General Public License version 2 only\|GNU GPL version 2 only' "$file"
done < <(find "$root/src" "$root/app" "$root/example" "$root/test" -type f \( -name '*.f90' -o -name '*.sh' \) | sort)
grep -q 'GNU GENERAL PUBLIC LICENSE' "$root/LICENSE"
grep -q 'license = "GPL-2.0-only"' "$root/fpm.toml"
echo 'GPL-2.0-only source license checks passed.'
