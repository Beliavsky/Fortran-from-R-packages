#!/usr/bin/env sh
set -eu

expected='! SPDX-License-Identifier: GPL-3.0-only'

for file in $(find src app example test -type f \( -name '*.f90' -o -name '*.F90' \) | sort); do
    first_line=$(sed -n '1p' "$file")
    if [ "$first_line" != "$expected" ]; then
        echo "Missing GPL-3.0-only SPDX header: $file" >&2
        exit 1
    fi
    if ! grep -q 'Distributed under the GNU General Public License, version 3 only.' "$file"; then
        echo "Missing GPL-3.0-only notice: $file" >&2
        exit 1
    fi
done

if ! grep -q 'GNU GENERAL PUBLIC LICENSE' LICENSE || ! grep -q 'Version 3, 29 June 2007' LICENSE; then
    echo 'LICENSE does not contain the GNU GPL version 3 text.' >&2
    exit 1
fi

if ! grep -q '^license = "GPL-3.0-only"$' fpm.toml; then
    echo 'fpm.toml does not declare GPL-3.0-only.' >&2
    exit 1
fi

echo 'GPL-3.0-only source license checks passed.'
