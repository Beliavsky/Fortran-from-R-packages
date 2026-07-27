#!/usr/bin/env sh
# SPDX-License-Identifier: GPL-2.0-or-later
set -eu

fpm build
fpm test
fpm run sde_demo
