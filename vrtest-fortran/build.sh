#!/usr/bin/env sh
# SPDX-License-Identifier: GPL-2.0-only
set -eu
fpm build
fpm test
fpm run vrtest_demo
