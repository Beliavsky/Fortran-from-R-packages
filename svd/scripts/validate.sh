#!/usr/bin/env sh
set -eu
fpm build
fpm test
fpm run svd_demo
fpm run --example dense_svd
fpm run --example matrix_free_svd
python scripts/audit.py
fpm clean --all
