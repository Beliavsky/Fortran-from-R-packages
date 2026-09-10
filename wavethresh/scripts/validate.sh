#!/usr/bin/env sh
set -eu
fpm build
fpm test
fpm run wavethresh_demo
fpm run --example basic_dwt
fpm run --example denoising
