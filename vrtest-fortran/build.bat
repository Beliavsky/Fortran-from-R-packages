@echo off
rem SPDX-License-Identifier: GPL-2.0-only
fpm build || exit /b 1
fpm test || exit /b 1
fpm run vrtest_demo || exit /b 1
