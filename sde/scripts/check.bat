@echo off
rem SPDX-License-Identifier: GPL-2.0-or-later
fpm build || exit /b 1
fpm test || exit /b 1
fpm run sde_demo || exit /b 1
