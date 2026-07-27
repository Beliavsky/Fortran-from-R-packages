@echo off
rem SPDX-License-Identifier: GPL-3.0-only

fpm build
if errorlevel 1 exit /b 1

fpm run
if errorlevel 1 exit /b 1

fpm test
if errorlevel 1 exit /b 1

fpm run --example fit_csv -- data\sample_returns.csv
if errorlevel 1 exit /b 1
