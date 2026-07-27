@echo off
setlocal
where fpm >nul 2>nul
if errorlevel 1 (
  echo fpm was not found on PATH.
  exit /b 1
)
fpm clean --all
fpm test
if errorlevel 1 exit /b 1
fpm run financialmath_demo
if errorlevel 1 exit /b 1
echo validation: PASS
