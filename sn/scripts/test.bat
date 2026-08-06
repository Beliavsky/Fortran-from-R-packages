@echo off
setlocal
cd /d %~dp0\..
if exist build rmdir /s /q build
make manifest || exit /b 1
make check || exit /b 1
make optimized || exit /b 1
make demo || exit /b 1
endlocal
