@echo off
setlocal
cd /d %~dp0\..
mingw32-make check || exit /b 1
mingw32-make optimized || exit /b 1
mingw32-make example || exit /b 1
endlocal
