@echo off
setlocal
cd /d %~dp0\..
mingw32-make clean || exit /b 1
mingw32-make check || exit /b 1
mingw32-make release || exit /b 1
mingw32-make MODE=release example || exit /b 1
endlocal
