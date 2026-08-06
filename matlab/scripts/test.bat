@echo off
setlocal
mingw32-make clean || exit /b 1
mingw32-make MODE=debug test || exit /b 1
mingw32-make clean || exit /b 1
mingw32-make MODE=release test || exit /b 1
endlocal
