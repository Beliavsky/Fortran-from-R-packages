@echo off
mingw32-make clean || exit /b 1
mingw32-make MODE=checked all || exit /b 1
