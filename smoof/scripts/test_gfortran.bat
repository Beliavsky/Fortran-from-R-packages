@echo off
setlocal
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
for %%F in (src\smoof_kinds.f90 src\smoof_single.f90 src\smoof_multi.f90 src\smoof_cec09.f90 src\smoof_cec2019.f90 src\smoof_ed.f90 src\smoof_nk.f90 src\smoof.f90) do (
  gfortran %FLAGS% -Ibuild_strict -Jbuild_strict -c %%F -o build_strict\%%~nF.o || exit /b 1
)
set OBJS=build_strict\smoof_kinds.o build_strict\smoof_single.o build_strict\smoof_multi.o build_strict\smoof_cec09.o build_strict\smoof_cec2019.o build_strict\smoof_ed.o build_strict\smoof_nk.o build_strict\smoof.o
for %%F in (test\*.f90) do (
  gfortran %FLAGS% -Ibuild_strict %%F %OBJS% -o build_strict\%%~nF.exe || exit /b 1
  build_strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %FLAGS% -Ibuild_strict %%F %OBJS% -o build_strict\%%~nF.exe || exit /b 1
  build_strict\%%~nF.exe || exit /b 1
)
endlocal
