@echo off
setlocal
cd /d %~dp0\..
if exist build_strict rmdir /s /q build_strict
mkdir build_strict\mod build_strict\obj build_strict\bin
set F=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
set S=src\coneproj_kinds.f90 src\coneproj_types.f90 src\coneproj_linalg.f90 src\coneproj_stats.f90 src\coneproj_core.f90 src\coneproj_shape.f90 src\coneproj_regression.f90 src\coneproj.f90
gfortran %F% -Jbuild_strict\mod -Ibuild_strict\mod -c %S% || exit /b 1
move /y *.o build_strict\obj\ >nul
for %%f in (test\*.f90) do (
  gfortran %F% -Ibuild_strict\mod %%f build_strict\obj\*.o -o build_strict\bin\%%~nf.exe || exit /b 1
  build_strict\bin\%%~nf.exe || exit /b 1
)
for %%f in (example\*.f90) do (
  gfortran %F% -Ibuild_strict\mod %%f build_strict\obj\*.o -o build_strict\bin\%%~nf.exe || exit /b 1
  build_strict\bin\%%~nf.exe || exit /b 1
)
endlocal
