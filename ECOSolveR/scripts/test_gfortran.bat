@echo off
setlocal
cd /d "%~dp0\.."
if exist build-strict rmdir /s /q build-strict
mkdir build-strict
set F=-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

gfortran %F% -J build-strict -I build-strict -c src\ecos_types.f90 src\ecos_sparse.f90 src\ecos_equilibration.f90 src\ecos_linalg.f90 src\ecos_sparse_cones.f90 src\ecos_cones.f90 src\ecos_sparse_solver.f90 src\ecos_solver.f90 src\ecos_bb.f90 src\ecos_api.f90 || exit /b 1

for %%F in (test\*.f90) do (
  gfortran %F% -J build-strict -I build-strict *.o %%F -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
for %%F in (example\*.f90) do (
  gfortran %F% -J build-strict -I build-strict *.o %%F -o build-strict\%%~nF.exe || exit /b 1
  build-strict\%%~nF.exe || exit /b 1
)
del /q *.o *.mod 2>nul
endlocal
