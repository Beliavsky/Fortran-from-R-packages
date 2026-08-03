@echo off
setlocal
cd /d %~dp0\..
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran
cd build-gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow
set SRC=..\src\fastcluster_kinds.f90 ..\src\fastcluster_types.f90 ..\src\fastcluster_distances.f90 ..\src\fastcluster_core.f90 ..\src\cluster_types.f90 ..\src\cluster_linalg.f90 ..\src\cluster_daisy.f90 ..\src\cluster_partition.f90 ..\src\cluster_hierarchy.f90 ..\src\cluster_diagnostics.f90 ..\src\cluster_ellipsoid.f90 ..\src\cluster.f90
gfortran %FLAGS% -J . -c %SRC% || exit /b 1
ar rcs libcluster.a *.o || exit /b 1
for %%F in (..\test\*.f90) do (
  gfortran %FLAGS% -I . "%%F" libcluster.a -o "%%~nF.exe" || exit /b 1
  "%%~nF.exe" || exit /b 1
)
endlocal
