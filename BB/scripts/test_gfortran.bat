@echo off
setlocal
if exist build-gfortran rmdir /s /q build-gfortran
mkdir build-gfortran\mod
mkdir build-gfortran\obj
mkdir build-gfortran\bin
pushd build-gfortran\obj
gfortran -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -J..\mod -I..\mod -c ..\..\src\quadprog_kinds.f90 ..\..\src\quadprog_core.f90 ..\..\src\quadprog.f90 ..\..\src\bb_kinds.f90 ..\..\src\bb_interfaces.f90 ..\..\src\bb_types.f90 ..\..\src\bb_projection.f90 ..\..\src\bb_spg.f90 ..\..\src\bb_aux_optim.f90 ..\..\src\bb_nonlinear.f90 ..\..\src\bb_drivers.f90 ..\..\src\bb.f90
if errorlevel 1 exit /b 1
ar rcs ..\libbb.a *.o
if errorlevel 1 exit /b 1
popd
for %%F in (test\*.f90) do (
  gfortran -std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -Ibuild-gfortran\mod %%F build-gfortran\libbb.a -o build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
  build-gfortran\bin\%%~nF.exe
  if errorlevel 1 exit /b 1
)
endlocal
