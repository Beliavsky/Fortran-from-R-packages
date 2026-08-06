@echo off
setlocal
set FFLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Wconversion-extra -Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0
if not exist build-validation mkdir build-validation
cd build-validation

gfortran %FFLAGS% -c ..\src\tvm_kinds.f90 ..\src\tvm_root.f90 ..\src\tvm_interpolation.f90 ..\src\tvm_cashflows.f90 ..\src\tvm_curves.f90 ..\src\tvm.f90 || exit /b 1
for %%F in (..\test\*.f90) do (
  gfortran %FFLAGS% *.o %%F -o %%~nF.exe || exit /b 1
  %%~nF.exe || exit /b 1
)
for %%F in (..\app\*.f90 ..\example\*.f90) do (
  gfortran %FFLAGS% *.o %%F -o %%~nF.exe || exit /b 1
  %%~nF.exe >nul || exit /b 1
)

echo validation: PASS
endlocal
