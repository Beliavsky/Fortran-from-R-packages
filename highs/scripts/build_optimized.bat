@echo off
setlocal EnableExtensions
for %%I in ("%~dp0..") do set "ROOT=%%~fI"
set "BUILD=%ROOT%\build\optimized"
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%"
if exist "%ROOT%\backend\bin\highs_fortran_bridge.dll" set "HIGHS_FORTRAN_BRIDGE=%ROOT%\backend\bin\highs_fortran_bridge.dll"

gcc -std=c11 -O3 -DNDEBUG -Wall -Wextra -c "%ROOT%\src\highs_dynamic_loader.c" -o "%BUILD%\highs_dynamic_loader.o"
if errorlevel 1 exit /b 1
pushd "%BUILD%"
gfortran -std=f2018 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -J . -I . -c ^
  "%ROOT%\src\highs_kinds.f90" "%ROOT%\src\highs_constants.f90" "%ROOT%\src\highs_sparse.f90" ^
  "%ROOT%\src\highs_types.f90" "%ROOT%\src\highs_c_bindings.f90" "%ROOT%\src\highs_solver_api.f90" ^
  "%ROOT%\src\highs_model_api.f90" "%ROOT%\src\highs.f90"
if errorlevel 1 (popd & exit /b 1)
for %%F in ("%ROOT%\test\*.f90" "%ROOT%\example\*.f90" "%ROOT%\app\*.f90") do (
  gfortran -std=f2018 -O3 -DNDEBUG -Wall -Wextra -Wpedantic -J . -I . "%%~fF" *.o -o "%%~nF.exe"
  if errorlevel 1 (popd & exit /b 1)
  pushd "%ROOT%"
  "%BUILD%\%%~nF.exe"
  if errorlevel 1 (popd & popd & exit /b 1)
  popd
)
popd
