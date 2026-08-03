@echo off
setlocal EnableExtensions EnableDelayedExpansion
for %%I in ("%~dp0..") do set "ROOT=%%~fI"
set "OSQP_SRC=%ROOT%\upstream\osqp-r\src\osqp_sources"
set "QDLDL_SRC=%ROOT%\upstream\osqp-r\src\qdldl_sources"
set "OSQP_BUILD=%ROOT%\backend\osqp-build"
set "OSQP_INSTALL=%ROOT%\backend\osqp-install"
set "BRIDGE_BUILD=%ROOT%\backend\bridge-build"
set "OUT_DIR=%ROOT%\backend\bin"

where cmake >nul 2>nul || (
  echo cmake is required to build the bundled OSQP backend.
  exit /b 1
)
where gcc >nul 2>nul || (
  echo gcc was not found. Use the same MinGW toolchain as gfortran.
  exit /b 1
)

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"
copy /y "%ROOT%\bridge\generic_printing.h" "%OSQP_SRC%\include\private\printing.h" >nul

set "GENERATOR="
where ninja >nul 2>nul && set "GENERATOR=-G Ninja"
if not defined GENERATOR (
  where mingw32-make >nul 2>nul && set GENERATOR=-G "MinGW Makefiles"
)
if not defined GENERATOR (
  echo Ninja or mingw32-make is required for the MinGW backend build.
  exit /b 1
)

cmake -S "%OSQP_SRC%" -B "%OSQP_BUILD%" %GENERATOR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc ^
  -DOSQP_VERSION=1.0.0 ^
  -DOSQP_BUILD_SHARED_LIB=OFF ^
  -DOSQP_BUILD_STATIC_LIB=ON ^
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON ^
  -DOSQP_ALGEBRA_BACKEND=builtin ^
  -DOSQP_ENABLE_PRINTING=ON ^
  -DOSQP_ENABLE_PROFILING=ON ^
  -DOSQP_ENABLE_INTERRUPT=ON ^
  -DOSQP_USE_LONG=OFF ^
  -DOSQP_USE_FLOAT=OFF ^
  -DOSQP_BUILD_DEMO_EXE=OFF ^
  -DOSQP_BUILD_UNITTESTS=OFF ^
  -DOSQP_CODEGEN=OFF ^
  -DOSQP_ENABLE_DERIVATIVES=OFF ^
  -DFETCHCONTENT_SOURCE_DIR_QDLDL="%QDLDL_SRC%" ^
  -DCMAKE_INSTALL_PREFIX="%OSQP_INSTALL%"
if errorlevel 1 exit /b 1
cmake --build "%OSQP_BUILD%" --target install --parallel
if errorlevel 1 exit /b 1

cmake -S "%ROOT%\bridge" -B "%BRIDGE_BUILD%" %GENERATOR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc ^
  -DOSQP_ROOT="%OSQP_INSTALL%" ^
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="%OUT_DIR%"
if errorlevel 1 exit /b 1
cmake --build "%BRIDGE_BUILD%" --config Release --parallel
if errorlevel 1 exit /b 1

set "BRIDGE_DLL="
for /r "%BRIDGE_BUILD%" %%F in (osqp_fortran_bridge.dll) do set "BRIDGE_DLL=%%~fF"
if not defined BRIDGE_DLL if exist "%OUT_DIR%\osqp_fortran_bridge.dll" set "BRIDGE_DLL=%OUT_DIR%\osqp_fortran_bridge.dll"
if not defined BRIDGE_DLL (
  echo osqp_fortran_bridge.dll was not produced.
  exit /b 1
)
copy /y "%BRIDGE_DLL%" "%OUT_DIR%\osqp_fortran_bridge.dll" >nul

for %%D in (libwinpthread-1.dll libgcc_s_seh-1.dll) do (
  for /f "delims=" %%P in ('where %%D 2^>nul') do if exist "%%P" copy /y "%%P" "%OUT_DIR%\%%D" >nul
)

echo Backend built: %OUT_DIR%\osqp_fortran_bridge.dll
cd /d "%ROOT%"
where fpm >nul 2>nul && (
  echo Verifying a real QP solve through the Fortran interface...
  fpm run --example basic_qp
  if errorlevel 1 exit /b 1
)
echo Run from the package root with: fpm run
exit /b 0
