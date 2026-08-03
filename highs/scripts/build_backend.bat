@echo off
setlocal EnableExtensions EnableDelayedExpansion
for %%I in ("%~dp0..") do set "ROOT=%%~fI"
set "HIGHS_SRC=%ROOT%\upstream\highs-r\inst\HiGHS"
set "HIGHS_BUILD=%ROOT%\backend\highs-build"
set "HIGHS_INSTALL=%ROOT%\backend\highs-install"
set "BRIDGE_BUILD=%ROOT%\backend\bridge-build"
set "OUT_DIR=%ROOT%\backend\bin"

where cmake >nul 2>nul || (
  echo cmake is required to build the bundled HiGHS backend.
  exit /b 1
)
where gcc >nul 2>nul || (
  echo gcc was not found. Use the same MinGW toolchain as gfortran.
  exit /b 1
)
where g++ >nul 2>nul || (
  echo g++ was not found. Install the C++ component of the MinGW toolchain.
  exit /b 1
)

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"
set "CPLUS_INCLUDE_PATH=%ROOT%\bridge\compat;%CPLUS_INCLUDE_PATH%"
set "GENERATOR="
where ninja >nul 2>nul && set "GENERATOR=-G Ninja"
if not defined GENERATOR (
  where mingw32-make >nul 2>nul && set GENERATOR=-G "MinGW Makefiles"
)
if not defined GENERATOR (
  echo Ninja or mingw32-make is required for the MinGW backend build.
  exit /b 1
)

cmake -S "%HIGHS_SRC%" -B "%HIGHS_BUILD%" %GENERATOR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc ^
  -DCMAKE_CXX_COMPILER=g++ ^
  -DCMAKE_INSTALL_PREFIX="%HIGHS_INSTALL%" ^
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON ^
  -DBUILD_SHARED_LIBS=OFF ^
  -DBUILD_TESTING=OFF ^
  -DBUILD_EXAMPLES=OFF ^
  -DZLIB=OFF ^
  -DFAST_BUILD=ON ^
  -DCUPDLP_GPU=OFF
if errorlevel 1 exit /b 1
cmake --build "%HIGHS_BUILD%" --target install --parallel
if errorlevel 1 exit /b 1

cmake -S "%ROOT%\bridge" -B "%BRIDGE_BUILD%" %GENERATOR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DCMAKE_C_COMPILER=gcc ^
  -DCMAKE_CXX_COMPILER=g++ ^
  -DHIGHS_ROOT="%HIGHS_INSTALL%" ^
  -DCMAKE_RUNTIME_OUTPUT_DIRECTORY="%OUT_DIR%"
if errorlevel 1 exit /b 1
cmake --build "%BRIDGE_BUILD%" --config Release --parallel
if errorlevel 1 exit /b 1

set "BRIDGE_DLL="
for /r "%BRIDGE_BUILD%" %%F in (highs_fortran_bridge.dll) do set "BRIDGE_DLL=%%~fF"
if not defined BRIDGE_DLL if exist "%OUT_DIR%\highs_fortran_bridge.dll" set "BRIDGE_DLL=%OUT_DIR%\highs_fortran_bridge.dll"
if not defined BRIDGE_DLL (
  echo highs_fortran_bridge.dll was not produced.
  exit /b 1
)
copy /y "%BRIDGE_DLL%" "%OUT_DIR%\highs_fortran_bridge.dll" >nul

for %%D in (libwinpthread-1.dll libgcc_s_seh-1.dll libstdc++-6.dll) do (
  for /f "delims=" %%P in ('where %%D 2^>nul') do if exist "%%P" copy /y "%%P" "%OUT_DIR%\%%D" >nul
)

echo Backend built: %OUT_DIR%\highs_fortran_bridge.dll
cd /d "%ROOT%"
where fpm >nul 2>nul && (
  echo Verifying a real LP solve through the Fortran interface...
  fpm run --example lp_example
  if errorlevel 1 exit /b 1
)
echo Run from the package root with: fpm run
exit /b 0
