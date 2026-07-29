@echo off
setlocal EnableExtensions EnableDelayedExpansion
set MODE=%1
if "%MODE%"=="" set MODE=strict
if not "%MODE%"=="strict" if not "%MODE%"=="optimized" (
  echo Usage: run_gfortran_tests.bat [strict^|optimized]
  exit /b 2
)
if "%FC%"=="" set FC=gfortran
set BUILD=build\%MODE%
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -ffree-line-length-none -Wall -Wextra -Werror -Wimplicit-interface
if "%MODE%"=="strict" set FLAGS=%FLAGS% -Wconversion-extra -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow
if "%MODE%"=="optimized" set FLAGS=%FLAGS% -O3
set SOURCES=acdm_kinds acdm_math acdm_distributions acdm_models acdm_fit acdm_data acdm_diagnostics acdm_profiles acdm acdm_api
set OBJECTS=
for %%S in (%SOURCES%) do (
  %FC% %FLAGS% -J "%BUILD%\mod" -I "%BUILD%\mod" -c "src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
  set OBJECTS=!OBJECTS! "%BUILD%\obj\%%S.o"
)
for %%D in (test app example) do (
  for %%F in (%%D\*.f90) do (
    %FC% %FLAGS% -I "%BUILD%\mod" "%%F" !OBJECTS! -o "%BUILD%\bin\%%~nF.exe" || exit /b 1
    echo ==^> %%D\%%~nF
    "%BUILD%\bin\%%~nF.exe" || exit /b 1
  )
)
echo All %MODE% targets passed.
endlocal
