@echo off
setlocal enabledelayedexpansion
if "%FC%"=="" set FC=gfortran
set ROOT=%~dp0
set BUILD=%ROOT%build-direct
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -Wall -Wextra -Wpedantic -Werror -fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace -O0 -g -J"%BUILD%\mod" -I"%BUILD%\mod"
for %%S in (jfe_kinds jfe_stats jfe_performance jfe) do (
  %FC% %FLAGS% -c "%ROOT%src\%%S.f90" -o "%BUILD%\obj\%%S.o" || exit /b 1
)
set OBJECTS="%BUILD%\obj\jfe_kinds.o" "%BUILD%\obj\jfe_stats.o" "%BUILD%\obj\jfe_performance.o" "%BUILD%\obj\jfe.o"
for %%T in ("%ROOT%test\*.f90") do (
  %FC% %FLAGS% "%%T" %OBJECTS% -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
endlocal
