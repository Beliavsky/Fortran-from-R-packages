@echo off
setlocal
cd /d "%~dp0.."
if not defined FC set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
if exist build_strict rmdir /s /q build_strict
mkdir build_strict
cd build_strict
%FC% %FLAGS% -J . -I . -c ..\src\soma_kinds.f90 ..\src\soma_types.f90 ..\src\soma_random.f90 ..\src\soma_optimizer.f90 ..\src\soma.f90
if errorlevel 1 exit /b 1
set OBJS=soma_kinds.o soma_types.o soma_random.o soma_optimizer.o soma.o
for %%F in (..\test\*.f90) do (
  %FC% %FLAGS% -J . -I . %OBJS% "%%F" -o "%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%%~nF.exe"
  if errorlevel 1 exit /b 1
)
for %%F in (..\example\*.f90) do (
  %FC% %FLAGS% -J . -I . %OBJS% "%%F" -o "%%~nF.exe"
  if errorlevel 1 exit /b 1
  "%%~nF.exe"
  if errorlevel 1 exit /b 1
)
endlocal
