@echo off
setlocal
if "%FC%"=="" set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran\mod build_gfortran\obj build_gfortran\bin
%FC% %FLAGS% -c src\cec2013_kinds.f90 -J build_gfortran\mod -o build_gfortran\obj\cec2013_kinds.o || exit /b 1
%FC% %FLAGS% -c src\cec2013_data.f90 -I build_gfortran\mod -J build_gfortran\mod -o build_gfortran\obj\cec2013_data.o || exit /b 1
%FC% %FLAGS% -c src\cec2013_functions.f90 -I build_gfortran\mod -J build_gfortran\mod -o build_gfortran\obj\cec2013_functions.o || exit /b 1
%FC% %FLAGS% -c src\cec2013.f90 -I build_gfortran\mod -J build_gfortran\mod -o build_gfortran\obj\cec2013.o || exit /b 1
set OBJS=build_gfortran\obj\cec2013_kinds.o build_gfortran\obj\cec2013_data.o build_gfortran\obj\cec2013_functions.o build_gfortran\obj\cec2013.o
for %%F in (test\*.f90) do (
  %FC% %FLAGS% "%%F" %OBJS% -I build_gfortran\mod -J build_gfortran\mod -o "build_gfortran\bin\%%~nF.exe" || exit /b 1
  "build_gfortran\bin\%%~nF.exe" || exit /b 1
)
for %%F in (example\*.f90) do (
  %FC% %FLAGS% "%%F" %OBJS% -I build_gfortran\mod -J build_gfortran\mod -o "build_gfortran\bin\%%~nF.exe" || exit /b 1
  "build_gfortran\bin\%%~nF.exe" >nul || exit /b 1
)
echo All strict tests and examples passed.
endlocal
