@echo off
setlocal
cd /d %~dp0\..
if "%FC%"=="" set FC=gfortran
set BUILD=build\checked
if exist %BUILD% rmdir /s /q %BUILD%
mkdir %BUILD%\mod %BUILD%\obj %BUILD%\bin
set FLAGS=-std=f2018 -Wall -Wextra -fcheck=all -fbacktrace -O0 -g
for %%S in (corpcor_kinds corpcor_types corpcor_linalg corpcor_weighted corpcor_matrix_tools corpcor_shrinkage corpcor) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c vendor\corpcor\src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
)
for %%S in (ren_kinds ren_types ren_linalg ren_random ren_regularization ren_portfolio ren_analysis ren) do (
  %FC% %FLAGS% -J %BUILD%\mod -I %BUILD%\mod -c src\%%S.f90 -o %BUILD%\obj\%%S.o || exit /b 1
)
ar rcs %BUILD%\libren.a %BUILD%\obj\*.o || exit /b 1
for %%T in (test_analysis test_portfolios test_prepare_data) do (
  %FC% %FLAGS% -I %BUILD%\mod test\%%T.f90 %BUILD%\libren.a -o %BUILD%\bin\%%T.exe || exit /b 1
  %BUILD%\bin\%%T.exe || exit /b 1
)
%FC% %FLAGS% -I %BUILD%\mod example\basic_portfolios.f90 %BUILD%\libren.a -o %BUILD%\bin\basic_portfolios.exe || exit /b 1
%FC% %FLAGS% -I %BUILD%\mod app\demo_ren.f90 %BUILD%\libren.a -o %BUILD%\bin\demo_ren.exe || exit /b 1
echo Checked GNU Fortran build: PASS
endlocal
