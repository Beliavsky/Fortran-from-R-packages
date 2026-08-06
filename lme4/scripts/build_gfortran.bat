@echo off
setlocal enabledelayedexpansion
cd /d "%~dp0\.."
if exist build_gfortran rmdir /s /q build_gfortran
mkdir build_gfortran\mod
mkdir build_gfortran\obj
set FLAGS=-std=f2018 -O2 -Wall -Wextra -Wpedantic
gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod -c dependencies\minqa\src\minqa_module.f90 -o build_gfortran\obj\minqa_module.o || exit /b 1
for %%F in (lme4_kinds lme4_types lme4_linalg lme4_covariance lme4_quadrature lme4_family lme4_lmm lme4_lmm_pls lme4_glmm lme4_custom_glmm lme4_aghq lme4_aghq_nd lme4_nlmm lme4_simulation lme4_diagnostics lme4_inference lme4_grouped lme4) do (
  gfortran %FLAGS% -J build_gfortran\mod -I build_gfortran\mod -c src\%%F.f90 -o build_gfortran\obj\%%F.o || exit /b 1
)
gfortran %FLAGS% -I build_gfortran\mod test\test_lme4.f90 build_gfortran\obj\*.o -o build_gfortran\test_lme4.exe || exit /b 1
build_gfortran\test_lme4.exe || exit /b 1
gfortran %FLAGS% -I build_gfortran\mod example\lme4_example.f90 build_gfortran\obj\*.o -o build_gfortran\lme4_example.exe || exit /b 1
build_gfortran\lme4_example.exe || exit /b 1

gfortran %FLAGS% -I build_gfortran\mod example\glmm_extensions_example.f90 build_gfortran\obj\*.o -o build_gfortran\glmm_extensions_example.exe || exit /b 1
build_gfortran\glmm_extensions_example.exe || exit /b 1

gfortran %FLAGS% -I build_gfortran\mod example\advanced_algorithms_example.f90 build_gfortran\obj\*.o -o build_gfortran\advanced_algorithms_example.exe || exit /b 1
build_gfortran\advanced_algorithms_example.exe || exit /b 1
