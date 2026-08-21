@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0..
set BUILD=%ROOT%\build-validate
if exist "%BUILD%" rmdir /s /q "%BUILD%"
mkdir "%BUILD%\mod" "%BUILD%\obj" "%BUILD%\bin"
set FLAGS=-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all -ffpe-trap=invalid,zero,overflow -J"%BUILD%\mod" -I"%BUILD%\mod"
set SRC=compositions_kinds.f90 bayesm_kinds.f90 bayesm_math.f90 bayesm_linalg.f90 bayesm_rng.f90 robustbase_kinds.f90 robustbase_sort.f90 robustbase_scale.f90 robustbase_probability.f90 robustbase_linalg.f90 robustbase_covariance.f90 robustbase_detmcd.f90 tensora_kinds.f90 tensora_types.f90 tensora_index.f90 tensora_core.f90 tensora_linalg.f90 tensora_stats.f90 tensora.f90 compositions_linalg.f90 compositions_geometry.f90 compositions_distributions.f90 compositions_zero.f90 compositions_imputation.f90 compositions_imputation_cache.f90 compositions_stats.f90 compositions_geostat.f90 compositions_counts.f90 compositions_gof.f90 compositions_energy_gof.f90 compositions_outliers.f90 compositions_tensor.f90 compositions.f90
set OBJS=
for %%S in (%SRC%) do (
  gfortran %FLAGS% -c "%ROOT%\src\%%S" -o "%BUILD%\obj\%%~nS.o" || exit /b 1
  set OBJS=!OBJS! "%BUILD%\obj\%%~nS.o"
)
for %%T in ("%ROOT%\test\*.f90") do (
  gfortran %FLAGS% "%%T" !OBJS! -llapack -lblas -o "%BUILD%\bin\%%~nT.exe" || exit /b 1
  "%BUILD%\bin\%%~nT.exe" || exit /b 1
)
gfortran %FLAGS% "%ROOT%\example\demo_compositions.f90" !OBJS! -llapack -lblas -o "%BUILD%\bin\demo_compositions.exe" || exit /b 1
"%BUILD%\bin\demo_compositions.exe" || exit /b 1
endlocal
