@echo off
setlocal enabledelayedexpansion
set ROOT=%~dp0
set FC=gfortran
set FLAGS=-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace -O0 -g -fcheck=all
set LIBS=-llapack -lblas
if not exist "%ROOT%build\windows" mkdir "%ROOT%build\windows"
cd /d "%ROOT%build\windows"
for %%F in (gogarch_kinds gogarch_linalg gogarch_optimizer gogarch_rng gogarch_distributions gogarch_types gogarch_univariate gogarch_orthogonal gogarch_core gogarch_ica gogarch_model gogarch_estimators gogarch) do (
  %FC% %FLAGS% -J. -I. -c "%ROOT%src\%%F.f90" || exit /b 1
)
%FC% %FLAGS% -J. -I. -c "%ROOT%test\test_helpers.f90" || exit /b 1
set OBJS=gogarch.o gogarch_estimators.o gogarch_model.o gogarch_ica.o gogarch_core.o gogarch_orthogonal.o gogarch_univariate.o gogarch_types.o gogarch_distributions.o gogarch_rng.o gogarch_optimizer.o gogarch_linalg.o gogarch_kinds.o
for %%T in (test_core test_univariate_extended test_estimators test_gogarch_extensions) do (
  %FC% %FLAGS% -J. -I. "%ROOT%test\%%T.f90" test_helpers.o %OBJS% %LIBS% -o %%T.exe || exit /b 1
  %%T.exe || exit /b 1
)
%FC% %FLAGS% -J. -I. "%ROOT%app\demo_gogarch.f90" %OBJS% %LIBS% -o demo_gogarch.exe || exit /b 1
%FC% %FLAGS% -J. -I. "%ROOT%example\fit_csv.f90" %OBJS% %LIBS% -o fit_csv.exe || exit /b 1
demo_gogarch.exe || exit /b 1
fit_csv.exe "%ROOT%data\sample_returns.csv" ica || exit /b 1
fit_csv.exe "%ROOT%data\sample_returns.csv" ica garch std 2 0 1 || exit /b 1
fit_csv.exe "%ROOT%data\sample_returns.csv" ica aparch sstd 1 1 1 || exit /b 1
echo Windows runtime-checked build and execution checks passed.
