@echo off
setlocal

pushd "%~dp0" >nul
if errorlevel 1 (
    echo ERROR: Could not enter the mixtools package directory.
    exit /b 1
)

where Rscript >nul 2>&1
if errorlevel 1 (
    echo ERROR: Rscript was not found on PATH.
    goto :failure
)

where fpm >nul 2>&1
if errorlevel 1 (
    echo ERROR: fpm was not found on PATH.
    goto :failure
)

if not exist "example\mvnormal_mixture_data.txt" (
    echo Generating multivariate-normal mixture data...
    Rscript "example\generate_mvnormal_mixture.R"
    if errorlevel 1 goto :failure
) else (
    echo Using existing data file "example\mvnormal_mixture_data.txt".
)

echo.
echo Fitting the mixture with R mixtools...
Rscript "example\fit_mvnormal_mixture.R"
if errorlevel 1 goto :failure

echo.
echo Fitting the mixture with the Fortran translation...
fpm run --example fit_mvnormal_mixture
if errorlevel 1 goto :failure

echo.
echo Comparison completed successfully.
echo Reports:
echo   example\mvnormal_mixture_fit_r.txt
echo   example\mvnormal_mixture_fit_fortran.txt
popd >nul
endlocal
exit /b 0

:failure
set "comparison_exit_code=%errorlevel%"
if "%comparison_exit_code%"=="0" set "comparison_exit_code=1"
echo.
echo ERROR: The comparison did not complete.
popd >nul
endlocal & exit /b %comparison_exit_code%
