@echo off
fpm run
if errorlevel 1 exit /b 1

fpm run --example fit_csv -- data\sample_returns.csv
if errorlevel 1 exit /b 1


fpm run --example complete_workflows
if errorlevel 1 exit /b 1
