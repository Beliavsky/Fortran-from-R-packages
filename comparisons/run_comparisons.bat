@echo off
setlocal
cd /d "%~dp0.."
python comparisons\run_comparisons.py
exit /b %errorlevel%
