@echo off
rem Runs the test cases with origin/main and with the current working tree,
rem and writes the comparison to docs\PR_results\<date>_<name>\
rem (see tools\compare_branches.py). Double-click to run.
cd /d "%~dp0.."

rem Use Anaconda/Miniconda if installed, otherwise python from the PATH
set "PY=python"
if exist "%LOCALAPPDATA%\anaconda3\python.exe" set "PY=%LOCALAPPDATA%\anaconda3\python.exe"
if exist "%USERPROFILE%\anaconda3\python.exe" set "PY=%USERPROFILE%\anaconda3\python.exe"
if exist "%USERPROFILE%\miniconda3\python.exe" set "PY=%USERPROFILE%\miniconda3\python.exe"

set /p "NAME=Short name of the PR (e.g. melting-fix): "
"%PY%" tools\compare_branches.py --name "%NAME%"
pause
