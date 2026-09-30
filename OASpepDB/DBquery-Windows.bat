@echo off
cd /d "%~dp0"
title Human Disease CDR3 Antibody Peptides Database

:: 1. Check if Rscript is in PATH
where Rscript >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    set "RSCRIPT_BIN=Rscript"
    goto :RUN_APP
)

:: 2. Search common R installation paths if not in system PATH
for /d %%D in ("C:\Program Files\R\R-*") do (
    if exist "%%D\bin\x64\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\x64\Rscript.exe"
        goto :RUN_APP
    )
)

for /d %%D in ("C:\PROGRA~1\R\R-*") do (
    if exist "%%D\bin\x64\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\x64\Rscript.exe"
        goto :RUN_APP
    )
)

echo Error: Rscript.exe could not be located in PATH or Program Files.
echo Please ensure R is installed and added to your system PATH.
echo.
pause
exit /b 1

:RUN_APP
echo Starting Shiny server and opening web browser...
echo (Closing the browser window will exit and close this window)
echo.

"%RSCRIPT_BIN%" -e "shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080)"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Application exited with status code %ERRORLEVEL%.
    pause
)
