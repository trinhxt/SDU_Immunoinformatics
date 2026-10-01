@echo off
cd /d "%~dp0"
title Human Disease CDR3 Antibody Peptides Database

:: 1. Check if Rscript is in PATH
where Rscript >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    set "RSCRIPT_BIN=Rscript"
    goto :CHECK_DEPS
)

:: 2. Search common R installation paths if not in system PATH
for /d %%D in ("C:\Program Files\R\R-*") do (
    if exist "%%D\bin\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\Rscript.exe"
        goto :CHECK_DEPS
    )
    if exist "%%D\bin\x64\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\x64\Rscript.exe"
        goto :CHECK_DEPS
    )
)

for /d %%D in ("C:\PROGRA~1\R\R-*") do (
    if exist "%%D\bin\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\Rscript.exe"
        goto :CHECK_DEPS
    )
    if exist "%%D\bin\x64\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\x64\Rscript.exe"
        goto :CHECK_DEPS
    )
)

for /d %%D in ("%LocalAppData%\Programs\R\R-*") do (
    if exist "%%D\bin\Rscript.exe" (
        set "RSCRIPT_BIN=%%D\bin\Rscript.exe"
        goto :CHECK_DEPS
    )
)

echo Error: Rscript.exe could not be located in PATH or standard installation folders.
echo Please ensure R (>= 4.2) is installed from https://cloud.r-project.org/
echo and optionally added to your system PATH.
echo.
pause
exit /b 1

:CHECK_DEPS
echo Initializing OASpepDB environment...
"%RSCRIPT_BIN%" -e "if (!requireNamespace('shiny', quietly = TRUE)) { message('Installing Shiny framework...'); install.packages('shiny', repos = 'https://cloud.r-project.org') }"

:RUN_APP
echo Starting Shiny server and opening web browser...
echo (Closing the browser or this console window will stop the application)
echo.

"%RSCRIPT_BIN%" -e "tryCatch(shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080), error = function(e) { message('Notice: Port 8080 is unavailable, selecting dynamic port...'); shiny::runApp('DBquery.R', launch.browser = TRUE) })"

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo Application exited with status code %ERRORLEVEL%.
    pause
)
