@echo off
REM ============================================================================
REM OASpepDB: Reference Proteome Download Wrapper (Windows)
REM Downloads UniProt Swiss-Prot and NCBI RefSeq Human Proteomes into Data/
REM ============================================================================
set SCRIPT_DIR=%~dp0
set PY_SCRIPT=%SCRIPT_DIR%..\Scripts\download_reference_proteomes.py

python "%PY_SCRIPT%" --out-dir "%SCRIPT_DIR%" %*
if errorlevel 1 (
    echo [ERROR] Download failed. Ensure Python is installed and on PATH.
    pause
)
