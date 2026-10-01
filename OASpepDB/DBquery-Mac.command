#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Human Disease CDR3 Antibody Peptides Database - macOS Launcher
# SDU Immunoinformatics
# ==============================================================================
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

if ! command -v Rscript &> /dev/null; then
    if [ -x "/usr/local/bin/Rscript" ]; then
        export PATH="/usr/local/bin:$PATH"
    elif [ -x "/opt/homebrew/bin/Rscript" ]; then
        export PATH="/opt/homebrew/bin:$PATH"
    elif [ -x "/Library/Frameworks/R.framework/Resources/bin/Rscript" ]; then
        export PATH="/Library/Frameworks/R.framework/Resources/bin:$PATH"
    else
        echo "Error: 'Rscript' was not found on your system."
        echo "Please install R for macOS from https://cloud.r-project.org/"
        echo ""
        read -n 1 -s -r -p "Press any key to exit..."
        exit 1
    fi
fi

echo "Checking required R packages..."
Rscript -e "req <- c('shiny','bslib','DT','duckdb','DBI','arrow','ggplot2','plotly','dplyr','htmlwidgets','zip'); missing <- req[!req %in% installed.packages()[,'Package']]; if (length(missing) > 0) { message('Installing missing R packages: ', paste(missing, collapse=', ')); install.packages(missing, repos='https://cloud.r-project.org') }"

echo "Starting Shiny server and opening web browser..."
echo "(Closing the browser or terminal will stop the application)"
echo ""

Rscript -e "tryCatch(shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080), error = function(e) { message('Notice: Port 8080 is unavailable, falling back to dynamic port...'); shiny::runApp('DBquery.R', launch.browser = TRUE) })"
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    osascript -e 'tell application "Terminal" to close (every window whose name contains "DBquery")' &> /dev/null || true
fi
