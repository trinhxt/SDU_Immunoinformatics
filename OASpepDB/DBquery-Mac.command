#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Human Disease CDR3 Antibody Peptides Database - macOS Launcher
# SDU Immunoinformatics
# ==============================================================================
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "Starting Shiny server and opening web browser..."
echo "(Closing the browser window will exit and close this window)"
echo ""

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

Rscript -e "shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080)"
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    osascript -e 'tell application "Terminal" to close (every window whose name contains "DBquery")' &> /dev/null || true
fi
