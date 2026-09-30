#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Human Disease CDR3 Antibody Peptides Database - Linux Launcher
# SDU Immunoinformatics
# ==============================================================================
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

echo "Starting Shiny server and opening web browser..."
echo "(Closing the browser window will exit and close this window)"
echo ""

if ! command -v Rscript &> /dev/null; then
    echo "Error: 'Rscript' was not found on your system PATH."
    echo "Please install R (e.g., sudo apt install r-base)."
    echo ""
    read -n 1 -s -r -p "Press any key to exit..."
    exit 1
fi

Rscript -e "shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080)"
