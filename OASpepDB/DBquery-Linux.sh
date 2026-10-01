#!/usr/bin/env bash
# ==============================================================================
# OASpepDB: Human Disease CDR3 Antibody Peptides Database - Linux Launcher
# SDU Immunoinformatics
# ==============================================================================
set -e

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

if ! command -v Rscript &> /dev/null; then
    echo "Error: 'Rscript' was not found on your system PATH."
    echo "Please install R (e.g., sudo apt install r-base)."
    echo ""
    read -n 1 -s -r -p "Press any key to exit..."
    exit 1
fi

echo "Checking required R packages..."
Rscript -e "req <- c('shiny','bslib','DT','duckdb','DBI','arrow','ggplot2','plotly','dplyr','htmlwidgets','zip'); missing <- req[!req %in% installed.packages()[,'Package']]; if (length(missing) > 0) { message('Installing missing R packages: ', paste(missing, collapse=', ')); install.packages(missing, repos='https://cloud.r-project.org') }"

echo "Starting Shiny server and opening web browser..."
echo "(Closing the browser or terminal will stop the application)"
echo ""

Rscript -e "tryCatch(shiny::runApp('DBquery.R', launch.browser = TRUE, port = 8080), error = function(e) { message('Notice: Port 8080 is unavailable, falling back to dynamic port...'); shiny::runApp('DBquery.R', launch.browser = TRUE) })"
