################################################################################
## scripts/01_download.R
##
## Part 1: Download OAS CSV.GZ files to RAW folder
##
## Description:
##   Parses 'bulk_download.sh' (wget commands) to extract URLs and downloads 
##   them directly to the 'raw' directory specified in 00_config.R.
##
## Outputs:
##   - Raw *.csv.gz files in P$raw_dir
##
## Notes:
##   - Resume capability: Skips files that already exist in P$raw_dir.
##   - Handles the specific format: "wget https://..."
################################################################################

# ==============================================================================
# Load config
# ==============================================================================
source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
})

# ==============================================================================
# Helper Functions
# ==============================================================================

# Parse the OAS shell script to extract URLs
read_bulk_urls <- function(sh_file) {
  if (!file.exists(sh_file)) return(character())
  lines <- readLines(sh_file, warn = FALSE)
  
  # Filter for lines starting with "wget" and containing "http"
  # This matches the format: wget https://...
  valid_lines <- grep("^wget.*http", lines, value = TRUE)
  
  if (length(valid_lines) == 0) return(character())
  
  # Extract the URL part (everything after 'wget ')
  # We assume the URL is the last distinct "word" or simply the second token
  # Simple regex: find 'http' and take everything until the end of the line (or space)
  urls <- sub(".*(https?://[^ ]+).*", "\\1", valid_lines)
  
  return(unique(urls))
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
  
  # 1. Setup paths from Config (P global)
  sh_file <- P$shell_script
  raw_dir <- P$raw_dir
  
  if (!file.exists(sh_file)) stop("Bulk download script not found: ", sh_file)
  if (!dir.exists(raw_dir)) dir.create(raw_dir, recursive = TRUE)
  
  # 2. Get Target List
  cat("Reading URLs from:", sh_file, "\n")
  urls <- read_bulk_urls(sh_file)
  
  if (length(urls) == 0) stop("No URLs extracted from ", sh_file)
  
  target_files <- basename(urls)
  n_total <- length(urls)
  
  cat("Found", n_total, "files to download.\n")
  cat("Target Directory:", raw_dir, "\n\n")
  
  # 3. Download Loop
  success_count <- 0
  skip_count    <- 0
  fail_count    <- 0
  
  for (i in seq_len(n_total)) {
    url <- urls[i]
    fn  <- target_files[i]
    dest_path <- file.path(raw_dir, fn)
    
    # Progress indicator (every 10 files or 10%)
    if (i %% 10 == 0 || i == n_total) {
      cat(sprintf("Progress: %d / %d (%.1f%%)\r", i, n_total, 100 * i / n_total))
      flush.console()
    }
    
    # Check Resume (Skip if exists and size > 0)
    if (file.exists(dest_path) && file.size(dest_path) > 0) {
      skip_count <- skip_count + 1
      next
    }
    
    # Download
    # mode="wb" is critical for binary (gzip) files on Windows
    tryCatch({
      download.file(url, dest_path, mode = "wb", quiet = TRUE)
      success_count <- success_count + 1
    }, error = function(e) {
      cat(sprintf("\nFAILED: %s - %s\n", fn, conditionMessage(e)))
      fail_count <<- fail_count + 1
    })
  }
  
  # 4. Report
  cat("\n\n==========================================\n")
  cat(" Download Summary\n")
  cat("==========================================\n")
  cat(" Total URLs:  ", n_total, "\n")
  cat(" Downloaded:  ", success_count, "\n")
  cat(" Skipped:     ", skip_count, "\n")
  cat(" Failed:      ", fail_count, "\n")
  cat("==========================================\n")
}

# Run
main()
