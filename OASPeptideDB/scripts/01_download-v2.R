# ==============================================================================
# R Script to Download OAS Data with Resume Capability
# ==============================================================================

# 1. Configuration
# ------------------------------------------------------------------------------
input_sh_file <- "D:/OAS/unpaired/bulk_download_human_unpaired_heavychain_disease.sh"
output_dir    <- "D:/OAS/unpaired/raw"  # Change this to where you want files saved

# Increase timeout limit for large files (default is often too short)
options(timeout = 600) 

# 2. Setup
# ------------------------------------------------------------------------------
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

if (!file.exists(input_sh_file)) {
  stop("Input file not found: ", input_sh_file)
}

# 3. Parse the Bash Script
# ------------------------------------------------------------------------------
cat("Reading download list from:", input_sh_file, "\n")
lines <- readLines(input_sh_file, warn = FALSE)

# Extract URLs: Find lines with http, clean off 'wget' and other flags
# This regex looks for 'http' or 'https' and grabs the string until a space or line end
raw_urls <- grep("http", lines, value = TRUE)
urls <- sub(".*(https?://[^ ]+).*", "\\1", raw_urls)
urls <- unique(urls) # Ensure no duplicates

total_files <- length(urls)
cat("Found", total_files, "unique files to download.\n\n")

# 4. Download Loop
# ------------------------------------------------------------------------------
success <- 0
skipped <- 0
failed  <- 0

for (i in seq_along(urls)) {
  url <- urls[i]
  file_name <- basename(url)
  dest_path <- file.path(output_dir, file_name)
  
  # --- RESUME LOGIC ---
  # Check if file exists AND has content (size > 0 bytes)
  if (file.exists(dest_path) && file.size(dest_path) > 0) {
    # Optional: You could add a progress print for skips, but it might clutter the console
    # cat(sprintf("[%d/%d] Skipping (Already exists): %s\n", i, total_files, file_name))
    skipped <- skipped + 1
    next
  }
  
  # --- DOWNLOAD LOGIC ---
  cat(sprintf("[%d/%d] Downloading: %s\n", i, total_files, file_name))
  
  tryCatch({
    # mode="wb" is CRITICAL for binary files (.gz) on Windows/WSL
    download.file(url, dest_path, mode = "wb", quiet = TRUE)
    
    # Validation: Check if the resulting file is actually there
    if (file.exists(dest_path) && file.size(dest_path) > 0) {
      success <- success + 1
    } else {
      cat("   -> Warning: Download finished but file is empty/missing.\n")
      failed <- failed + 1
      if(file.exists(dest_path)) unlink(dest_path) # Delete empty file
    }
    
  }, error = function(e) {
    cat(sprintf("   -> FAILED: %s\n", e$message))
    failed <<- failed + 1
    # Clean up partial download if it exists
    if(file.exists(dest_path)) unlink(dest_path)
  })
}

# 5. Summary
# ------------------------------------------------------------------------------
cat("\n===============================================\n")
cat("Download Complete\n")
cat("===============================================\n")
cat("Total Files: ", total_files, "\n")
cat("Skipped:     ", skipped, " (Already existed)\n")
cat("Downloaded:  ", success, "\n")
cat("Failed:      ", failed, "\n")
cat("===============================================\n")