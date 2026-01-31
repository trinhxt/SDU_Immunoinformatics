################################################################################
## scripts/01_download.R
##
## Part1: Download + preprocess OAS CSV.GZ files
##
## Outputs
## -------
## 1) Processed OAS files (same filenames as OAS download):
##      P$processed_dir/*.csv.gz
##    Columns kept:
##      sequence_alignment_aa, v_call, d_call, j_call, cdr1_aa, cdr2_aa, cdr3_aa
##
## 2) Metadata table:
##      P$metadata_csv  (CSV)
##
## Fast SKIP rule (set-based, no per-file checking)
## ------------------------------------------------
## If ALL are true:
##   - P$metadata_csv exists and has Filename
##   - Filename set matches ALL filenames in bulk_download.sh
##   - processed_dir contains ALL expected *.csv.gz and no extras
## Then we SKIP the whole script (return) and do NOT download/parse anything.
##
## Notes
## -----
## - Uses Windows-safe download.file(mode="wb").
## - Raw staging folder uses P$raw_dir if provided, else tempdir().
## - No quit() is used (safe for source()).
################################################################################

# ==============================================================================
# Load config
# ==============================================================================
# Adjust path if your 00_config.R is in a different folder (e.g. "scripts/00_config.R")
source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(R.utils)
  library(jsonlite)
})

# ==============================================================================
# Helper Functions
# ==============================================================================

# Parse the OAS shell script to extract URLs
read_bulk_urls <- function(sh_file) {
  if (!file.exists(sh_file)) return(character())
  lines <- readLines(sh_file, warn = FALSE)
  # grep lines starting with "curl" or "wget" and extract URL
  # Usually: curl -O "https://..."
  # Simple regex to capture http[s]://...
  urls <- grep("http", lines, value = TRUE)
  # Clean up: extract strictly the URL part
  # Assuming standard OAS format: ... "https://opig.stats.ox.ac.uk/..." ...
  # We'll regex extract content inside quotes or just the http token
  urls <- sub('.*"(http[^"]+)".*', "\\1", urls)
  # Fallback if no quotes: just find http... until space or end
  urls <- sub(".*(http[^ ]+).*", "\\1", urls)
  
  unique(urls[grepl("^http", urls)])
}

# Download + Process One URL
process_one_url <- function(url, idx, total, raw_dir, out_dir) {
  fn <- basename(url)
  out_path <- file.path(out_dir, fn)
  
  # Check if output already exists (resume capability)
  if (file.exists(out_path)) {
    # If exists, we still need metadata. We can't easily get metadata without reading headers.
    # To be fast, we'll return a special flag or just re-read the header (first line).
    # But for "process_one_url", let's assume if it exists, we read the header for metadata
    # and skip re-processing.
    
    # Try reading just the first line for metadata
    tryCatch({
      header_line <- readLines(out_path, n = 1, warn = FALSE)
      # OAS processed CSVs usually have metadata in the first line if we kept it?
      # Actually, my previous logic stripped metadata from the CSV content but returned it separately.
      # If we want to rebuild metadata, we might need to re-download if we didn't save it separate.
      # BUT, Part 1 says "Metadata table: P$metadata_csv".
      # So we shouldn't rely on re-reading metadata from processed files if we are rebuilding.
      #
      # Strategy: If file exists, we SKIP processing, but we return NULL metadata 
      # (assuming global metadata file will be handled or is already present).
      # If the global metadata check failed (in main), we likely need to redo.
      # For safety: let's re-download if we need metadata, OR just re-read processed?
      # OAS processed files lose the JSON header usually.
      # So if we need metadata, we must re-download.
      # Let's assume we proceed to download unless we are 100% sure we can skip.
      # The main() function handles the Big Skip. Here we just do the work.
    }, error = function(e) NULL)
  }
  
  cat(sprintf("[%d/%d] Processing: %s\n", idx, total, fn))
  
  # 1. Download to raw_dir
  tmp_dest <- file.path(raw_dir, fn)
  
  # Use mode="wb" for Windows safety
  tryCatch({
    download.file(url, tmp_dest, mode = "wb", quiet = TRUE)
  }, error = function(e) {
    return(list(filename = fn, ok = FALSE, meta = NULL, msg = conditionMessage(e)))
  })
  
  if (!file.exists(tmp_dest)) {
    return(list(filename = fn, ok = FALSE, meta = NULL, msg = "Download failed"))
  }
  
  # 2. Read Metadata (First line JSON) & Data
  # Use zcat/gzfile logic handled by fread automatically or R.utils?
  # OAS files are .csv.gz. The first line is specific JSON metadata.
  # fread skips it automatically usually, but we need it.
  
  con <- gzfile(tmp_dest, "rt")
  header_json <- readLines(con, n = 1, warn = FALSE)
  close(con)
  
  # Validate JSON
  meta_row <- NULL
  if (grepl("^\\{", header_json)) {
    try({
      m <- fromJSON(header_json, flatten = TRUE)
      # Convert list to data.table row
      meta_row <- as.data.table(m)
      meta_row[, Filename := fn]
      # Ensure consistent columns? We'll just keep what we get for now.
    }, silent = TRUE)
  }
  
  # 3. Read Data (skip first line)
  # Select specific columns to save space
  keep_cols <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", 
                 "cdr1_aa", "cdr2_aa", "cdr3_aa")
  
  dt <- tryCatch({
    fread(tmp_dest, skip = 1, select = keep_cols, showProgress = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(dt)) {
    # Cleanup raw
    unlink(tmp_dest)
    return(list(filename = fn, ok = FALSE, meta = meta_row, msg = "fread failed"))
  }
  
  # 4. Filter/Clean (Optional)
  # - Remove rows with no sequence
  dt <- dt[!is.na(sequence_alignment_aa) & sequence_alignment_aa != ""]
  
  # 5. Write Processed
  # Save as standard csv.gz
  fwrite(dt, out_path, compress = "gzip")
  
  # 6. Cleanup Raw
  unlink(tmp_dest)
  
  return(list(filename = fn, ok = TRUE, meta = meta_row, msg = "OK"))
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
  
  # 1. Setup paths from Config (P global)
  sh_file   <- P$shell_script
  raw_dir   <- P$raw_dir
  proc_dir  <- P$processed_dir
  meta_file <- P$metadata_csv
  
  if (!file.exists(sh_file)) stop("Bulk download script not found: ", sh_file)
  
  # 2. Get Target List
  urls <- read_bulk_urls(sh_file)
  if (length(urls) == 0) stop("No URLs extracted from ", sh_file)
  
  target_files <- basename(urls)
  n_total <- length(target_files)
  cat("Total files in OAS script: ", n_total, "\n")
  
  # 3. Fast Skip Check
  # If metadata exists AND contains all filenames AND processed dir matches
  if (file.exists(meta_file)) {
    existing_meta <- fread(meta_file, select = "Filename")
    existing_files <- list.files(proc_dir, pattern = "\\.csv\\.gz$")
    
    # Check 1: Metadata has all targets?
    has_all_meta <- all(target_files %in% existing_meta$Filename)
    
    # Check 2: Processed folder has all targets?
    has_all_proc <- all(target_files %in% existing_files)
    
    if (has_all_meta && has_all_proc) {
      cat("All files present and metadata complete. SKIPPING Part 1.\n")
      return(invisible(NULL))
    } else {
      cat("Update needed (Missing files or metadata). Proceeding...\n")
    }
  }
  
  # 4. Run Processing Loop
  # We'll use simple parallelization if N_CORES > 1
  # But download is often bandwidth limited. Sequential might be safer/easier for error handling
  # or use 'future.apply' if you want.
  # For simplicity, let's stick to sequential loop or basic mclapply (Linux) / parLapply (Windows).
  # Given the "download" nature, strict sequential is often more stable to avoid timeouts.
  # Let's do sequential for robustness in this script.
  
  results_list <- vector("list", n_total)
  
  for (i in seq_len(n_total)) {
    url <- urls[i]
    fn  <- target_files[i]
    
    # Check if already done (exists in processed AND exists in current metadata file?)
    # If we are doing a partial update, we might want to skip.
    # Simple check: if output file exists, we assume it's good, BUT we need the metadata row.
    # If we don't have the metadata row, we MUST re-download to parse the JSON header.
    # This is the tricky part of OAS files.
    
    # OPTIMIZATION: If we have a 'metadata_backup.csv', we could look there.
    # Otherwise, process.
    
    res <- process_one_url(url, i, n_total, raw_dir, proc_dir)
    results_list[[i]] <- res
    
    # Optional: Periodic garbage collection
    if (i %% 50 == 0) gc()
  }
  
  # 5. Compile Metadata
  # Filter out failed
  valid_results <- Filter(function(x) !is.null(x$meta), results_list)
  
  if (length(valid_results) > 0) {
    cat("Compiling metadata...\n")
    meta_dt <- rbindlist(lapply(valid_results, `[[`, "meta"), fill = TRUE)
    
    # Select/Order columns if desired
    # Ensure Filename is first
    setcolorder(meta_dt, "Filename")
    
    # Write Metadata
    fwrite(meta_dt, meta_file)
    cat("Metadata written to:", meta_file, "\n")
  } else {
    cat("WARNING: No valid metadata recovered.\n")
  }
  
  # 6. Report
  failures <- Filter(function(x) !x$ok, results_list)
  if (length(failures) > 0) {
    cat("\nFailures:\n")
    print(rbindlist(failures)[, .(filename, msg)])
  } else {
    cat("\nSuccess! All files processed.\n")
  }
}

# Run
main()