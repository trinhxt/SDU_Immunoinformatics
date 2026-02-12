################################################################################
## scripts/01_download.R
##
## Part 1: Download + Preprocess OAS Data
##
## Goal:
##   1. Extract URLs from shell script.
##   2. Download to raw folder.
##   3. Parse Metadata (Handling NaN, quotes, escaping).
##   4. Save Processed Data to 'antibody' folder.
##   5. Append Metadata to 'OAS_metadata.csv' incrementally.
##   6. FINALIZE: Generate 'Patient' column (Author + Subject) at the end.
##
## RESUME LOGIC:
##   - Checks P$processed_dir. If file exists -> SKIP.
################################################################################

# ------------------------------------------------------------------------------
# 0. Safety Check
# ------------------------------------------------------------------------------
if (!exists("P")) {
  stop("CRITICAL ERROR: Configuration 'P' not found.\n", 
       "Please run this script by executing 'source(\"DBbuild.R\")' instead.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(R.utils)
  library(dplyr)
})

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

read_bulk_urls <- function(sh_file) {
  if (!file.exists(sh_file)) return(character())
  lines <- readLines(sh_file, warn = FALSE)
  url_lines <- grep("http", lines, value = TRUE)
  urls <- sub('.*"(http[^"]+)".*', "\\1", url_lines)
  urls <- sub(".*(http[^ ]+).*", "\\1", urls)
  unique(urls[grepl("^http", urls)])
}

# --- CUSTOM METADATA PARSER ---
parse_first_line_metadata <- function(filepath) {
  # Open connection to GZ file
  con <- gzfile(filepath, "rt")
  on.exit(close(con))
  
  # Read first line
  first_line <- readLines(con, n = 1, warn = FALSE)
  
  if (length(first_line) == 0) return(NULL)
  
  # 1. Preprocess to replace NaN with null
  clean_line <- gsub('\\bNaN\\b', 'null', first_line)
  
  # 2. Remove leading and trailing quotes
  clean_line <- gsub('^"|"$', '', clean_line)
  
  # 3. Replace double double-quotes with single double-quotes
  clean_line <- gsub('""', '"', clean_line)
  
  # 4. Parse
  tryCatch({
    jsonlite::fromJSON(clean_line, flatten = TRUE)
  }, error = function(e) {
    return(NULL)
  })
}

process_one_url <- function(url, idx, total, raw_dir, out_dir) {
  fn       <- basename(url)
  out_path <- file.path(out_dir, fn)
  tmp_dest <- file.path(raw_dir, fn)
  
  cat(sprintf("[%d/%d] Downloading: %s ... ", idx, total, fn))
  
  # 1. Download
  tryCatch({
    download.file(url, tmp_dest, mode = "wb", quiet = TRUE)
  }, error = function(e) {
    cat("FAILED (Network error)\n")
    return(list(filename = fn, ok = FALSE, meta = NULL, msg = conditionMessage(e)))
  })
  
  if (!file.exists(tmp_dest)) {
    cat("FAILED (File not found)\n")
    return(list(filename = fn, ok = FALSE, meta = NULL, msg = "Download failed"))
  }
  
  # 2. Extract Metadata using Custom Logic
  meta_list <- parse_first_line_metadata(tmp_dest)
  
  # Convert to data.table row
  meta_row <- NULL
  if (!is.null(meta_list)) {
    meta_row <- as.data.table(meta_list)
    meta_row[, Filename := fn]
  }
  
  # 3. Read Data (Skip header)
  cols_to_keep <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", 
                    "cdr1_aa", "cdr2_aa", "cdr3_aa")
  
  dt <- tryCatch({
    data.table::fread(tmp_dest, skip = 1, select = cols_to_keep, showProgress = FALSE)
  }, error = function(e) NULL)
  
  # Cleanup Raw
  unlink(tmp_dest)
  
  if (is.null(dt)) {
    cat("FAILED (Corrupt CSV)\n")
    return(list(filename = fn, ok = FALSE, meta = meta_row, msg = "fread failed"))
  }
  
  # 4. Filter Empty Sequences
  dt <- dt[!is.na(sequence_alignment_aa) & sequence_alignment_aa != ""]
  
  # 5. Save Processed Data
  data.table::fwrite(dt, out_path, compress = "gzip")
  
  cat("OK\n")
  return(list(filename = fn, ok = TRUE, meta = meta_row, msg = "OK"))
}

# ==============================================================================
# 2. Main Execution
# ==============================================================================

main <- function() {
  
  sh_file   <- P$shell_script
  raw_dir   <- P$raw_dir
  proc_dir  <- P$processed_dir
  meta_file <- P$metadata_csv
  
  if (!file.exists(sh_file)) stop("Bulk download script not found at: ", sh_file)
  
  urls <- read_bulk_urls(sh_file)
  if (length(urls) == 0) stop("No URLs found!")
  
  target_files <- basename(urls)
  n_total      <- length(target_files)
  
  existing_files <- list.files(proc_dir, pattern = "\\.csv\\.gz$")
  
  cat(sprintf("Found %d targets. Starting pipeline...\n", n_total))
  
  skipped_count   <- 0
  processed_count <- 0
  failed_count    <- 0
  
  # --- Loop Start ---
  for (i in seq_len(n_total)) {
    url <- urls[i]
    fn  <- target_files[i]
    
    # 1. Skip Check
    if (fn %in% existing_files) {
      if (skipped_count %% 100 == 0) cat(sprintf("[%d/%d] Skipping %s\n", i, n_total, fn))
      skipped_count <- skipped_count + 1
      next
    }
    
    # 2. Process File
    res <- process_one_url(url, i, n_total, raw_dir, proc_dir)
    
    if (isTRUE(res$ok)) {
      processed_count <- processed_count + 1
      
      # 3. Append Metadata Immediately
      if (!is.null(res$meta)) {
        
        # Flatten list columns to characters
        res$meta <- res$meta[, lapply(.SD, function(col) {
          if (is.list(col)) return(sapply(col, toString))
          return(col)
        })]
        
        # Append to CSV
        if (file.exists(meta_file)) {
          fwrite(res$meta, meta_file, append = TRUE)
        } else {
          fwrite(res$meta, meta_file, append = FALSE)
        }
        
      } else {
        cat(sprintf("  [WARN] Metadata missing for %s\n", fn))
      }
      
    } else {
      failed_count <- failed_count + 1
      cat(sprintf("  [ERROR] Failed %s: %s\n", fn, res$msg))
    }
    
    if (i %% 50 == 0) gc()
  }
  
  # --- Finalize Metadata (Add Patient Column) ---
  # We do this at the end so ID generation is consistent across the whole dataset
  if (file.exists(meta_file)) {
    cat("\nFinalizing Metadata (Generating Patient IDs)...\n")
    metadata_df <- fread(meta_file)
    
    # Check if necessary columns exist
    if ("Author" %in% names(metadata_df) && "Subject" %in% names(metadata_df)) {
      
      # 1. Combine Author and Subject
      metadata_df[, Patient := paste(Author, Subject, sep = "_")]
      
      # 2. Create unique numeric identifier (P1, P2...)
      #    factor() assigns levels alphabetically, ensuring consistent IDs for same Author/Subject
      metadata_df[, Patient := paste0("P", as.numeric(factor(Patient)))]
      
      # Save back to CSV
      fwrite(metadata_df, meta_file)
      cat("Metadata successfully updated with 'Patient' column.\n")
      
    } else {
      cat("[WARN] 'Author' or 'Subject' columns missing. Skipping Patient ID generation.\n")
    }
  }
  
  cat("\n============================================================\n")
  cat("Pipeline Finished.\n")
  cat(sprintf(" - Processed:     %d\n", processed_count))
  cat(sprintf(" - Skipped:       %d\n", skipped_count))
  cat(sprintf(" - Failed:        %d\n", failed_count))
  cat("============================================================\n")
}

main()