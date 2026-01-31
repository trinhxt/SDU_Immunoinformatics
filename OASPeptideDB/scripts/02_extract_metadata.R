################################################################################
## scripts/02_extract_metadata.R
##
## Part 2: Extract JSON Metadata from Raw CSV.GZ files
##
## Logic based on provided 'Part1_Download.R':
## 1. Scan P$raw_dir for all *.csv.gz files.
## 2. Read ONLY the first line of each file (gzfile).
## 3. Apply specific Regex cleaning (NaN -> null, fix quotes).
## 4. Parse JSON -> Data Table.
## 5. Generate 'Patient' ID from Author + Subject.
## 6. Flatten list columns to strings.
## 7. Save to P$metadata_csv.
################################################################################

# ==============================================================================
# Load Config
# ==============================================================================
source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(future)
  library(future.apply)
})

# ==============================================================================
# Helper Functions (Based on your provided logic)
# ==============================================================================

# Function to read and parse the first line of a .csv.gz file
# Logic adapted strictly from user provided Part1_Download.R
get_metadata_row <- function(file_path) {
  
  # Safety check
  if (!file.exists(file_path)) return(NULL)
  
  tryCatch({
    # Open connection to gzip file
    con <- gzfile(file_path, "rt")
    first_line <- readLines(con, n = 1, warn = FALSE)
    close(con)
    
    if (length(first_line) == 0) return(NULL)
    
    # --- User Logic Start ---
    # Preprocess to replace NaN with null
    clean_json <- gsub('\\bNaN\\b', 'null', first_line)
    
    # Remove leading and trailing quotes, and replace double double-quotes with single double-quotes
    clean_json <- gsub('^"|"$', '', clean_json)
    clean_json <- gsub('""', '"', clean_json)
    # --- User Logic End ---
    
    # Parse JSON
    # flatten = TRUE helps convert nested JSON to columns immediately
    meta_list <- fromJSON(clean_json, flatten = TRUE)
    
    # Convert to data.table
    dt <- as.data.table(meta_list)
    
    # Add Filename (using the actual file on disk)
    dt[, Filename := basename(file_path)]
    
    return(dt)
    
  }, error = function(e) {
    # If a file is corrupt or empty, return NULL
    return(NULL)
  })
}

# ==============================================================================
# Main
# ==============================================================================

main <- function() {
  
  # 1. Identify Files
  raw_dir <- P$raw_dir
  if (!dir.exists(raw_dir)) stop("Raw directory not found: ", raw_dir)
  
  files <- list.files(raw_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  n_total <- length(files)
  
  if (n_total == 0) stop("No .csv.gz files found in ", raw_dir)
  
  cat("Found", n_total, "files in raw folder. Extracting metadata...\n")
  
  # 2. Extract Metadata (Parallel)
  # Using future_lapply for speed as we have ~13,000 files
  
  # Set up parallel plan based on config
  if (.Platform$OS.type == "windows") {
    plan(multisession, workers = N_CORES)
  } else {
    plan(multicore, workers = N_CORES)
  }
  
  # Run extraction
  meta_list <- future_lapply(seq_along(files), function(i) {
    get_metadata_row(files[i])
  })
  
  # Reset parallel plan
  plan(sequential)
  
  # 3. Combine Data
  cat("Combining results...\n")
  
  # Remove NULLs (failed files)
  meta_list <- meta_list[!sapply(meta_list, is.null)]
  
  if (length(meta_list) == 0) stop("No valid metadata extracted.")
  
  metadata_df <- rbindlist(meta_list, fill = TRUE)
  
  cat("Extracted metadata for", nrow(metadata_df), "files.\n")
  
  # 4. Post-Processing (User Logic)
  
  # Combine Author and Subject into a single identifier
  # Check if columns exist first to avoid crashing
  if (all(c("Author", "Subject") %in% names(metadata_df))) {
    cat("Generating Patient IDs...\n")
    metadata_df[, Patient_Temp := paste(Author, Subject, sep = "_")]
    # Create unique numeric identifier "P1", "P2", etc.
    metadata_df[, Patient := paste0("P", as.numeric(factor(Patient_Temp)))]
    metadata_df[, Patient_Temp := NULL]
  } else {
    cat("WARNING: 'Author' or 'Subject' columns missing. Skipping Patient ID generation.\n")
  }
  
  # Flatten List Columns
  # (Data.table cannot write list columns to CSV)
  cat("Flattening list columns...\n")
  
  # Identify columns that are lists
  list_cols <- names(which(sapply(metadata_df, is.list)))
  
  if (length(list_cols) > 0) {
    for (col in list_cols) {
      # Replace list with string representation
      # Using set() is faster than metadata_df[, col := ...]
      set(metadata_df, j = col, value = sapply(metadata_df[[col]], toString))
    }
  }
  
  # 5. Save Output
  out_file <- P$metadata_csv
  
  cat("Saving to:", out_file, "\n")
  fwrite(metadata_df, out_file)
  
  cat("Done.\n")
}

# Run
main()