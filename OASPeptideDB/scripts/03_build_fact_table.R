################################################################################
## scripts/03_build_fact_table.R
##
## Goal: Build the Fact Table (db_stage1)
##       - Aggregates "Disease" Parquet files from Stage 2.
##       - Partitions data by Clinical Variables for fast querying.
##       - EXCLUDES Healthy (Background) data to maintain schema consistency.
##
## Optimization:
##   - Direct Parquet->DuckDB->Parquet pipeline (Zero R memory overhead).
##   - Bulk loading of file lists.
################################################################################

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dplyr)
  library(duckdb)
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # 1. Configuration & Paths
  # ---------------------------------------------------------------------------
  INPUT_DIR   <- P$digest_csv_dir      # Now contains .parquet files from step 02
  DB_STAGE1   <- file.path(P$parquet_dir, "db_stage1")
  
  # Partition Keys (How the data will be organized on disk)
  # partitioning by Disease/BSource makes querying specific cohorts instant.
  PARTITION_COLS <- c("Disease", "BSource", "BType", "Isotype")
  
  cat("=== STEP 3: Building Fact Table (db_stage1) ===\n")
  cat("Input Dir: ", INPUT_DIR, "\n")
  cat("Target DB: ", DB_STAGE1, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 2. Availability Check
  # ---------------------------------------------------------------------------
  # If DB exists, check if it looks valid to avoid accidental overwrites
  if (dir.exists(DB_STAGE1) && length(list.files(DB_STAGE1, recursive = TRUE)) > 0) {
    cat("[INFO] Database directory exists. Checking schema...\n")
    tryCatch({
      ds <- arrow::open_dataset(DB_STAGE1)
      # Check for a critical column that only exists in processed Disease data
      if ("CDR3_spanning_count" %in% names(ds)) {
        cat("[SUCCESS] Database appears valid. Rows:", format(nrow(ds), big.mark=","), "\n")
        cat("To force rebuild, delete the 'db_stage1' folder manually.\n")
        return(invisible(NULL))
      }
    }, error = function(e) cat("[WARNING] DB validation failed. Proceeding with rebuild.\n"))
  }
  
  # ---------------------------------------------------------------------------
  # 3. File Selection (CRITICAL: Filter out Healthy)
  # ---------------------------------------------------------------------------
  cat("Loading Metadata to select DISEASE files...\n")
  
  if (!file.exists(P$metadata_csv)) stop("Metadata missing")
  
  # Load Metadata
  meta <- fread(P$metadata_csv, select = c("Filename", "Disease"))
  
  # 1. Get all parquet files in the digest folder
  all_files <- list.files(INPUT_DIR, pattern = "\\.parquet$", full.names = TRUE)
  file_dt   <- data.table(path = all_files, Filename = gsub("\\.parquet$", ".csv.gz", basename(all_files)))
  
  # 2. Join with Metadata
  #    Note: filenames in 'file_dt' derived from parquet need to match metadata 'Filename'
  #    Adjust logic if your naming convention in 02 changed slightly.
  #    Assuming 02 output: "filename.parquet" matches "filename.csv.gz" in metadata
  tasks <- merge(file_dt, meta, by = "Filename")
  
  # 3. FILTER: Keep ONLY Disease Samples
  #    Healthy samples (Disease == "None") have a different schema (Peptide only)
  #    and cannot be merged into the Fact Table.
  disease_tasks <- tasks[Disease != "None"]
  
  cat("Total Files Found:     ", nrow(tasks), "\n")
  cat("Healthy Files (Skip):  ", nrow(tasks) - nrow(disease_tasks), "\n")
  cat("Disease Files (Merge): ", nrow(disease_tasks), "\n")
  
  if (nrow(disease_tasks) == 0) stop("No Disease files found to process!")
  
  # ---------------------------------------------------------------------------
  # 4. DuckDB Construction
  # ---------------------------------------------------------------------------
  cat("\nInitializing DuckDB for assembly...\n")
  
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE))
  
  # Memory & Thread Management
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCKDB_MEMORY_LIMIT, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  
  # ---------------------------------------------------------------------------
  # 5. Execution: Bulk Copy
  # ---------------------------------------------------------------------------
  # Strategy: We pass the vector of filenames directly to DuckDB's read_parquet.
  # This is much faster than creating a view or looping.
  
  # Create a temporary table containing the list of files to read
  # DuckDB can read a list of files provided as a string list
  # We format it as explicit SQL list: ['file1', 'file2', ...]
  
  # Windows path safety
  safe_paths <- gsub("\\\\", "/", disease_tasks$path)
  
  # Create the output directory fresh
  if (dir.exists(DB_STAGE1)) unlink(DB_STAGE1, recursive = TRUE)
  dir.create(DB_STAGE1, recursive = TRUE)
  
  cat("Streaming", length(safe_paths), "files into Partitioned Parquet DB...\n")
  
  # We use a temporary view to inspect schema if needed, but COPY is most efficient.
  # Note: partition_by(...) requires the columns to be present in the input.
  # 02_digestion.R ensures 'Disease', 'BSource', etc. are in the parquet files.
  
  query <- sprintf("
    COPY (
      SELECT * FROM read_parquet([%s])
    )
    TO '%s'
    (FORMAT PARQUET, PARTITION_BY (%s), COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000)
  ", 
                   paste(paste0("'", safe_paths, "'"), collapse = ", "), # 'path1', 'path2'
                   gsub("\\\\", "/", DB_STAGE1),
                   paste(PARTITION_COLS, collapse = ", ")
  )
  
  tryCatch({
    # Execute the massive copy
    dbExecute(con, query)
    
    cat("\n[SUCCESS] Fact Table Built Successfully.\n")
    cat("Location: ", DB_STAGE1, "\n")
    
  }, error = function(e) {
    cat("\n[FATAL] DuckDB Copy failed: ", conditionMessage(e), "\n")
    stop("Build failed.")
  })
}

main()
