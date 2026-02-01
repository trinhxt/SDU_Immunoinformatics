################################################################################
## scripts/04_build_dimension_table.R
##
## Goal: Build the Dimension Table (db_stage2)
##       - Aggregates stats per Peptide
##       - Columns: Peptide, N_Diseases, N_Patients, N_Antibodies, partition_id
##
## Improvements:
##   - Correctness: Uses exact distinct counts for Antibodies (no hashing).
##   - Safety: Uses bitwise partitioning to avoid integer overflows.
##   - Config: Uses centralized memory limits.
################################################################################

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(duckdb)
  library(arrow)
  library(dplyr)
})

main <- function() {
  
  INPUT_DB  <- file.path(P$parquet_dir, "db_stage1")
  OUTPUT_DB <- file.path(P$parquet_dir, "db_stage2")
  
  # Columns strictly as requested
  EXPECTED_COLS <- c(
    "Peptide", 
    "N_Diseases", 
    "N_Patients", 
    "N_Antibodies", 
    "partition_id"
  )
  
  cat("=== STEP 2: Building Dimension Table (db_stage2) ===\n")
  cat("Target DB:", OUTPUT_DB, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 1. Validation Checks
  # ---------------------------------------------------------------------------
  
  # Check Input Availability
  if (!dir.exists(INPUT_DB) || length(list.files(INPUT_DB, recursive = TRUE)) == 0) {
    stop("Input 'db_stage1' does not exist or is empty. Run scripts/03_build_fact_table.R first.")
  }
  
  # Check Output Existence & Schema
  if (dir.exists(OUTPUT_DB) && length(list.files(OUTPUT_DB, recursive = TRUE)) > 0) {
    cat("[INFO] Database directory exists. Checking schema...\n")
    
    tryCatch({
      ds <- arrow::open_dataset(OUTPUT_DB)
      actual_cols <- names(ds)
      missing_cols <- setdiff(EXPECTED_COLS, actual_cols)
      
      if (length(missing_cols) == 0) {
        cat("[SUCCESS] Database is valid. Rows:", format(nrow(ds), big.mark=","), ". Column:", names(ds), "\n")
        cat("Skipping rebuild.\n")
        return(invisible(NULL))
      } else {
        cat("[WARNING] Database exists but is missing columns:\n")
        cat(paste(missing_cols, collapse = ", "), "\n")
        cat("Proceeding with rebuild...\n\n")
      }
    }, error = function(e) {
      cat("[WARNING] Database directory exists but could not be opened.\n")
      cat("Proceeding with rebuild...\n\n")
    })
  }
  
  # ---------------------------------------------------------------------------
  # 2. Build Process (DuckDB)
  # ---------------------------------------------------------------------------
  
  # Temp space for DuckDB spilling
  TEMP_DIR  <- file.path(P$work_dir, "duckdb_spill")
  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # Persistent DB file (Vital for 2B rows to allow spilling)
  build_db_file <- file.path(P$work_dir, "build.duckdb")
  
  # Clean slate
  if (file.exists(build_db_file)) unlink(build_db_file)
  
  con <- dbConnect(duckdb(), dbdir = build_db_file)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
    if (file.exists(build_db_file)) unlink(build_db_file)
    if (dir.exists(TEMP_DIR)) unlink(TEMP_DIR, recursive = TRUE)
  }, add = TRUE)
  
  # Performance Settings
  # Use config memory limit if available, else default to 48GB
  mem_limit <- if(exists("DUCKDB_MEMORY_LIMIT")) DUCKDB_MEMORY_LIMIT else "48GB"
  
  dbExecute(con, paste0("PRAGMA memory_limit='", mem_limit, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  dbExecute(con, paste0("PRAGMA temp_directory='", gsub("\\\\", "/", TEMP_DIR), "'"))
  dbExecute(con, "PRAGMA preserve_insertion_order=false") # Improves parallelism
  
  input_glob <- paste0(gsub("\\\\", "/", INPUT_DB), "/**/*.parquet")
  output_dest <- gsub("\\\\", "/", OUTPUT_DB)
  
  cat("Aggregating 2 Billion rows... (Calculating Distinct Counts)\n")
  
  # SQL Logic:
  # 1. Exact Counts: We use COUNT(DISTINCT Antibody) instead of hash() to ensure scientific accuracy.
  # 2. Safe Partitioning: (hash(Peptide) & 255) ensures 0-255 range without integer overflow risks.
  
  sql <- sprintf("
    COPY (
      SELECT 
        Peptide,
        
        -- Aggregations
        COUNT(DISTINCT Disease)  AS N_Diseases,
        COUNT(DISTINCT Patient)  AS N_Patients,
        COUNT(DISTINCT Antibody) AS N_Antibodies, 
        
        -- Partitioning ID for Join Optimization (Bitwise AND is safer/faster than ABS/MOD)
        (hash(Peptide) & 255) AS partition_id
        
      FROM read_parquet('%s')
      GROUP BY Peptide
    )
    TO '%s'
    (FORMAT PARQUET, PARTITION_BY (partition_id), COMPRESSION ZSTD)
  ", input_glob, output_dest)
  
  dbExecute(con, sql)
  cat("Dimension Table Built: ", OUTPUT_DB, "\n")
}

main()
