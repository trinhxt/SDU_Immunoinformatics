################################################################################
## scripts/05_combine_db.R
##
## Goal: Create db_stage3 (The "Fast" Database)
##
## Features:
##   - Validation: Verifies stage1/stage2 columns before starting.
##   - Resumable: If db_stage3 exists but is partial, it resumes (skips done).
##   - Robust: Handles crashes or interruptions without data loss.
##
## Method: DENORMALIZATION (Merging Fact + Dimension)
##   1. Joins 'db_stage1' (Observations) with 'db_stage2' (Metrics).
##   2. Writes to 'db_stage3' (Single Flat Table).
################################################################################

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(arrow)
  library(dplyr)
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # 1. Setup & Configuration
  # ---------------------------------------------------------------------------
  INPUT_FACT  <- file.path(P$parquet_dir, "db_stage1")
  INPUT_DIM   <- file.path(P$parquet_dir, "db_stage2")
  OUTPUT_DB   <- file.path(P$parquet_dir, "db_stage3")
  
  # Memory Limit (Safe default or Config)
  DUCK_MEM <- if(exists("DUCKDB_MEMORY_LIMIT")) DUCKDB_MEMORY_LIMIT else "24GB"
  
  cat("=== STEP 5: COMBINING DATABASES (Creating db_stage3) ===\n")
  cat("Input Fact:  ", INPUT_FACT, "\n")
  cat("Input Dim:   ", INPUT_DIM, "\n")
  cat("Output DB:   ", OUTPUT_DB, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 2. Input Validation (Check Stage 1 & 2 Columns)
  # ---------------------------------------------------------------------------
  cat("[1/4] Validating Inputs...\n")
  
  if (!dir.exists(INPUT_FACT)) stop("FATAL: db_stage1 missing!")
  if (!dir.exists(INPUT_DIM))  stop("FATAL: db_stage2 missing!")
  
  # Define Expected Columns
  COLS_FACT <- c("Peptide", "Antibody", "Disease", "CDR3_spanning_count", "Patient", "Isotype")
  COLS_DIM  <- c("Peptide", "N_Patients", "N_Diseases", "N_Antibodies")
  
  # Check Fact Table
  ds_fact <- arrow::open_dataset(INPUT_FACT)
  if (!all(COLS_FACT %in% names(ds_fact))) {
    stop("FATAL: db_stage1 is missing required columns: ", 
         paste(setdiff(COLS_FACT, names(ds_fact)), collapse=", "))
  }
  
  # Check Dimension Table
  ds_dim <- arrow::open_dataset(INPUT_DIM)
  if (!all(COLS_DIM %in% names(ds_dim))) {
    stop("FATAL: db_stage2 is missing required columns: ", 
         paste(setdiff(COLS_DIM, names(ds_dim)), collapse=", "))
  }
  
  cat("      Inputs are valid.\n")
  
  # ---------------------------------------------------------------------------
  # 3. Setup DuckDB & Views
  # ---------------------------------------------------------------------------
  cat("[2/4] Initializing Database Engine...\n")
  
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE))
  
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCK_MEM, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  
  # Register Views (Recursive Glob)
  fact_glob <- paste0(gsub("\\\\", "/", INPUT_FACT), "/**/*.parquet")
  dim_glob  <- paste0(gsub("\\\\", "/", INPUT_DIM), "/**/*.parquet")
  
  dbExecute(con, sprintf("CREATE VIEW obs AS SELECT * FROM read_parquet('%s')", fact_glob))
  dbExecute(con, sprintf("CREATE VIEW met AS SELECT * FROM read_parquet('%s')", dim_glob))
  
  # Get Source of Truth: List of Diseases from Fact Table
  source_diseases <- dbGetQuery(con, "SELECT DISTINCT Disease FROM obs WHERE Disease IS NOT NULL ORDER BY Disease")[[1]]
  cat("      Found", length(source_diseases), "unique Diseases in source.\n")
  
  # ---------------------------------------------------------------------------
  # 4. Check Existing Output (Resume Logic)
  # ---------------------------------------------------------------------------
  cat("[3/4] Checking Status of db_stage3...\n")
  
  completed_diseases <- character(0)
  
  if (dir.exists(OUTPUT_DB)) {
    tryCatch({
      # Use Arrow to read partition keys from disk
      ds_out <- arrow::open_dataset(OUTPUT_DB)
      
      # Check if output has expected columns (if not, force rebuild)
      # We check for merged columns (Observation + Metric)
      if (!all(c("N_Patients", "Peptide", "Disease") %in% names(ds_out))) {
        cat("      [WARN] Existing db_stage3 corrupted or missing columns. Forcing Rebuild.\n")
        unlink(OUTPUT_DB, recursive = TRUE)
      } else {
        completed_diseases <- ds_out %>% 
          select(Disease) %>% 
          distinct() %>% 
          collect() %>% 
          pull(Disease)
      }
    }, error = function(e) {
      cat("      [WARN] Could not read existing db_stage3 (likely empty). Starting fresh.\n")
    })
  }
  
  # Calculate what is left to do
  missing_diseases <- setdiff(source_diseases, completed_diseases)
  
  if (length(missing_diseases) == 0) {
    cat("\n[SUCCESS] db_stage3 is already complete! (All", length(source_diseases), "diseases present).\n")
    cat("Skipping build.\n")
    return(invisible(NULL))
  }
  
  cat("      Already Done: ", length(completed_diseases), "\n")
  cat("      Remaining:    ", length(missing_diseases), "\n")
  
  # ---------------------------------------------------------------------------
  # 5. Build Process (Resume Loop)
  # ---------------------------------------------------------------------------
  cat("[4/4] Starting/Resuming Merge Process...\n")
  
  for (i in seq_along(missing_diseases)) {
    curr_disease <- missing_diseases[i]
    
    cat(sprintf("\n [%d/%d] Processing: %s ... ", i, length(missing_diseases), curr_disease))
    t1 <- Sys.time()
    
    # The Query:
    # 1. Select ALL columns from Observations (o.*)
    # 2. Select ONLY metric columns from Metrics (m.*)
    # 3. Filter by current Disease (Partition Pruning)
    # 4. Left Join on Peptide
    
    # OVERWRITE_OR_IGNORE handles re-runs safely
    
    sql <- sprintf("
      COPY (
        SELECT 
          o.*, 
          m.N_Patients,
          m.N_Diseases,
          m.N_Antibodies
        FROM obs o
        LEFT JOIN met m ON o.Peptide = m.Peptide
        WHERE o.Disease = '%s'
      ) 
      TO '%s' 
      (FORMAT PARQUET, PARTITION_BY (Disease, BSource, BType, Isotype), COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000, OVERWRITE_OR_IGNORE)
    ", curr_disease, gsub("\\\\", "/", OUTPUT_DB))
    
    tryCatch({
      dbExecute(con, sql)
      t2 <- Sys.time()
      cat("Done in", round(difftime(t2, t1, units="mins"), 2), "min.")
    }, error = function(e) {
      cat("\n[FAILED] Error processing", curr_disease, ":", conditionMessage(e), "\n")
      cat("Note: This error is isolated. You can re-run this script to retry just this disease.\n")
    })
    
    # Garbage Collection to keep RAM stable over long runs
    gc()
  }
  
  cat("\n================================================================\n")
  cat("SUCCESS! db_stage3 build complete.\n")
  cat("Location: ", OUTPUT_DB, "\n")
  cat("================================================================\n")
}

main()