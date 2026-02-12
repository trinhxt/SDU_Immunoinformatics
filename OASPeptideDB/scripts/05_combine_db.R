################################################################################
## scripts/05_combine_db.R
##
## Goal: Create db_stage3 (The "Fast" Database)
##
## Workflow:
##   1. Scans db_stage1 (Facts) and db_stage2 (Dimensions).
##   2. Identifies all unique Diseases in db_stage1.
##   3. Iterates through each Disease:
##      a. Loads Facts for that Disease (Partition Pruning).
##      b. Joins with db_stage2 on 'Peptide' to attach stats & flags.
##      c. Writes Denormalized result to db_stage3.
##
## Output Schema (db_stage3):
##   - All Clinical Metadata (Patient, Disease, BSource, etc.)
##   - All Molecular Data (Peptide, CDR3_spanning_count, etc.)
##   - Global Stats (N_Patients, N_Diseases, etc.)
##   - Healthy Background Flag (is_healthy_background)
##
## Robustness:
##   - Resume Capability: Skips Diseases that already exist in output.
##   - Memory Safety: Processes one Disease cohort at a time.
################################################################################

if (!exists("P")) stop("Configuration 'P' not found. Please run this via DBbuild.R")

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(arrow)
  library(dplyr)
  library(utils) # for URLdecode
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # 1. Setup & Configuration
  # ---------------------------------------------------------------------------
  INPUT_FACT  <- file.path(P$parquet_dir, "db_stage1")
  INPUT_DIM   <- file.path(P$parquet_dir, "db_stage2")
  OUTPUT_DB   <- file.path(P$parquet_dir, "db_stage3")
  
  # Memory Limit (Safe default or Config)
  DUCK_MEM <- if(exists("DUCKDB_MEMORY_LIMIT")) DUCKDB_MEMORY_LIMIT else "48GB"
  
  cat("=== STEP 5: COMBINING DATABASES (Creating db_stage3) ===\n")
  cat("Input Fact:  ", INPUT_FACT, "\n")
  cat("Input Dim:   ", INPUT_DIM, "\n")
  cat("Output DB:   ", OUTPUT_DB, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 2. Input Validation (Check Stage 1 & 2 Columns)
  # ---------------------------------------------------------------------------
  cat("[1/4] Validating Inputs...\n")
  
  if (!dir.exists(INPUT_FACT)) stop("FATAL: db_stage1 missing! Run script 03.")
  if (!dir.exists(INPUT_DIM))  stop("FATAL: db_stage2 missing! Run script 04.")
  
  # Define Expected Columns based on previous steps
  # Facts (from 03_build_fact_table.R)
  COLS_FACT <- c("Peptide", "Antibody", "Disease", "Patient", "BSource", "BType", "Isotype", 
                 "CDR3_spanning_count", "CDR3_spanning_pct", "v_call", "d_call", "j_call", "cdr3_aa")
  
  # Dimensions (from 04_build_dimension_table.R)
  COLS_DIM  <- c("Peptide", "N_Patients", "N_Diseases", "N_Antibodies", "is_healthy_background")
  
  # Check Fact Table
  ds_fact <- arrow::open_dataset(INPUT_FACT)
  missing_fact <- setdiff(COLS_FACT, names(ds_fact))
  if (length(missing_fact) > 0) {
    stop("FATAL: db_stage1 is missing columns: ", paste(missing_fact, collapse=", "))
  }
  
  # Check Dimension Table
  ds_dim <- arrow::open_dataset(INPUT_DIM)
  missing_dim <- setdiff(COLS_DIM, names(ds_dim))
  if (length(missing_dim) > 0) {
    stop("FATAL: db_stage2 is missing columns: ", paste(missing_dim, collapse=", "))
  }
  
  cat("      Inputs are valid.\n")
  
  # ---------------------------------------------------------------------------
  # 3. Setup DuckDB & Views
  # ---------------------------------------------------------------------------
  cat("[2/4] Initializing Database Engine...\n")
  
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE))
  
  # Optimize for high-throughput joins
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCK_MEM, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  dbExecute(con, "PRAGMA preserve_insertion_order=false")
  
  # Register Views (Recursive Glob to catch all partitioned files)
  fact_glob <- paste0(gsub("\\\\", "/", INPUT_FACT), "/**/*.parquet")
  dim_glob  <- paste0(gsub("\\\\", "/", INPUT_DIM), "/**/*.parquet")
  
  dbExecute(con, sprintf("CREATE OR REPLACE VIEW obs AS SELECT * FROM read_parquet('%s')", fact_glob))
  dbExecute(con, sprintf("CREATE OR REPLACE VIEW met AS SELECT * FROM read_parquet('%s')", dim_glob))
  
  # Get Source of Truth: List of Diseases from Fact Table
  source_diseases <- dbGetQuery(con, "SELECT DISTINCT Disease FROM obs WHERE Disease IS NOT NULL ORDER BY Disease")[[1]]
  cat("      Found", length(source_diseases), "unique Diseases in source.\n")
  
  # ---------------------------------------------------------------------------
  # 4. Check Existing Output (Resume Logic)
  # ---------------------------------------------------------------------------
  cat("[3/4] Checking Status of db_stage3...\n")
  
  # We check the directory structure of OUTPUT_DB for "Disease=X" folders
  completed_diseases <- character(0)
  
  if (dir.exists(OUTPUT_DB)) {
    # List directories like "Disease=Covid", "Disease=Healthy", etc.
    dirs <- list.dirs(OUTPUT_DB, recursive = TRUE, full.names = FALSE)
    # Extract "X" from "Disease=X"
    disease_dirs <- grep("^Disease=", basename(dirs), value = TRUE)
    encoded_names <- gsub("^Disease=", "", disease_dirs)
    
    # FIX: Decode URL-encoded names (e.g. "Allergy%2FNoSIT" -> "Allergy/NoSIT")
    # This ensures they match the actual disease names in 'source_diseases'
    completed_diseases <- sapply(encoded_names, function(x) utils::URLdecode(x), USE.NAMES = FALSE)
  }
  
  # Calculate what is left to do
  missing_diseases <- setdiff(source_diseases, completed_diseases)
  
  if (length(missing_diseases) == 0) {
    cat("\n[SUCCESS] db_stage3 is already complete! (All", length(source_diseases), "diseases present).\n")
    cat("Skipping build.\n")
    return(invisible(NULL))
  }
  
  # List Status
  if (length(completed_diseases) > 0) {
    cat("      Already Done (", length(completed_diseases), "):\n", sep="")
    cat(paste0("       - ", completed_diseases), sep = "\n")
  } else {
    cat("      Already Done: 0\n")
  }
  
  if (length(missing_diseases) > 0) {
    cat("      Remaining (", length(missing_diseases), "):\n", sep="")
    cat(paste0("       - ", missing_diseases), sep = "\n")
  }
  
  # ---------------------------------------------------------------------------
  # 5. Build Process (Resume Loop)
  # ---------------------------------------------------------------------------
  cat("[4/4] Starting/Resuming Merge Process...\n")
  
  # Loop through diseases one by one to keep memory usage stable
  for (i in seq_along(missing_diseases)) {
    curr_disease <- missing_diseases[i]
    
    cat(sprintf("\n [%d/%d] Processing: %s ... \n", i, length(missing_diseases), curr_disease))
    t1 <- Sys.time()
    
    # The Join Query:
    # 1. Filter 'obs' (Facts) to just the current disease.
    # 2. Left Join 'met' (Dimensions) on Peptide.
    # 3. Select all relevant columns.
    
    # Note on Output:
    # We use PARTITION_BY to maintain the clinical hierarchy.
    # DuckDB will automatically URL-encode special chars in partition folders (e.g. "/" -> "%2F").
    
    sql <- sprintf("
      COPY (
        SELECT 
          o.*, 
          m.N_Patients,
          m.N_Diseases,
          m.N_Antibodies,
          m.is_healthy_background
        FROM obs o
        LEFT JOIN met m ON o.Peptide = m.Peptide
        WHERE o.Disease = '%s'
      ) 
      TO '%s' 
      (FORMAT PARQUET, 
       PARTITION_BY (Disease, BSource, BType, Isotype), 
       COMPRESSION ZSTD, 
       ROW_GROUP_SIZE 1000000, 
       OVERWRITE_OR_IGNORE)
    ", curr_disease, gsub("\\\\", "/", OUTPUT_DB))
    
    tryCatch({
      dbExecute(con, sql)
      t2 <- Sys.time()
      cat("Done in", round(difftime(t2, t1, units="mins"), 2), "min.")
      
      # Explicit Checkpoint to free WAL memory
      dbExecute(con, "CHECKPOINT")
      
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
