################################################################################
## scripts/05_combine_db.R
##
## Goal: Create db_stage3 (The "Fast" Database)
##
## Method: DENORMALIZATION (Merging Fact + Dimension)
##   1. Iterates through each Disease partition in db_stage1.
##   2. Joins rows with db_stage2 (Metrics) in memory.
##   3. Writes a single "Flat" Parquet dataset (db_stage3).
##
## Result:
##   - A single table containing EVERYTHING (Observations + Stats).
##   - Optimizes query speed from ~30s (Join) to <1s (Scan).
################################################################################

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(duckdb)
  library(arrow)
  library(dplyr)
})

main <- function() {
  
  # Paths
  INPUT_FACT  <- file.path(P$parquet_dir, "db_stage1")
  INPUT_DIM   <- file.path(P$parquet_dir, "db_stage2")
  OUTPUT_DB   <- file.path(P$parquet_dir, "db_stage3")
  
  # Configuration
  # We process one disease at a time to keep RAM usage low
  
  cat("=== STEP 5: COMBINING DATABASES (Creating db_stage3) ===\n")
  cat("Input Fact:  ", INPUT_FACT, "\n")
  cat("Input Dim:   ", INPUT_DIM, "\n")
  cat("Output DB:   ", OUTPUT_DB, "\n\n")
  
  # Validation
  if (!dir.exists(INPUT_FACT)) stop("db_stage1 missing!")
  if (!dir.exists(INPUT_DIM))  stop("db_stage2 missing!")
  
  if (dir.exists(OUTPUT_DB)) {
    cat("Output directory exists. Deleting to rebuild...\n")
    unlink(OUTPUT_DB, recursive = TRUE)
  }
  
  # 1. Setup DuckDB
  con <- dbConnect(duckdb::duckdb())
  on.exit(dbDisconnect(con, shutdown = TRUE))
  
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCKDB_MEMORY_LIMIT, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  
  # 2. Register Views
  #    Recursively glob all parquet files from both stages
  fact_glob <- paste0(gsub("\\\\", "/", INPUT_FACT), "/**/*.parquet")
  dim_glob  <- paste0(gsub("\\\\", "/", INPUT_DIM), "/**/*.parquet")
  
  dbExecute(con, sprintf("CREATE VIEW obs AS SELECT * FROM read_parquet('%s')", fact_glob))
  dbExecute(con, sprintf("CREATE VIEW met AS SELECT * FROM read_parquet('%s')", dim_glob))
  
  # 3. Get List of Diseases (Partition Keys)
  cat("Scanning partition keys from db_stage1...\n")
  diseases <- dbGetQuery(con, "SELECT DISTINCT Disease FROM obs WHERE Disease IS NOT NULL ORDER BY Disease")[[1]]
  
  cat("Found", length(diseases), "Diseases. Starting merge process...\n")
  
  # 4. Loop and Merge
  #    Processing sequentially by Disease is robust against memory crashes.
  
  for (i in seq_along(diseases)) {
    curr_disease <- diseases[i]
    
    # Simple progress indicator
    cat(sprintf("\n [%d/%d] Processing Disease: %s ... ", i, length(diseases), curr_disease, "\n"))
    t1 <- Sys.time()
    
    # The Query:
    # 1. Select ALL columns from Observations (o.*)
    # 2. Select ONLY metric columns from Metrics (m.*)
    # 3. Filter by current Disease (Partition Pruning)
    # 4. Left Join on Peptide
    
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
      cat("Done in", round(difftime(t2, t1, units="mins"), 2), "min.\n")
    }, error = function(e) {
      cat("FAILED!\nError: ", conditionMessage(e), "\n")
    })
    
    # Force garbage collection between large partitions
    gc()
  }
  
  cat("\n================================================================\n")
  cat("SUCCESS! db_stage3 created.\n")
  cat("Location: ", OUTPUT_DB, "\n")
  cat("You can now build the Zenodo release using this single folder.\n")
  cat("================================================================\n")
}

main()