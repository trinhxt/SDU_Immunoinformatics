################################################################################
## scripts/03_build_fact_table.R
##
## Goal: Build the Fact Table (db_stage1)
##       - Joins Digest Data with Metadata
##       - Partitions by Clinical Variables
##       - Calculates CDR3 Overlap (Amino Acid Count) per row
##
## Optimizations:
##   - Vectorized Overlap Calc (stringi)
##   - Flat Staging (Reduced I/O overhead)
##   - Defensive Metadata Joins
################################################################################

source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(dplyr)
  library(duckdb)
  library(future)
  library(future.apply)
  library(stringi) # C-optimized string operations
})

# Increase future globals limit for large objects (just in case)
options(future.globals.maxSize = 2000 * 1024^2) # 2GB

main <- function() {
  
  # Paths
  STAGING_DIR <- file.path(P$parquet_dir, "_staging")
  DB_STAGE1   <- file.path(P$parquet_dir, "db_stage1")
  
  # Partition Keys (Applied by DuckDB at the end)
  PARTITION_COLS <- c("Disease", "BSource", "BType", "Isotype")
  
  # Expected Columns
  EXPECTED_COLS <- c(
    "Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa", 
    "filename", "Patient", "Disease", "BSource", "BType", "Isotype",
    "CDR3_spanning_count"
  )
  
  cat("=== STEP 1: Building Fact Table (db_stage1) ===\n")
  cat("Target DB:", DB_STAGE1, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 1. Availability & Schema Check
  # ---------------------------------------------------------------------------
  if (dir.exists(DB_STAGE1) && length(list.files(DB_STAGE1, recursive = TRUE)) > 0) {
    cat("[INFO] Database directory exists. Checking schema...\n")
    tryCatch({
      ds <- arrow::open_dataset(DB_STAGE1)
      missing_cols <- setdiff(EXPECTED_COLS, names(ds))
      if (length(missing_cols) == 0) {
        cat("[SUCCESS] Database is valid. Row: ", format(nrow(ds), big.mark=","), ". Column:", names(ds), "\n" , "Skipping rebuild.")
        return(invisible(NULL))
      } else {
        cat("[WARNING] Missing columns:", paste(missing_cols, collapse=", "), "\n")
      }
    }, error = function(e) cat("[WARNING] DB Check failed. Rebuilding.\n"))
  }
  
  # ---------------------------------------------------------------------------
  # 2. Setup
  # ---------------------------------------------------------------------------
  
  # Load Metadata
  if (!file.exists(P$metadata_csv)) stop("Metadata missing")
  
  # Load and Key Metadata
  meta <- fread(P$metadata_csv, select = c("Filename", "Patient", "Disease", "BSource", "BType", "Isotype"))
  meta <- meta[!is.na(Filename)]
  setkey(meta, Filename)
  
  # Prepare Task List
  files <- list.files(P$digest_csv_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  tasks <- data.table(path = files, filename = basename(files))
  tasks <- tasks[filename %in% meta$Filename]
  
  cat("Processing", nrow(tasks), "files...\n")
  
  # Clean/Create Staging
  if (dir.exists(STAGING_DIR)) unlink(STAGING_DIR, recursive = TRUE)
  dir.create(STAGING_DIR, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 3. Worker Function (Optimized)
  # ---------------------------------------------------------------------------
  stage_worker <- function(fp, fn, meta_row, out_dir) {
    # Fast Read
    dt <- fread(fp)
    if (nrow(dt) == 0) return(NULL)
    
    # A. Defensive Metadata Merge
    if (nrow(meta_row) != 1) {
      stop(sprintf("Metadata error: %s maps to %d rows (must be 1)", fn, nrow(meta_row)))
    }
    
    # B. Vectorized Overlap Calculation (stringi)
    #    This is orders of magnitude faster than regexpr loop
    if (all(c("Antibody", "Peptide", "cdr3_aa") %in% names(dt))) {
      
      # Get location matrices [start, end]
      # Returns NA if not found
      p_loc <- stringi::stri_locate_first_fixed(dt$Antibody, dt$Peptide)
      c_loc <- stringi::stri_locate_first_fixed(dt$Antibody, dt$cdr3_aa)
      
      p_starts <- p_loc[,1]; p_ends <- p_loc[,2]
      c_starts <- c_loc[,1]; c_ends <- c_loc[,2]
      
      # Identify invalid rows (missing patterns)
      invalid <- is.na(p_starts) | is.na(c_starts)
      
      # Calculate Overlap: Intersection Length
      # Formula: max(0, min(EndA, EndB) - max(StartA, StartB) + 1)
      overlap <- pmax(0L, pmin(p_ends, c_ends) - pmax(p_starts, c_starts) + 1L)
      
      # Zero out invalids
      overlap[invalid] <- 0L
      
      dt[, CDR3_spanning_count := overlap]
    } else {
      dt[, CDR3_spanning_count := 0L]
    }
    
    # C. Attach Metadata
    dt[, `:=`(
      Patient = as.character(meta_row$Patient),
      Disease = as.character(meta_row$Disease),
      BSource = as.character(meta_row$BSource),
      BType   = as.character(meta_row$BType),
      Isotype = as.character(meta_row$Isotype),
      filename = fn
    )]
    
    # D. Flat Write (Faster I/O)
    #    Write 1 parquet file per input. No partitioning yet.
    #    We let DuckDB handle the partitioning shuffle later.
    out_file <- file.path(out_dir, paste0(fn, ".parquet"))
    arrow::write_parquet(dt, out_file, compression = "zstd")
    
    return(TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # 4. Execution (Parallel)
  # ---------------------------------------------------------------------------
  # Note: If memory issues occur, reduce workers = N_CORES to workers = 4
  plan(multisession, workers = N_CORES) 
  
  res <- future_lapply(seq_len(nrow(tasks)), function(i) {
    # Defensive lookup
    m <- meta[J(tasks$filename[i])]
    stage_worker(tasks$path[i], tasks$filename[i], m, STAGING_DIR)
  })
  
  plan(sequential)
  
  # ---------------------------------------------------------------------------
  # 5. Compaction & Partitioning (DuckDB)
  # ---------------------------------------------------------------------------
  cat("Compacting and Partitioning into Final DB...\n")
  
  con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
  on.exit(dbDisconnect(con, shutdown = TRUE))
  
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCKDB_MEMORY_LIMIT, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  
  # Read from FLAT staging files
  src_glob <- paste0(gsub("\\\\", "/", STAGING_DIR), "/*.parquet")
  dest_path <- gsub("\\\\", "/", DB_STAGE1)
  
  # DuckDB handles the heavy lifting of partitioning here
  sql <- sprintf("
    COPY (SELECT * FROM read_parquet('%s'))
    TO '%s'
    (FORMAT PARQUET, PARTITION_BY (%s), COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000)
  ", src_glob, dest_path, paste(PARTITION_COLS, collapse = ", "))
  
  tryCatch({
    dbExecute(con, sql)
    
    # Cleanup only on success
    cat("Success. Removing staging files...\n")
    unlink(STAGING_DIR, recursive = TRUE)
    cat("Fact Table Built: ", DB_STAGE1, "\n")
    
  }, error = function(e) {
    cat("[FATAL] DuckDB Compaction failed: ", conditionMessage(e), "\n")
    cat("Staging files preserved in: ", STAGING_DIR, "\n")
  })
}

main()