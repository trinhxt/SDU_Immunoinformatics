################################################################################
## scripts/04_build_dimension_table.R
##
## Goal: Build the Dimension Table (db_stage2) as partitioned Parquet
##
## Required output columns (one row per DISTINCT disease peptide):
##   - Peptide (string)
##   - N_Diseases   = COUNT(DISTINCT Disease)
##   - N_Patients   = COUNT(DISTINCT Patient)
##   - N_Antibodies = COUNT(DISTINCT Antibody)
##   - is_healthy_background (TRUE/FALSE)  [STRICT / EXACT membership]
##   - partition_id (0–127) derived from peptide hash
##
## STRICT healthy label means:
##   is_healthy_background = TRUE  iff the exact Peptide string occurs in healthy
##   background files. No hash-based approximation.
##
## High-level approach (strict & scalable):
##   1) disease_stats: GROUP BY Peptide from disease_facts, compute counts + partition_id.
##      Write as partitioned parquet: disease_stats_part/partition_id=...
##
##   2) healthy batches: stream healthy parquet files in batches.
##      For each batch, write DISTINCT Peptide + partition_id to its OWN output dir:
##        healthy_batches_part/batch_<b>/batch_id=<b>/partition_id=...
##      (Important: DuckDB COPY refuses writing into a non-empty dataset directory;
##       per-batch subdirectories avoid that.)
##
##   3) healthy index per partition: for each partition_id, read all batch outputs
##      for that partition, run SELECT DISTINCT Peptide, write one parquet file:
##        healthy_index_part/partition_id=<pid>/part-<pid>.parquet
##
##   4) final output per partition: strict string join within each partition_id:
##        disease_stats(pid) LEFT JOIN healthy_index(pid) ON Peptide
##      Write db_stage2/partition_id=<pid>/part-<pid>.parquet
##
## Resume support:
##   - Set RESUME <- TRUE to skip work that already exists on disk.
##   - Checks final DB partitions and schema validity before starting.
################################################################################

if (!exists("P")) stop("Configuration 'P' not found. Please run this via DBbuild.R")

suppressPackageStartupMessages({
  library(duckdb)
  library(arrow)
  library(dplyr)
  library(data.table)
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # 0. Resume / Restart controls
  # ---------------------------------------------------------------------------
  # RESUME = TRUE  -> keep intermediate outputs; skip completed chunks
  # RESUME = FALSE -> delete intermediate outputs and rebuild everything
  RESUME <- TRUE
  
  # If you want a "clean rebuild" but keep stage1, set RESUME <- FALSE.
  
  # ---------------------------------------------------------------------------
  # 1. Configuration & Paths
  # ---------------------------------------------------------------------------
  FACT_DB     <- file.path(P$parquet_dir, "db_stage1")
  DIGEST_DIR  <- P$digest_csv_dir
  OUTPUT_DB   <- file.path(P$parquet_dir, "db_stage2")
  META_FILE   <- P$metadata_csv
  
  # Healthy scan batching constraints
  MAX_BATCH_FILES <- 2000
  MAX_BATCH_SIZE  <- 1 * 1024^3  # 1GB parquet bytes per batch (tune if needed)
  
  # Partitioning requirement: 0–127 (128 partitions)
  PART_MASK <- 127L
  N_PARTS   <- 128L
  
  cat("=== STEP 4: Building Dimension Table (db_stage2) [STRICT HEALTHY LABEL] ===\n")
  cat("Fact Table:   ", FACT_DB, "\n")
  cat("Healthy Src:  ", DIGEST_DIR, "\n")
  cat("Target DB:    ", OUTPUT_DB, "\n")
  cat("Partitions:   0–", PART_MASK, " (", N_PARTS, " total)\n", sep = "")
  cat("Resume mode:  ", ifelse(RESUME, "ON", "OFF (clean rebuild)"), "\n\n", sep = "")
  
  # ---------------------------------------------------------------------------
  # 2. Validation & Fast Skip
  # ---------------------------------------------------------------------------
  if (!dir.exists(FACT_DB)) stop("Fact Table missing!")
  if (!file.exists(META_FILE)) stop("Metadata CSV missing!")
  
  # Check if final DB is already complete and valid (Fast Skip)
  if (dir.exists(OUTPUT_DB) && RESUME) {
    cat("Checking existing partitions and schema in db_stage2...\n")
    
    # A. Partition Count Check
    # We expect 128 partitions: partition_id=0 ... partition_id=127
    # And inside each, a part-XXX.parquet file.
    missing_parts <- 0
    for (i in 0:(N_PARTS - 1L)) {
      pfile <- file.path(OUTPUT_DB, paste0("partition_id=", i), sprintf("part-%03d.parquet", i))
      if (!file.exists(pfile)) {
        missing_parts <- missing_parts + 1
      }
    }
    
    # B. Schema Check (only if partitions look okay)
    schema_ok <- FALSE
    if (missing_parts == 0) {
      # Open the first partition to validate schema
      tryCatch({
        test_file <- file.path(OUTPUT_DB, "partition_id=0", "part-000.parquet")
        ds <- arrow::open_dataset(test_file)
        actual_cols <- names(ds)
        
        # Required columns for db_stage2
        required_cols <- c("Peptide", "N_Diseases", "N_Patients", "N_Antibodies", 
                           "is_healthy_background", "partition_id")
        
        if (all(required_cols %in% actual_cols)) {
          schema_ok <- TRUE
        } else {
          cat("[WARNING] Schema mismatch in existing db_stage2.\n")
          cat("Expected:", paste(required_cols, collapse=", "), "\n")
          cat("Found:   ", paste(actual_cols, collapse=", "), "\n")
        }
      }, error = function(e) {
        cat("[WARNING] Failed to read existing db_stage2 for validation:", conditionMessage(e), "\n")
      })
    }
    
    # C. Decision
    if (missing_parts == 0 && schema_ok) {
      cat("[SUCCESS] db_stage2 is complete and valid (128 partitions, correct schema).\n")
      cat("Skipping Step 4.\n")
      return(invisible(NULL))
    } else {
      cat(sprintf("[INFO] db_stage2 incomplete or invalid. Missing partitions: %d. Schema OK: %s. Resuming build...\n\n", 
                  missing_parts, schema_ok))
    }
  }
  
  # ---------------------------------------------------------------------------
  # 3. Initialize DuckDB (temporary build DB)
  # ---------------------------------------------------------------------------
  build_db <- file.path(P$work_dir, "build_dim.duckdb")
  if (file.exists(build_db)) unlink(build_db)
  
  con <- dbConnect(duckdb::duckdb(), dbdir = build_db)
  on.exit({
    dbDisconnect(con, shutdown = TRUE)
    if (file.exists(build_db)) unlink(build_db)
  }, add = TRUE)
  
  dbExecute(con, paste0("PRAGMA memory_limit='", DUCKDB_MEMORY_LIMIT, "'"))
  dbExecute(con, paste0("PRAGMA threads=", N_CORES))
  dbExecute(con, "PRAGMA preserve_insertion_order=false")
  
  # Temp spill directory (fast NVMe recommended)
  temp_spill <- file.path(P$work_dir, "duckdb_tmp")
  dir.create(temp_spill, showWarnings = FALSE, recursive = TRUE)
  dbExecute(con, paste0("PRAGMA temp_directory='", gsub("\\\\", "/", temp_spill), "'"))
  
  # ---------------------------------------------------------------------------
  # 4. Register Disease Facts input (view over stage1 parquet)
  # ---------------------------------------------------------------------------
  fact_glob <- paste0(gsub("\\\\", "/", FACT_DB), "/**/*.parquet")
  dbExecute(con, sprintf(
    "CREATE VIEW disease_facts AS SELECT * FROM read_parquet('%s')",
    fact_glob
  ))
  
  # ---------------------------------------------------------------------------
  # 5. Identify Healthy Background parquet files
  # ---------------------------------------------------------------------------
  cat("Identifying Healthy Background files...\n")
  meta <- fread(META_FILE, select = c("Filename", "Disease"))
  healthy_files <- meta[Disease == "None", Filename]
  healthy_paths <- file.path(DIGEST_DIR, gsub("\\.csv\\.gz$", ".parquet", healthy_files))
  healthy_paths <- healthy_paths[file.exists(healthy_paths)]
  
  if (length(healthy_paths) == 0) {
    cat("[WARNING] No Healthy background files found.\n")
    cat("          is_healthy_background will be FALSE for all peptides.\n\n")
  } else {
    cat("Found ", length(healthy_paths), " healthy background files.\n\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 6. Intermediate working directories (on disk)
  # ---------------------------------------------------------------------------
  # We store intermediates under work_dir so resume works across restarts.
  tmp_root <- gsub("\\\\", "/", file.path(P$work_dir, "dim_build_tmp_128parts"))
  
  disease_stats_dir   <- gsub("\\\\", "/", file.path(tmp_root, "disease_stats_part"))
  healthy_batches_dir <- gsub("\\\\", "/", file.path(tmp_root, "healthy_batches_part"))
  healthy_index_dir   <- gsub("\\\\", "/", file.path(tmp_root, "healthy_index_part"))
  
  if (!RESUME) {
    # Clean rebuild: delete tmp_root and final output before starting
    if (dir.exists(tmp_root)) unlink(tmp_root, recursive = TRUE)
    if (dir.exists(OUTPUT_DB)) unlink(OUTPUT_DB, recursive = TRUE)
  }
  
  dir.create(disease_stats_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(healthy_batches_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(healthy_index_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Helper: check if disease_stats seems completed
  has_disease_stats <- function() {
    any(grepl("partition_id=", list.files(disease_stats_dir, recursive = TRUE, full.names = TRUE)))
  }
  
  # Helper: check if final output seems completed (per partition)
  has_final_partition <- function(pid) {
    out_file <- file.path(OUTPUT_DB, paste0("partition_id=", pid), sprintf("part-%03d.parquet", pid))
    file.exists(out_file)
  }
  
  # ---------------------------------------------------------------------------
  # 7. Build disease_stats (or reuse if present)
  # ---------------------------------------------------------------------------
  cat("[1/4] Building disease_stats (distinct disease peptides + counts)...\n")
  
  if (RESUME && has_disease_stats()) {
    cat("disease_stats already exists at: ", disease_stats_dir, "  -> skipping.\n\n", sep = "")
  } else {
    t1 <- Sys.time()
    
    sql_disease_stats <- sprintf("
      COPY (
        SELECT
          Peptide,
          COUNT(DISTINCT Disease)  AS N_Diseases,
          COUNT(DISTINCT Patient)  AS N_Patients,
          COUNT(DISTINCT Antibody) AS N_Antibodies,
          (hash(Peptide) & %d)     AS partition_id
        FROM disease_facts
        GROUP BY Peptide
      )
      TO '%s'
      (FORMAT PARQUET,
       PARTITION_BY (partition_id),
       COMPRESSION ZSTD,
       ROW_GROUP_SIZE 1000000)
    ", PART_MASK, disease_stats_dir)
    
    tryCatch({
      dbExecute(con, sql_disease_stats)
      dbExecute(con, "CHECKPOINT")
    }, error = function(e) {
      stop(paste("Failed building disease_stats:", conditionMessage(e)))
    })
    
    t2 <- Sys.time()
    cat("disease_stats written to: ", disease_stats_dir, "\n", sep = "")
    cat("Time elapsed: ", round(difftime(t2, t1, units = "mins"), 2), " minutes.\n\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 8. Stream healthy background into per-batch directories (FIXED)
  # ---------------------------------------------------------------------------
  cat("[2/4] Streaming healthy background into partitioned batch parquet...\n")
  
  if (length(healthy_paths) == 0) {
    cat("No healthy files -> skipping batch streaming.\n\n")
  } else {
    
    # Plan batches by file count and total size
    file_info <- data.table(path = healthy_paths, size = file.size(healthy_paths))
    file_info[, batch_id := NA_integer_]
    
    curr_batch <- 1L
    curr_size  <- 0
    curr_cnt   <- 0L
    
    for (i in 1:nrow(file_info)) {
      f_size <- file_info$size[i]
      if ((curr_cnt >= MAX_BATCH_FILES) || ((curr_size + f_size) > MAX_BATCH_SIZE)) {
        curr_batch <- curr_batch + 1L
        curr_size  <- 0
        curr_cnt   <- 0L
      }
      file_info$batch_id[i] <- curr_batch
      curr_size <- curr_size + f_size
      curr_cnt  <- curr_cnt + 1L
    }
    
    n_batches <- max(file_info$batch_id)
    cat(sprintf("Planned %d batches.\n", n_batches))
    
    t1 <- Sys.time()
    
    for (b in 1:n_batches) {
      
      # Batch output directory (unique per batch) -> avoids "directory not empty" error.
      batch_out_dir <- file.path(healthy_batches_dir, sprintf("batch_%03d", b))
      batch_out_dir <- gsub("\\\\", "/", batch_out_dir)
      dir.create(batch_out_dir, showWarnings = FALSE, recursive = TRUE)
      
      # Resume: if this batch directory already has any parquet files, skip it.
      if (RESUME) {
        existing <- list.files(batch_out_dir, recursive = TRUE, pattern = "\\.parquet$", full.names = TRUE)
        if (length(existing) > 0) {
          cat(sprintf("  Healthy batch %d/%d already exists -> skipping.\n", b, n_batches))
          next
        }
      } else {
        # Clean rebuild: ensure batch output dir is empty
        if (dir.exists(batch_out_dir)) unlink(batch_out_dir, recursive = TRUE)
        dir.create(batch_out_dir, showWarnings = FALSE, recursive = TRUE)
      }
      
      batch_files   <- file_info[batch_id == b, path]
      batch_size_mb <- sum(file_info[batch_id == b, size]) / 1024^2
      
      cat(sprintf("  Healthy batch %d/%d (%d files, %.1f MB)...\n",
                  b, n_batches, length(batch_files), batch_size_mb))
      
      safe_paths <- paste0("'", gsub("\\\\", "/", batch_files), "'", collapse = ",")
      
      # DISTINCT within batch; write partitioned by (partition_id)
      # We also include batch_id as a column for traceability, but NOT as partition;
      # the batch directory already isolates batches, which avoids overwrite issues.
      sql_healthy_batch <- sprintf("
        COPY (
          SELECT DISTINCT
            Peptide,
            %d AS batch_id,
            (hash(Peptide) & %d) AS partition_id
          FROM read_parquet([%s])
        )
        TO '%s'
        (FORMAT PARQUET,
         PARTITION_BY (partition_id),
         COMPRESSION ZSTD)
      ", b, PART_MASK, safe_paths, batch_out_dir)
      
      tryCatch({
        dbExecute(con, sql_healthy_batch)
        dbExecute(con, "CHECKPOINT")
      }, error = function(e) {
        stop(paste("Healthy batch failed:", conditionMessage(e)))
      })
      
      gc()
    }
    
    t2 <- Sys.time()
    cat("Healthy batches written under: ", healthy_batches_dir, "\n", sep = "")
    cat("Time elapsed: ", round(difftime(t2, t1, units = "mins"), 2), " minutes.\n\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 9. Build a deduplicated healthy index per partition_id (resume-capable)
  # ---------------------------------------------------------------------------
  cat("[3/4] Building deduplicated healthy index per partition_id...\n")
  
  if (length(healthy_paths) == 0) {
    cat("No healthy files -> skipping healthy index build.\n\n")
  } else {
    
    t1 <- Sys.time()
    
    for (pid in 0:(N_PARTS - 1L)) {
      
      out_dir  <- gsub("\\\\", "/", file.path(healthy_index_dir, paste0("partition_id=", pid)))
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      
      out_file <- gsub("\\\\", "/", file.path(out_dir, sprintf("part-%03d.parquet", pid)))
      
      # Resume: if healthy index file already exists, skip it.
      if (RESUME && file.exists(out_file)) {
        if (pid %% 16 == 0) {
          cat(sprintf("  Healthy index partition %d/%d already exists -> skipping\n", pid, N_PARTS - 1L))
        }
        next
      }
      
      # Read all batch outputs for this partition_id:
      in_glob <- paste0(healthy_batches_dir, "/batch_*/partition_id=", pid, "/*.parquet")
      pid_files <- Sys.glob(in_glob)
      
      if (length(pid_files) == 0) {
        # Create empty parquet so downstream join doesn't fail.
        dbExecute(con, sprintf("
          COPY (
            SELECT CAST(NULL AS VARCHAR) AS Peptide
            WHERE FALSE
          )
          TO '%s'
          (FORMAT PARQUET, COMPRESSION ZSTD)
        ", out_file))
        next
      }
      
      # Global dedup within this partition only.
      sql_dedup_pid <- sprintf("
        COPY (
          SELECT DISTINCT Peptide
          FROM read_parquet('%s')
        )
        TO '%s'
        (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000)
      ", in_glob, out_file)
      
      tryCatch({
        dbExecute(con, sql_dedup_pid)
      }, error = function(e) {
        stop(paste("Healthy index dedup failed for partition", pid, ":", conditionMessage(e)))
      })
      
      if (pid %% 16 == 0) {
        cat(sprintf("  Healthy index deduped for partition %d/%d\n", pid, N_PARTS - 1L))
      }
    }
    
    t2 <- Sys.time()
    cat("Healthy index written to: ", healthy_index_dir, "\n", sep = "")
    cat("Time elapsed: ", round(difftime(t2, t1, units = "mins"), 2), " minutes.\n\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 10. Build final db_stage2 per partition with STRICT string join (resume-capable)
  # ---------------------------------------------------------------------------
  cat("[4/4] Writing final db_stage2 (strict join, 128 partitions)...\n")
  
  dir.create(OUTPUT_DB, showWarnings = FALSE, recursive = TRUE)
  
  t1 <- Sys.time()
  
  for (pid in 0:(N_PARTS - 1L)) {
    
    # Resume: skip if final output file exists
    if (RESUME && has_final_partition(pid)) {
      if (pid %% 16 == 0) {
        cat(sprintf("  Final partition %d/%d already exists -> skipping\n", pid, N_PARTS - 1L))
      }
      next
    }
    
    disease_glob_pid <- paste0(disease_stats_dir, "/partition_id=", pid, "/*.parquet")
    disease_files <- Sys.glob(disease_glob_pid)
    
    out_dir  <- gsub("\\\\", "/", file.path(OUTPUT_DB, paste0("partition_id=", pid)))
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    out_file <- gsub("\\\\", "/", file.path(out_dir, sprintf("part-%03d.parquet", pid)))
    
    if (length(disease_files) == 0) {
      # Empty disease partition -> write empty output partition file
      dbExecute(con, sprintf("
        COPY (
          SELECT
            CAST(NULL AS VARCHAR) AS Peptide,
            CAST(NULL AS BIGINT)  AS N_Diseases,
            CAST(NULL AS BIGINT)  AS N_Patients,
            CAST(NULL AS BIGINT)  AS N_Antibodies,
            CAST(NULL AS BOOLEAN) AS is_healthy_background,
            %d                    AS partition_id
          WHERE FALSE
        )
        TO '%s'
        (FORMAT PARQUET, COMPRESSION ZSTD)
      ", pid, out_file))
      next
    }
    
    healthy_file_pid <- gsub("\\\\", "/", file.path(
      healthy_index_dir,
      paste0("partition_id=", pid),
      sprintf("part-%03d.parquet", pid)
    ))
    
    # Strict string join (exact correctness)
    sql_final_pid <- sprintf("
      COPY (
        SELECT
          d.Peptide,
          d.N_Diseases,
          d.N_Patients,
          d.N_Antibodies,
          CASE WHEN h.Peptide IS NOT NULL THEN TRUE ELSE FALSE END AS is_healthy_background,
          %d AS partition_id
        FROM read_parquet('%s') d
        LEFT JOIN read_parquet('%s') h
          ON d.Peptide = h.Peptide
      )
      TO '%s'
      (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000)
    ", pid, disease_glob_pid, healthy_file_pid, out_file)
    
    tryCatch({
      dbExecute(con, sql_final_pid)
    }, error = function(e) {
      stop(paste("Final join failed for partition", pid, ":", conditionMessage(e)))
    })
    
    if (pid %% 16 == 0) {
      cat(sprintf("  Wrote db_stage2 partition %d/%d\n", pid, N_PARTS - 1L))
    }
  }
  
  t2 <- Sys.time()
  cat("\n[SUCCESS] db_stage2 built at: ", OUTPUT_DB, "\n", sep = "")
  cat("Final phase time elapsed: ", round(difftime(t2, t1, units="mins"), 2), " minutes.\n", sep = "")
  cat("Intermediate build data kept at: ", tmp_root, "\n", sep = "")
  cat("Resume is enabled; rerun safely to continue after interruptions.\n")
}

main()