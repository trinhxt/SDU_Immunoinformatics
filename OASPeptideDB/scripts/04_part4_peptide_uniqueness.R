################################################################################
## scripts/04_part4_peptide_uniqueness.R
##
## Part4: Peptide uniqueness
## Goal: For each Peptide, count how many DISTINCT Antibodies it appears in.
##
## Uses config/paths.yml (via helpers_paths.R):
##   - P$final_part2_dir : final Parquet dataset from Part2
##   - P$part4_root      : output root for Part4
##
## Outputs (under P$part4_root):
##   - base_distinct/peptide_antibody/  (bucket-partitioned distinct map)
##   - counts/peptide/                  (bucket-partitioned counts)
##
## Resume/Skip
## ----------
## - Each output folder is considered "done" if it contains any parquet files.
## - Set FORCE_REBUILD=TRUE to rebuild outputs even if they exist.
##
## Notes
## -----
## - Safe for source(): wrapped in main(); no quit()
## - Windows path safe via helpers (forward slashes for DuckDB)
################################################################################

# ==============================================================================
# Load config + helpers
# ==============================================================================
source("R/helpers_paths.R")
cfg <- load_config()
P   <- get_paths(cfg)

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # Settings
  # ---------------------------------------------------------------------------
  FORCE_REBUILD  <- FALSE
  N_BUCKETS      <- 256L
  
  DUCKDB_THREADS <- 6L
  DUCKDB_MEM     <- "44GB"
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  part_done <- function(out_dir) {
    if (!dir.exists(out_dir)) return(FALSE)
    length(list.files(out_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)) > 0
  }
  
  # ---------------------------------------------------------------------------
  # INPUT (Part2 final dataset)
  # ---------------------------------------------------------------------------
  final_root   <- norm(P$final_part2_dir)
  parquet_glob <- paste0(final_root, "/**/*.parquet")
  
  if (!dir.exists(final_root)) {
    stop("Part2 final dataset folder not found: ", final_root)
  }
  assert_no_backslash(parquet_glob, "parquet_glob")
  
  # ---------------------------------------------------------------------------
  # OUTPUT (Part4 root)
  # ---------------------------------------------------------------------------
  work_root <- norm(P$part4_root)
  dir.create(work_root, showWarnings = FALSE, recursive = TRUE)
  
  duckdb_file <- pjoin(work_root, "build.duckdb")
  duckdb_temp <- pjoin(work_root, "duckdb_temp")
  dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)
  
  map_dir   <- pjoin(work_root, "base_distinct", "peptide_antibody")
  count_dir <- pjoin(work_root, "counts", "peptide")
  
  dir.create(map_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(count_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Sanity checks
  assert_no_backslash(work_root, "work_root")
  assert_no_backslash(duckdb_file, "duckdb_file")
  assert_no_backslash(duckdb_temp, "duckdb_temp")
  assert_no_backslash(map_dir, "map_dir")
  assert_no_backslash(count_dir, "count_dir")
  
  cat("Part4 input Parquet dataset:\n  ", final_root, "\n", sep = "")
  cat("Part4 output root:\n  ", work_root, "\n\n", sep = "")
  
  # ---------------------------------------------------------------------------
  # DUCKDB
  # ---------------------------------------------------------------------------
  con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  
  dbExecute(con, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
  dbExecute(con, "PRAGMA enable_progress_bar=true;")
  
  # ============================================================
  # PART 1/2: Distinct (Peptide, Antibody) mapping
  # ============================================================
  if (!FORCE_REBUILD && part_done(map_dir)) {
    cat("[1/2] Distinct peptide-antibody map already exists -> SKIP\n")
  } else {
    cat("[1/2] Writing distinct peptide-antibody map...\n")
    
    sql_map <- sprintf("
      COPY (
        WITH base AS (
          SELECT Peptide, Antibody
          FROM read_parquet('%s')
          WHERE Peptide  IS NOT NULL AND trim(Peptide)  <> ''
            AND Antibody IS NOT NULL AND trim(Antibody) <> ''
        )
        SELECT DISTINCT
          Peptide,
          Antibody,
          (hash(Peptide) %% %d)::INTEGER AS bucket
        FROM base
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", parquet_glob, N_BUCKETS, map_dir)
    
    assert_no_backslash(sql_map, "sql_map")
    dbExecute(con, sql_map)
    
    if (!part_done(map_dir)) stop("Part 1 failed: no parquet written in ", map_dir)
  }
  
  # ============================================================
  # PART 2/2: Count distinct antibodies per peptide
  # ============================================================
  if (!FORCE_REBUILD && part_done(count_dir)) {
    cat("[2/2] Peptide uniqueness counts already exist -> SKIP\n")
  } else {
    cat("[2/2] Writing peptide uniqueness counts (n_distinct_antibodies)...\n")
    
    map_glob <- paste0(map_dir, "/**/*.parquet")
    assert_no_backslash(map_glob, "map_glob")
    
    sql_counts <- sprintf("
      COPY (
        SELECT
          Peptide,
          (hash(Peptide) %% %d)::INTEGER AS bucket,
          COUNT(DISTINCT Antibody) AS n_distinct_antibodies
        FROM read_parquet('%s')
        GROUP BY Peptide, bucket
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, map_glob, count_dir)
    
    assert_no_backslash(sql_counts, "sql_counts")
    dbExecute(con, sql_counts)
    
    if (!part_done(count_dir)) stop("Part 2 failed: no parquet written in ", count_dir)
  }
  
  cat("\nDONE (Part4).\n")
  cat("Distinct peptide-antibody map:\n  ", map_dir, "\n", sep = "")
  cat("Peptide uniqueness counts:\n  ", count_dir, "\n", sep = "")
  
  invisible(TRUE)
}

main()
