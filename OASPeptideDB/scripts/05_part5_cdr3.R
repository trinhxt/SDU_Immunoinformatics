################################################################################
## scripts/05_part5_cdr3.R
##
## Part5: Peptide-in-CDR3 flag
## Goal:
##   For each Peptide, compute peptide_in_cdr3:
##     1 if peptide is contained in cdr3_aa in ANY observed row context
##     0 otherwise
##
## Scalable plan:
##   1) DISTINCT (Peptide, cdr3_aa) pairs (dedupe early) + bucket by Peptide
##   2) Per-pair flag via instr(cdr3_aa, Peptide) > 0
##   3) Aggregate to one row per Peptide: MAX(flag) -> 1/0
##
## Uses config/paths.yml (via helpers_paths.R):
##   - P$final_part2_dir : final Parquet dataset from Part2
##   - P$part5_root      : output root for Part5
##
## Outputs (under P$part5_root):
##   - base_distinct/peptide_cdr3_pair/  (bucket-partitioned distinct pairs)
##   - answer/peptide/                   (bucket-partitioned peptide_in_cdr3)
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
  # OUTPUT (Part5 root)
  # ---------------------------------------------------------------------------
  work_root <- norm(P$part5_root)
  dir.create(work_root, showWarnings = FALSE, recursive = TRUE)
  
  duckdb_file <- pjoin(work_root, "build.duckdb")
  duckdb_temp <- pjoin(work_root, "duckdb_temp")
  dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)
  
  pair_dir   <- pjoin(work_root, "base_distinct", "peptide_cdr3_pair")
  answer_dir <- pjoin(work_root, "answer", "peptide")
  
  dir.create(pair_dir,   showWarnings = FALSE, recursive = TRUE)
  dir.create(answer_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Sanity checks
  assert_no_backslash(work_root, "work_root")
  assert_no_backslash(duckdb_file, "duckdb_file")
  assert_no_backslash(duckdb_temp, "duckdb_temp")
  assert_no_backslash(pair_dir, "pair_dir")
  assert_no_backslash(answer_dir, "answer_dir")
  
  cat("Part5 input Parquet dataset:\n  ", final_root, "\n", sep = "")
  cat("Part5 output root:\n  ", work_root, "\n\n", sep = "")
  
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
  # PART 1/2: DISTINCT (Peptide, cdr3_aa) pairs
  # ============================================================
  if (!FORCE_REBUILD && part_done(pair_dir)) {
    cat("[1/2] Distinct (Peptide, cdr3_aa) pairs already exist -> SKIP\n")
  } else {
    cat("[1/2] Writing distinct (Peptide, cdr3_aa) pairs...\n")
    
    sql_pairs <- sprintf("
      COPY (
        WITH base AS (
          SELECT Peptide, cdr3_aa
          FROM read_parquet('%s')
          WHERE Peptide  IS NOT NULL AND trim(Peptide)  <> ''
            AND cdr3_aa  IS NOT NULL AND trim(cdr3_aa)  <> ''
        )
        SELECT DISTINCT
          Peptide,
          cdr3_aa,
          (hash(Peptide) %% %d)::INTEGER AS bucket
        FROM base
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", parquet_glob, N_BUCKETS, pair_dir)
    
    assert_no_backslash(sql_pairs, "sql_pairs")
    dbExecute(con, sql_pairs)
    
    if (!part_done(pair_dir)) stop("Part 1 failed: no parquet written in ", pair_dir)
  }
  
  # ============================================================
  # PART 2/2: Peptide-level answer (1/0)
  # ============================================================
  if (!FORCE_REBUILD && part_done(answer_dir)) {
    cat("[2/2] Peptide CDR3 answer table already exists -> SKIP\n")
  } else {
    cat("[2/2] Writing peptide CDR3 answer table (peptide_in_cdr3)...\n")
    
    pair_glob <- paste0(pair_dir, "/**/*.parquet")
    assert_no_backslash(pair_glob, "pair_glob")
    
    sql_answer <- sprintf("
      COPY (
        SELECT
          Peptide,
          (hash(Peptide) %% %d)::INTEGER AS bucket,
          MAX(CASE WHEN instr(cdr3_aa, Peptide) > 0 THEN 1 ELSE 0 END) AS peptide_in_cdr3
        FROM read_parquet('%s')
        GROUP BY Peptide, bucket
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, pair_glob, answer_dir)
    
    assert_no_backslash(sql_answer, "sql_answer")
    dbExecute(con, sql_answer)
    
    if (!part_done(answer_dir)) stop("Part 2 failed: no parquet written in ", answer_dir)
  }
  
  cat("\nDONE (Part5).\n")
  cat("Distinct pairs:\n  ", pair_dir, "\n", sep = "")
  cat("Peptide answer (1/0):\n  ", answer_dir, "\n", sep = "")
  
  invisible(TRUE)
}

main()
