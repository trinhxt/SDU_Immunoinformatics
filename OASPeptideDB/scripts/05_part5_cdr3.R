################################################################################
## scripts/05_part5_cdr3.R
##
## Part5: Peptide-in-CDR3 flag (CDR3-SPANNING by overlap > 3 AA)
##
## Update requested:
## - First dedupe DISTINCT (Antibody, Peptide, cdr3_aa)
## - Then mark peptide_in_cdr3 = 1 if Peptide overlaps the CDR3 region by >3 AA
##   using positional overlap in the Antibody sequence:
##
##   (Peptide_Start <= CDR3_End) AND (Peptide_End >= CDR3_Start)
##   overlap_len = min(Peptide_End, CDR3_End) - max(Peptide_Start, CDR3_Start) + 1
##   peptide_in_cdr3 = 1 if overlap_len > 3 else 0
##
## Assumptions:
## - "Antibody" = Antibody (full AA sequence) in Part2 parquet.
## - Peptide_Start and CDR3_Start are located by instr(Antibody, substring).
##   (instr() returns first match; if there are repeats, this uses the first.)
##
## Uses config/paths.yml (via helpers_paths.R):
##   - P$final_part2_dir : final Parquet dataset from Part2
##   - P$part5_root      : output root for Part5
##
## Outputs (under P$part5_root):
##   - base_distinct/antibody_peptide_cdr3/   (bucket-partitioned distinct triples)
##   - answer/peptide/                        (bucket-partitioned peptide_in_cdr3)
##
## Resume/Skip
## ----------
## - Each output folder is considered "done" if it contains any parquet files.
## - Set FORCE_REBUILD=TRUE to rebuild outputs even if they exist.
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
  
  # Updated: distinct triples (Antibody, Peptide, cdr3_aa)
  triple_dir <- pjoin(work_root, "base_distinct", "antibody_peptide_cdr3")
  answer_dir <- pjoin(work_root, "answer", "peptide")
  
  dir.create(triple_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(answer_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Sanity checks
  assert_no_backslash(work_root,  "work_root")
  assert_no_backslash(duckdb_file,"duckdb_file")
  assert_no_backslash(duckdb_temp,"duckdb_temp")
  assert_no_backslash(triple_dir, "triple_dir")
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
  # PART 1/2: DISTINCT (Antibody, Peptide, cdr3_aa) triples
  # ============================================================
  if (!FORCE_REBUILD && part_done(triple_dir)) {
    cat("[1/2] Distinct (Antibody, Peptide, cdr3_aa) already exist -> SKIP\n")
  } else {
    cat("[1/2] Writing distinct (Antibody, Peptide, cdr3_aa) triples...\n")
    
    sql_triples <- sprintf("
      COPY (
        WITH base AS (
          SELECT
            Antibody AS Antibody,
            Peptide,
            cdr3_aa
          FROM read_parquet('%s')
          WHERE Antibody IS NOT NULL AND trim(Antibody) <> ''
            AND Peptide              IS NOT NULL AND trim(Peptide)              <> ''
            AND cdr3_aa              IS NOT NULL AND trim(cdr3_aa)              <> ''
        )
        SELECT DISTINCT
          Antibody,
          Peptide,
          cdr3_aa,
          (hash(Peptide) %% %d)::INTEGER AS bucket
        FROM base
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", parquet_glob, N_BUCKETS, triple_dir)
    
    assert_no_backslash(sql_triples, "sql_triples")
    dbExecute(con, sql_triples)
    
    if (!part_done(triple_dir)) stop("Part 1 failed: no parquet written in ", triple_dir)
  }
  
  # ============================================================
  # PART 2/2: Peptide-level answer (1/0) based on overlap > 3 AA
  # ============================================================
  if (!FORCE_REBUILD && part_done(answer_dir)) {
    cat("[2/2] Peptide CDR3 answer table already exists -> SKIP\n")
  } else {
    cat("[2/2] Writing peptide CDR3 answer table (peptide_in_cdr3 by overlap > 3 AA)...\n")
    
    triple_glob <- paste0(triple_dir, "/**/*.parquet")
    assert_no_backslash(triple_glob, "triple_glob")
    
    sql_answer <- sprintf("
      COPY (
        WITH t AS (
          SELECT
            Antibody,
            Peptide,
            cdr3_aa,
            (hash(Peptide) %% %d)::INTEGER AS bucket,
            instr(Antibody, Peptide) AS pep_start,
            instr(Antibody, cdr3_aa) AS cdr3_start,
            length(Peptide) AS pep_len,
            length(cdr3_aa) AS cdr3_len
          FROM read_parquet('%s')
        ),
        pos AS (
          SELECT
            Peptide,
            bucket,
            pep_start,
            (pep_start  + pep_len  - 1) AS pep_end,
            cdr3_start,
            (cdr3_start + cdr3_len - 1) AS cdr3_end
          FROM t
          WHERE pep_start  > 0
            AND cdr3_start > 0
        ),
        scored AS (
          SELECT
            Peptide,
            bucket,
            CASE
              WHEN (pep_start <= cdr3_end) AND (pep_end >= cdr3_start)
              THEN
                -- overlap_len = min(pep_end,cdr3_end) - max(pep_start,cdr3_start) + 1
                CASE
                  WHEN (LEAST(pep_end, cdr3_end) - GREATEST(pep_start, cdr3_start) + 1) > 3
                  THEN 1 ELSE 0
                END
              ELSE 0
            END AS span_flag
          FROM pos
        )
        SELECT
          Peptide,
          bucket,
          MAX(span_flag) AS peptide_in_cdr3
        FROM scored
        GROUP BY Peptide, bucket
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, triple_glob, answer_dir)
    
    assert_no_backslash(sql_answer, "sql_answer")
    dbExecute(con, sql_answer)
    
    if (!part_done(answer_dir)) stop("Part 2 failed: no parquet written in ", answer_dir)
  }
  
  cat("\nDONE (Part5).\n")
  cat("Distinct triples:\n  ", triple_dir, "\n", sep = "")
  cat("Peptide answer (1/0) peptide_in_cdr3:\n  ", answer_dir, "\n", sep = "")
  
  invisible(TRUE)
}

main()
