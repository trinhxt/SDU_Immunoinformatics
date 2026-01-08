# Part5_CDR3.R  (UPDATED)
# Goal:
#   For each Peptide, answer: "Is the peptide contained in cdr3_aa?" as:
#     Yes = 1, No = 0
#
# Interpretation (row-consistent):
#   We check containment within the *same row* (Peptide vs that row's cdr3_aa),
#   then aggregate to Peptide-level:
#     peptide_in_cdr3 = 1 if the peptide is found in ANY observed cdr3_aa context
#                     = 0 otherwise
#
# Scalable plan:
#   1) DISTINCT (Peptide, cdr3_aa) pairs (dedupe early) + bucket by Peptide
#   2) Compute per-pair flag via instr(cdr3_aa, Peptide) > 0
#   3) Aggregate to one row per Peptide: MAX(flag) -> 1/0

library(DBI)
library(duckdb)

FORCE_REBUILD <- FALSE  # set TRUE to rebuild even if outputs exist

# ---- Helpers: always use forward slashes + absolute paths for DuckDB ----
norm  <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
pjoin <- function(...) norm(paste(..., sep = "/"))

assert_no_backslash <- function(x, name) {
  if (grepl("\\\\", x)) stop("Backslash found in ", name, ":\n", x)
}

part_done <- function(out_dir) {
  if (!dir.exists(out_dir)) return(FALSE)
  length(list.files(out_dir, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)) > 0
}

# -----------------------------
# INPUT (your big parquet DB)
# -----------------------------
final_root   <- norm("D:/OAS/OAS_human_heavychain_disease_tryptic/parquet_db_partitioned/final")
parquet_glob <- paste0(final_root, "/**/*.parquet")

# -----------------------------
# OUTPUT WORK AREA
# -----------------------------
work_root <- norm("OAS_human_heavychain_disease_tryptic/CDR3_check")  # <-- CHANGE if needed
dir.create(work_root, showWarnings = FALSE, recursive = TRUE)

duckdb_file <- pjoin(work_root, "build.duckdb")
duckdb_temp <- pjoin(work_root, "duckdb_temp")
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

N_BUCKETS <- 256L

# Output dirs
pair_dir    <- pjoin(work_root, "base_distinct", "peptide_cdr3_pair")  # distinct pairs
answer_dir  <- pjoin(work_root, "answer", "peptide")                   # peptide -> 1/0

dir.create(pair_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(answer_dir, showWarnings = FALSE, recursive = TRUE)

# Sanity checks
assert_no_backslash(parquet_glob, "parquet_glob")
assert_no_backslash(work_root, "work_root")
assert_no_backslash(duckdb_file, "duckdb_file")
assert_no_backslash(duckdb_temp, "duckdb_temp")
assert_no_backslash(pair_dir, "pair_dir")
assert_no_backslash(answer_dir, "answer_dir")

# -----------------------------
# DUCKDB
# -----------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

dbExecute(con, "PRAGMA threads=6;")
dbExecute(con, "PRAGMA memory_limit='44GB';")
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
dbExecute(con, "PRAGMA enable_progress_bar=true;")

# ============================================================
# PART 1/2: DISTINCT (Peptide, cdr3_aa) pairs
# ============================================================
if (!FORCE_REBUILD && part_done(pair_dir)) {
  cat("\n[1/2] Distinct (Peptide, cdr3_aa) pairs already exist -> SKIP\n")
} else {
  cat("\n[1/2] Writing distinct (Peptide, cdr3_aa) pairs...\n")
  
  sql_pairs <- sprintf("
  COPY (
    WITH base AS (
      SELECT Peptide, cdr3_aa
      FROM read_parquet('%s')
      WHERE Peptide IS NOT NULL AND trim(Peptide) <> ''
        AND cdr3_aa IS NOT NULL AND trim(cdr3_aa) <> ''
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
#   If peptide appears in ANY cdr3_aa context -> 1 else 0
# ============================================================
if (!FORCE_REBUILD && part_done(answer_dir)) {
  cat("\n[2/2] Peptide CDR3 answer table already exists -> SKIP\n")
} else {
  cat("\n[2/2] Writing peptide CDR3 answer table (1/0)...\n")
  
  pair_glob <- paste0(pair_dir, "/**/*.parquet")
  
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

cat("\nDONE.\n")
cat("Distinct pairs:", pair_dir, "\n")
cat("Peptide answer (1/0):", answer_dir, "\n")

# Optional quick peek:
# answer_glob <- paste0(answer_dir, "/**/*.parquet")
# print(dbGetQuery(con, sprintf("
#   SELECT peptide_in_cdr3, COUNT(*) AS n_peptides
#   FROM read_parquet('%s')
#   GROUP BY 1
#   ORDER BY 1 DESC;", answer_glob)))
