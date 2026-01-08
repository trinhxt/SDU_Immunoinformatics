# Part4_Peptide_uniqueness.R
# Goal: For each Peptide, count how many DISTINCT Antibodies it appears in.
# Output: partitioned parquet (by bucket) with columns:
#   Peptide, bucket, n_distinct_antibodies

library(DBI)
library(duckdb)

# ============================================================
# Resume/skip logic:
# - This script writes into its own output folders.
# - If output folder already has parquet files, we SKIP on re-run.
# - Force rebuild by deleting output folders or setting FORCE_REBUILD=TRUE.
# ============================================================

FORCE_REBUILD <- FALSE  # set TRUE to rebuild even if outputs exist

# helper: absolute + forward slashes (Windows-safe for DuckDB)
norm <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

# helper: join without file.path() (prevents backslashes on Windows)
pjoin <- function(...) norm(paste(..., sep = "/"))

# Optional: catch backslashes before running
assert_no_backslash <- function(x, name) {
  if (grepl("\\\\", x)) stop("Backslash found in ", name, ":\n", x)
}

# Check whether a part is already done (folder exists and has parquet files)
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
# OUTPUT WORK AREA (BIG DRIVE)
# -----------------------------
work_root <- norm("OAS_human_heavychain_disease_tryptic/Peptide_uniqueness")  # <-- CHANGE if needed
dir.create(work_root, showWarnings = FALSE, recursive = TRUE)

duckdb_file <- pjoin(work_root, "build.duckdb")
duckdb_temp <- pjoin(work_root, "duckdb_temp")
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

# Buckets (match your other scripts if you want)
N_BUCKETS <- 256L

# Output dirs
map_dir   <- pjoin(work_root, "base_distinct", "peptide_antibody")  # distinct (Peptide, Antibody)
count_dir <- pjoin(work_root, "counts", "peptide")                  # peptide -> #antibodies

dir.create(map_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(count_dir, showWarnings = FALSE, recursive = TRUE)

# Sanity checks (must be all forward slashes)
assert_no_backslash(parquet_glob, "parquet_glob")
assert_no_backslash(work_root, "work_root")
assert_no_backslash(duckdb_file, "duckdb_file")
assert_no_backslash(duckdb_temp, "duckdb_temp")
assert_no_backslash(map_dir, "map_dir")
assert_no_backslash(count_dir, "count_dir")

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
# PART 1: Distinct Peptide-Antibody mapping (dedupe early)
# ============================================================
if (!FORCE_REBUILD && part_done(map_dir)) {
  cat("\n[1/2] Distinct peptide-antibody map already exists -> SKIP\n")
} else {
  cat("\n[1/2] Writing distinct peptide-antibody map...\n")
  
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
# PART 2: Count distinct antibodies per peptide
# ============================================================
if (!FORCE_REBUILD && part_done(count_dir)) {
  cat("\n[2/2] Peptide uniqueness table already exists -> SKIP\n")
} else {
  cat("\n[2/2] Writing peptide uniqueness table (n_distinct_antibodies)...\n")
  
  map_glob <- paste0(map_dir, "/**/*.parquet")
  
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

cat("\nDONE.\n")
cat("Distinct peptide-antibody map:", map_dir, "\n")
cat("Peptide uniqueness counts:", count_dir, "\n")

# Optional: quick example query (comment out if not needed)
# counts_glob <- paste0(count_dir, "/**/*.parquet")
# head(dbGetQuery(con, sprintf("SELECT * FROM read_parquet('%s') LIMIT 10;", counts_glob)))
