# Script: Part3_Disease_presence.R
library(DBI)
library(duckdb)

# ============================================================
# Resume/skip logic:
# - Each of the 4 parts writes into its own output folder.
# - If that folder already contains parquet files, we assume that part succeeded
#   and we SKIP it on re-run.
# - You can force a rebuild by deleting that folder (or set FORCE_REBUILD=TRUE).
# ============================================================

FORCE_REBUILD <- FALSE  # set TRUE to rebuild all parts even if outputs exist

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
# INPUT
# -----------------------------
final_root   <- norm("D:/OAS/OAS_human_heavychain_disease_tryptic/parquet_db_partitioned/final")
parquet_glob <- paste0(final_root, "/**/*.parquet")

# -----------------------------
# OUTPUT WORK AREA (BIG DRIVE)
# -----------------------------
work_root <- norm("OAS_human_heavychain_disease_tryptic/Disease_presence")  # <-- CHANGE if needed
dir.create(work_root, showWarnings = FALSE, recursive = TRUE)

duckdb_file <- pjoin(work_root, "build.duckdb")
duckdb_temp <- pjoin(work_root, "duckdb_temp")
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

# Buckets
N_BUCKETS <- 256L

# Output dirs
ab_base_dir <- pjoin(work_root, "base_distinct", "antibody")
pp_base_dir <- pjoin(work_root, "base_distinct", "peptide")

ab_pres_dir <- pjoin(work_root, "presence", "antibody")
pp_pres_dir <- pjoin(work_root, "presence", "peptide")

dir.create(ab_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pp_base_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ab_pres_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pp_pres_dir, showWarnings = FALSE, recursive = TRUE)

# Sanity checks (must be all forward slashes)
assert_no_backslash(parquet_glob, "parquet_glob")
assert_no_backslash(work_root, "work_root")
assert_no_backslash(duckdb_file, "duckdb_file")
assert_no_backslash(duckdb_temp, "duckdb_temp")
assert_no_backslash(ab_base_dir, "ab_base_dir")
assert_no_backslash(pp_base_dir, "pp_base_dir")
assert_no_backslash(ab_pres_dir, "ab_pres_dir")
assert_no_backslash(pp_pres_dir, "pp_pres_dir")

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
# PART 1: Antibody DISTINCT base
# ============================================================
if (!FORCE_REBUILD && part_done(ab_base_dir)) {
  cat("\n[1/4] Antibody DISTINCT base already exists -> SKIP\n")
} else {
  cat("\n[1/4] Writing antibody DISTINCT base...\n")
  
  sql_ab_base <- sprintf("
  COPY (
    WITH base AS (
      SELECT Disease, Isotype, BSource, BType, Antibody
      FROM read_parquet('%s')
      WHERE Disease IS NOT NULL AND trim(Disease) <> '' AND Disease <> 'None'
        AND Antibody IS NOT NULL AND trim(Antibody) <> ''
    )
    SELECT DISTINCT
      Disease, Isotype, BSource, BType, Antibody,
      (hash(Antibody) %% %d)::INTEGER AS bucket
    FROM base
  ) TO '%s'
  (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
  ", parquet_glob, N_BUCKETS, ab_base_dir)
  
  assert_no_backslash(sql_ab_base, "sql_ab_base")
  dbExecute(con, sql_ab_base)
  
  if (!part_done(ab_base_dir)) stop("Part 1 failed: no parquet written in ", ab_base_dir)
}

# ============================================================
# PART 2: Peptide DISTINCT base
# ============================================================
if (!FORCE_REBUILD && part_done(pp_base_dir)) {
  cat("\n[2/4] Peptide DISTINCT base already exists -> SKIP\n")
} else {
  cat("\n[2/4] Writing peptide DISTINCT base...\n")
  
  sql_pp_base <- sprintf("
  COPY (
    WITH base AS (
      SELECT Disease, Isotype, BSource, BType, Peptide
      FROM read_parquet('%s')
      WHERE Disease IS NOT NULL AND trim(Disease) <> '' AND Disease <> 'None'
        AND Peptide IS NOT NULL AND trim(Peptide) <> ''
    )
    SELECT DISTINCT
      Disease, Isotype, BSource, BType, Peptide,
      (hash(Peptide) %% %d)::INTEGER AS bucket
    FROM base
  ) TO '%s'
  (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
  ", parquet_glob, N_BUCKETS, pp_base_dir)
  
  assert_no_backslash(sql_pp_base, "sql_pp_base")
  dbExecute(con, sql_pp_base)
  
  if (!part_done(pp_base_dir)) stop("Part 2 failed: no parquet written in ", pp_base_dir)
}

# ============================================================
# PART 3: Antibody presence table
# ============================================================
if (!FORCE_REBUILD && part_done(ab_pres_dir)) {
  cat("\n[3/4] Antibody presence table already exists -> SKIP\n")
} else {
  cat("\n[3/4] Writing antibody presence table...\n")
  
  ab_base_glob <- paste0(ab_base_dir, "/**/*.parquet")
  
  sql_ab_pres <- sprintf("
  COPY (
    SELECT
      Antibody,
      (hash(Antibody) %% %d)::INTEGER AS bucket,
      COUNT(DISTINCT Disease) AS Disease_presence_antibody
    FROM read_parquet('%s')
    GROUP BY Antibody, bucket
  ) TO '%s'
  (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
  ", N_BUCKETS, ab_base_glob, ab_pres_dir)
  
  assert_no_backslash(sql_ab_pres, "sql_ab_pres")
  dbExecute(con, sql_ab_pres)
  
  if (!part_done(ab_pres_dir)) stop("Part 3 failed: no parquet written in ", ab_pres_dir)
}

# ============================================================
# PART 4: Peptide presence table
# ============================================================
if (!FORCE_REBUILD && part_done(pp_pres_dir)) {
  cat("\n[4/4] Peptide presence table already exists -> SKIP\n")
} else {
  cat("\n[4/4] Writing peptide presence table...\n")
  
  pp_base_glob <- paste0(pp_base_dir, "/**/*.parquet")
  
  sql_pp_pres <- sprintf("
  COPY (
    SELECT
      Peptide,
      (hash(Peptide) %% %d)::INTEGER AS bucket,
      COUNT(DISTINCT Disease) AS Disease_presence_peptide
    FROM read_parquet('%s')
    GROUP BY Peptide, bucket
  ) TO '%s'
  (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
  ", N_BUCKETS, pp_base_glob, pp_pres_dir)
  
  assert_no_backslash(sql_pp_pres, "sql_pp_pres")
  dbExecute(con, sql_pp_pres)
  
  if (!part_done(pp_pres_dir)) stop("Part 4 failed: no parquet written in ", pp_pres_dir)
}

cat("\nDONE.\n")
cat("Antibody base:", ab_base_dir, "\n")
cat("Antibody presence:", ab_pres_dir, "\n")
cat("Peptide base:", pp_base_dir, "\n")
cat("Peptide presence:", pp_pres_dir, "\n")




################################################################################
# Part 5 # Get stats tables for antibody and peptides (Windows path-safe)
library(DBI)
library(duckdb)
library(data.table)

# ---- Helpers: always use forward slashes + absolute paths for DuckDB ----
norm  <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
pjoin <- function(...) norm(paste(..., sep = "/"))

assert_no_backslash <- function(x, name) {
  if (grepl("\\\\", x)) stop("Backslash found in ", name, ":\n", x)
}

# -----------------------------
# PATHS (MUST match build script)
# -----------------------------
work_root  <- norm("OAS_human_heavychain_disease_tryptic/Disease_presence")   # same as above

duckdb_file <- pjoin(work_root, "query.duckdb")
duckdb_temp <- pjoin(work_root, "duckdb_temp_query")
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

ab_base_dir <- pjoin(work_root, "base_distinct", "antibody")
pp_base_dir <- pjoin(work_root, "base_distinct", "peptide")
ab_pres_dir <- pjoin(work_root, "presence", "antibody")
pp_pres_dir <- pjoin(work_root, "presence", "peptide")

# DuckDB globs (always explicit forward slashes)
ab_base_glob <- paste0(ab_base_dir, "/**/*.parquet")
pp_base_glob <- paste0(pp_base_dir, "/**/*.parquet")
ab_pres_glob <- paste0(ab_pres_dir, "/**/*.parquet")
pp_pres_glob <- paste0(pp_pres_dir, "/**/*.parquet")

# Sanity checks
assert_no_backslash(work_root, "work_root")
assert_no_backslash(duckdb_file, "duckdb_file")
assert_no_backslash(duckdb_temp, "duckdb_temp")
assert_no_backslash(ab_base_glob, "ab_base_glob")
assert_no_backslash(pp_base_glob, "pp_base_glob")
assert_no_backslash(ab_pres_glob, "ab_pres_glob")
assert_no_backslash(pp_pres_glob, "pp_pres_glob")

# -----------------------------
# DUCKDB
# -----------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)

dbExecute(con, "PRAGMA threads=6;")
dbExecute(con, "PRAGMA memory_limit='44GB';")
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
dbExecute(con, "PRAGMA enable_progress_bar=true;")

# -----------------------------
# Antibody sankey stats (fast)
# -----------------------------
sql_ab_stats <- sprintf("
SELECT
  'Human' AS Species,
  p.Disease_presence_antibody,
  b.Disease,
  b.Isotype,
  b.BSource,
  b.BType,
  COUNT(DISTINCT b.Antibody) AS n_unique_antibodies
FROM read_parquet('%s') b
JOIN read_parquet('%s') p
  ON b.Antibody = p.Antibody
GROUP BY 1,2,3,4,5,6
ORDER BY 3,4,5,6,2;
", ab_base_glob, ab_pres_glob)

assert_no_backslash(sql_ab_stats, "sql_ab_stats")
tbl_unique_antibodies <- as.data.table(dbGetQuery(con, sql_ab_stats))

# -----------------------------
# Peptide sankey stats (fast)
# -----------------------------
sql_pp_stats <- sprintf("
SELECT
  'Human' AS Species,
  p.Disease_presence_peptide,
  b.Disease,
  b.Isotype,
  b.BSource,
  b.BType,
  COUNT(DISTINCT b.Peptide) AS n_unique_peptides
FROM read_parquet('%s') b
JOIN read_parquet('%s') p
  ON b.Peptide = p.Peptide
GROUP BY 1,2,3,4,5,6
ORDER BY 3,4,5,6,2;
", pp_base_glob, pp_pres_glob)

assert_no_backslash(sql_pp_stats, "sql_pp_stats")
tbl_unique_peptides <- as.data.table(dbGetQuery(con, sql_pp_stats))

# -----------------------------
# Save
# -----------------------------
setcolorder(tbl_unique_antibodies,
            c("Species","Disease_presence_antibody","Disease","Isotype","BSource","BType","n_unique_antibodies"))
setcolorder(tbl_unique_peptides,
            c("Species","Disease_presence_peptide","Disease","Isotype","BSource","BType","n_unique_peptides"))

write.csv(tbl_unique_antibodies, pjoin(work_root, "sankey_data_antibody.csv"), row.names = FALSE)
write.csv(tbl_unique_peptides,  pjoin(work_root, "sankey_data_peptide.csv"),  row.names = FALSE)

cat("\nDONE.\n")
cat("Wrote:\n")
cat(" - ", pjoin(work_root, "sankey_data_antibody.csv"), "\n", sep = "")
cat(" - ", pjoin(work_root, "sankey_data_peptide.csv"),  "\n", sep = "")

