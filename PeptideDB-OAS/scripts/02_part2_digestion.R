# ==============================================================================
# Script: Part2_Digestion_OAS_human_heavychain_disease_tryptic_2sections_v1c.R
# Date: 2026-01
#
# GOAL
# ----
# Build a scalable “database-like” dataset from many OAS CSV.GZ files.
#
# Section 1 (PARALLEL, 6 workers)
#   For each disease-mapped input CSV.GZ file:
#     1) Read only required columns for digestion + annotations.
#     2) In-silico tryptic digestion (0–1 missed cleavages).
#     3) Remove peptides found in reference proteomes (UniProt/NCBI tryptic set).
#     4) Attach per-antibody annotations (v_call, d_call, j_call, cdr3_aa).
#     5) Write one compressed CSV.GZ output per input file.
#   Output columns ONLY:
#     Peptide, Antibody, v_call, d_call, j_call, cdr3_aa
#   Output filenames EXACTLY match input filenames (same *.csv.gz name).
#   This is critical because Section 2 joins metadata by Filename.
#
# Section 2 (PARALLEL staging + SEQUENTIAL compaction)
#   Stage step (parallel, safe):
#     1) Read each Section 1 digest CSV.GZ.
#     2) Add metadata constants by filename (Disease/BSource/BType/Isotype).
#     3) Write into a per-file staging folder as a partitioned Arrow dataset.
#        (Per-file folder avoids write contention between workers.)
#   Compaction step (DuckDB, sequential):
#     1) Read all staging parquet fragments across all per-file folders.
#     2) Rewrite into ONE final partitioned dataset:
#        partitioned by Disease/BSource/BType/Isotype.
#        This compacts small pieces into fewer larger Parquet files and
#        is faster for later DuckDB scans.
#
# Notes
# -----
# - Arrow 22.0.0.1 compatible:
#     write_dataset() supports existing_data_behavior = "overwrite" | "error" | "delete_matching"
#     No "file_options" argument.
# - Section 1 uses an atomic-ish write pattern: write tmp file -> rename.
#   This reduces the chance of leaving half-written files after crashes.
# ==============================================================================

library(data.table)
library(arrow)
library(Biostrings)
library(cleaver)
library(DBI)
library(duckdb)
library(future)
library(future.apply)
library(dplyr)

cat("Arrow version:", as.character(packageVersion("arrow")), "\n\n")

# ==============================================================================
# 0) PATHS + LOGS (pipeline is resumable)
# ==============================================================================
# Input directory containing original OAS *.csv.gz files (non-recursive listing below).
in_dir <- normalizePath("OAS_human_heavychain", winslash = "/", mustWork = FALSE)

# Root output directory for this pipeline
work_root <- normalizePath("OAS_human_heavychain_disease_tryptic", winslash = "/", mustWork = FALSE)

# Section 1 output: one digest CSV.GZ per input file (same filename as input)
digest_csv_dir <- file.path(work_root, "digest_csv_gz")
dir.create(digest_csv_dir, showWarnings = FALSE, recursive = TRUE)

# Section 2 output root:
#   parquet_db_dir/_staging : per-file staging folders (safe parallel writes)
#   parquet_db_dir/final    : final compacted partitioned dataset (query this)
parquet_db_dir <- file.path(work_root, "parquet_db_partitioned")
dir.create(parquet_db_dir, showWarnings = FALSE, recursive = TRUE)

staging_root <- file.path(parquet_db_dir, "_staging")
final_root   <- file.path(parquet_db_dir, "final")
dir.create(staging_root, showWarnings = FALSE, recursive = TRUE)
dir.create(final_root,   showWarnings = FALSE, recursive = TRUE)

# Logs allow safe reruns: processed files are skipped if outputs exist
log_dir <- file.path(work_root, "Log_files")
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)

processed_log_s1 <- file.path(log_dir, "Section1_digestion_processed_files.txt")
failed_log_s1    <- file.path(log_dir, "Section1_digestion_failed_files.txt")
processed_log_s2 <- file.path(log_dir, "Section2_parquet_processed_files.txt")
failed_log_s2    <- file.path(log_dir, "Section2_parquet_failed_files.txt")

for (f in c(processed_log_s1, failed_log_s1, processed_log_s2, failed_log_s2)) {
  if (!file.exists(f)) file.create(f)
}

processed_s1 <- readLines(processed_log_s1, warn = FALSE)
processed_s2 <- readLines(processed_log_s2, warn = FALSE)

# ==============================================================================
# 1) REFERENCE PEPTIDES (filter set)
# ==============================================================================
# Any peptide that exists in the reference tryptic peptide set is removed.
# This reduces peptides that may come from common proteome background.
load("UniProtNCBI_Tryptic.RData") # creates UniProtNCBI_Tryptic
ref_dt <- unique(data.table(Peptide = as.character(UniProtNCBI_Tryptic)))
setkey(ref_dt, Peptide)          # enables fast anti-join: out_dt[!ref_dt, on="Peptide"]
rm(UniProtNCBI_Tryptic)
gc()

# ==============================================================================
# 2) METADATA (used to select “disease-mapped” files + Section 2 partitioning)
# ==============================================================================
# Only keep what we need for dataset partitioning (and to decide which files to process).
meta_cols <- c("Filename", "Disease", "BSource", "BType", "Isotype")
metadata  <- fread("OAS_metadata.csv", select = meta_cols)

# Force consistent types for matching + partition folder creation
metadata[, (meta_cols) := lapply(.SD, as.character), .SDcols = meta_cols]

# Disease-only rows + require all partition values (avoid NA partitions)
metadata <- metadata[
  !is.na(Disease) & Disease != "None" & Disease != "" &
    !is.na(BSource) & BSource != "" &
    !is.na(BType)   & BType   != "" &
    !is.na(Isotype) & Isotype != ""
]
setkey(metadata, Filename)

# ==============================================================================
# 3) INPUT DISCOVERY (files to process are limited to those present in metadata)
# ==============================================================================
all_files <- list.files(in_dir, full.names = TRUE, recursive = FALSE)
if (!length(all_files)) stop("No files found in: ", in_dir)

task_dt <- data.table(file_path = all_files, filename = basename(all_files))
task_dt <- task_dt[filename %in% metadata$Filename]
setorder(task_dt, filename)

cat("Total disease-mapped input files:", nrow(task_dt), "\n\n")

# ==============================================================================
# SECTION 1: DIGESTION -> individual CSV.GZ files (PARALLEL)
# ==============================================================================
# Input columns needed:
# - sequence_alignment_aa is the antibody AA sequence to digest
# - v/d/j/cdr3 are per-antibody annotations that we attach to each peptide
KEEP_INPUT_COLS_S1 <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr3_aa")

# Output schema for Section 1
S1_OUT_COLS <- c("Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa")

digest_one_to_csv_gz <- function(file_path, fn, ref_dt, out_dir, keep_cols) {
  
  # ---- Read minimal columns (Arrow) ----
  tab <- arrow::read_csv_arrow(file_path, col_select = keep_cols)
  if (!("sequence_alignment_aa" %in% names(tab))) stop("Missing sequence_alignment_aa in: ", fn)
  
  dt <- as.data.table(tab)
  dt[, sequence_alignment_aa := as.character(sequence_alignment_aa)]
  
  # Some inputs may be missing certain annotation columns; create them as NA for schema stability.
  for (cc in setdiff(keep_cols, "sequence_alignment_aa")) {
    if (!cc %in% names(dt)) dt[, (cc) := NA_character_]
  }
  dt[, `:=`(
    v_call  = as.character(v_call),
    d_call  = as.character(d_call),
    j_call  = as.character(j_call),
    cdr3_aa = as.character(cdr3_aa)
  )]
  
  # ---- Build unique antibody list (digest unique sequences only) ----
  seq_vec <- unique(dt$sequence_alignment_aa[!is.na(dt$sequence_alignment_aa) & nzchar(dt$sequence_alignment_aa)])
  
  # Output filenames must match inputs (critical for Section 2 metadata join)
  out_file <- file.path(out_dir, fn)
  tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
  
  # If no sequences, write header-only gz CSV (still a valid artifact for “processed” status)
  if (!length(seq_vec)) {
    out_dt <- data.table(Peptide=character(), Antibody=character(), v_call=character(),
                         d_call=character(), j_call=character(), cdr3_aa=character())
    fwrite(out_dt, tmp_file, sep = ",", compress = "gzip")
    file.rename(tmp_file, out_file)
    return(invisible(TRUE))
  }
  
  # ---- Tryptic digestion (0–1 missed cleavages) ----
  aa  <- Biostrings::AAStringSet(stats::setNames(seq_vec, seq_vec))
  dig <- cleaver::cleave(aa, enzym = "trypsin", missedCleavages = 0:1, unique = TRUE)
  
  # Expand digest result into long form (Peptide, Antibody)
  lens    <- lengths(dig)
  pep_seq <- as.character(unlist(dig, use.names = FALSE))
  pep_ab  <- rep.int(names(dig), lens)
  out_dt  <- data.table(Peptide = pep_seq, Antibody = pep_ab)
  
  # ---- Filters ----
  # 1) remove very short peptides
  # 2) ensure unique (Peptide, Antibody) pairs
  out_dt <- unique(out_dt[nchar(Peptide) > 4], by = c("Peptide", "Antibody"))
  
  # 3) remove peptides found in reference proteomes
  setkey(out_dt, Peptide)
  out_dt <- out_dt[!ref_dt, on = "Peptide"]
  setkey(out_dt, NULL)
  
  # ---- Attach per-antibody annotations ----
  # Build a unique annotation table per antibody sequence (sequence_alignment_aa)
  ann_dt <- unique(dt[
    !is.na(sequence_alignment_aa) & nzchar(sequence_alignment_aa),
    .(sequence_alignment_aa, v_call, d_call, j_call, cdr3_aa)
  ], by = "sequence_alignment_aa")
  setkey(ann_dt, sequence_alignment_aa)
  
  # Join annotations onto peptide rows using Antibody == sequence_alignment_aa
  out_dt <- ann_dt[out_dt, on = c(sequence_alignment_aa = "Antibody")]
  out_dt[, `:=`(Antibody = sequence_alignment_aa, sequence_alignment_aa = NULL)]
  
  # Enforce fixed column order (stable schema across all output files)
  for (cc in S1_OUT_COLS) if (!cc %in% names(out_dt)) out_dt[, (cc) := NA_character_]
  out_dt <- out_dt[, ..S1_OUT_COLS]
  
  # Write gzipped CSV (standard gzip, 7zip/WinRAR friendly), using tmp -> rename to avoid partial files
  fwrite(out_dt, tmp_file, sep = ",", compress = "gzip")
  file.rename(tmp_file, out_file)
  
  invisible(TRUE)
}

cat("========== SECTION 1 (PARALLEL): DIGESTION -> CSV.GZ ==========\n")
cat("Output dir:", digest_csv_dir, "\n\n")

# Build Section 1 task list:
# Skip if:
#   - filename is in processed log AND
#   - output file exists
s1_tasks <- task_dt[, .(file_path, filename)]
s1_tasks[, out_path := file.path(digest_csv_dir, filename)]
s1_tasks <- s1_tasks[!(filename %in% processed_s1 & file.exists(out_path))]

cat("Section 1 files to digest now:", nrow(s1_tasks), "\n\n")

if (nrow(s1_tasks) == 0) {
  
  cat("Section 1: nothing to do (all files already processed and outputs exist).\n\n")
  s1_res_dt <- data.table(filename = character(), ok = logical(), msg = character())
  s1_ok  <- character()
  s1_bad <- character()
  
} else {
  
  n_workers <- 6L
  if (.Platform$OS.type == "windows") plan(multisession, workers = n_workers) else plan(multicore, workers = n_workers)
  options(future.globals.maxSize = 10 * 1024^3)
  
  s1_results <- future_lapply(seq_len(nrow(s1_tasks)), function(i) {
    fp <- s1_tasks$file_path[i]
    fn <- s1_tasks$filename[i]
    
    ok <- TRUE
    msg <- "ok"
    tryCatch({
      digest_one_to_csv_gz(fp, fn, ref_dt, digest_csv_dir, KEEP_INPUT_COLS_S1)
    }, error = function(e) {
      ok  <<- FALSE
      msg <<- conditionMessage(e)
    })
    list(filename = fn, ok = ok, msg = msg)
  })
  
  plan(sequential)
  
  s1_res_dt <- rbindlist(s1_results, fill = TRUE)
  s1_ok  <- s1_res_dt[ok == TRUE,  filename]
  s1_bad <- s1_res_dt[ok == FALSE, filename]
  
  if (length(s1_ok))  cat(paste(s1_ok,  collapse = "\n"), "\n", file = processed_log_s1, append = TRUE, sep = "")
  if (length(s1_bad)) cat(paste(s1_bad, collapse = "\n"), "\n", file = failed_log_s1,    append = TRUE, sep = "")
  
  cat(sprintf("\nSECTION 1 done. OK: %d | Failed: %d\n\n", length(s1_ok), length(s1_bad)))
  if (length(s1_bad)) {
    cat("First few Section 1 failures:\n")
    print(s1_res_dt[ok == FALSE][1:min(10, .N)])
  }
}

# ==============================================================================
# SECTION 2: BUILD PARTITIONED PARQUET DB
#   Stage (parallel): per-file staging folders (no contention)
#   Compact (duckdb): rewrite to one final partitioned dataset (fast to query)
# ==============================================================================
cat("========== SECTION 2 (PARALLEL STAGING + COMPACT): BUILD PARTITIONED PARQUET DB ==========\n")
cat("Input digest CSV dir:", digest_csv_dir, "\n")
cat("Parquet DB dir      :", parquet_db_dir, "\n\n")

PARTITION_COLS <- c("Disease", "BSource", "BType", "Isotype")
PARQUET_WRITE_ARGS <- list(
  compression = "zstd",
  use_dictionary = TRUE,
  write_statistics = TRUE,
  chunk_size = 100000
)

S2_OUT_COLS <- c("Peptide","Antibody","v_call","d_call","j_call","cdr3_aa",
                 "filename","Disease","BSource","BType","Isotype")

coerce_schema_s2 <- function(dt) {
  missing <- setdiff(S2_OUT_COLS, names(dt))
  if (length(missing)) for (cc in missing) dt[, (cc) := NA_character_]
  extra <- setdiff(names(dt), S2_OUT_COLS)
  if (length(extra)) dt[, (extra) := NULL]
  for (cc in S2_OUT_COLS) dt[, (cc) := as.character(dt[[cc]])]
  dt[, ..S2_OUT_COLS]
}

# Windows-safe folder naming (illegal characters replaced)
safe_folder_name <- function(x) gsub("[\\\\/:*?\"<>|]", "_", x)

# Build work list from Section 1 outputs
digest_files <- list.files(digest_csv_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
digest_dt <- data.table(file_path = digest_files, filename = basename(digest_files))
setorder(digest_dt, filename)

cat("Section 1 digest CSV.GZ files found:", nrow(digest_dt), "\n\n")

# Process only files not already logged as staged
s2_tasks <- digest_dt[!(filename %in% processed_s2)]
cat("Section 2 files to STAGE now:", nrow(s2_tasks), "\n\n")

# ==============================================================================
# 2A) PARALLEL STAGING
# ==============================================================================
if (nrow(s2_tasks) == 0) {
  
  cat("Section 2 staging: nothing to do (all files already staged per processed log).\n\n")
  s2_res_dt <- data.table(filename = character(), ok = logical(), msg = character())
  
} else {
  
  n_workers <- 6L
  if (.Platform$OS.type == "windows") plan(multisession, workers = n_workers) else plan(multicore, workers = n_workers)
  options(future.globals.maxSize = 10 * 1024^3)
  
  s2_results <- future_lapply(seq_len(nrow(s2_tasks)), function(i) {
    fp <- s2_tasks$file_path[i]
    fn <- s2_tasks$filename[i]
    
    # Look up metadata (partition constants for this file)
    m <- metadata[J(fn)]
    if (nrow(m) == 0) {
      return(list(filename = fn, ok = FALSE, msg = paste0("No metadata row found (cannot partition): ", fn)))
    }
    m <- m[1]
    
    # Each file writes into its own stage_dir => safe parallel writes
    stage_dir <- file.path(staging_root, paste0("file_", safe_folder_name(fn)))
    dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
    
    # Skip if parquet already exists in stage_dir (useful if logs were reset)
    if (length(list.files(stage_dir, pattern = "\\.parquet$", recursive = TRUE)) > 0) {
      return(list(filename = fn, ok = TRUE, msg = "staged_exists"))
    }
    
    ok  <- TRUE
    msg <- "ok"
    tryCatch({
      tab <- arrow::read_csv_arrow(fp)
      dt  <- as.data.table(tab)
      
      # Attach filename and metadata constants
      dt[, filename := fn]
      dt[, `:=`(
        Disease = as.character(m$Disease),
        BSource = as.character(m$BSource),
        BType   = as.character(m$BType),
        Isotype = as.character(m$Isotype)
      )]
      
      dt <- coerce_schema_s2(dt)
      
      # Write a partitioned dataset *inside this stage_dir*
      # existing_data_behavior="error" protects against accidental overwrite of non-empty stage_dir
      do.call(arrow::write_dataset, c(
        list(
          dataset = dt,
          path = stage_dir,
          format = "parquet",
          partitioning = PARTITION_COLS,
          existing_data_behavior = "error"
        ),
        PARQUET_WRITE_ARGS
      ))
    }, error = function(e) {
      ok  <<- FALSE
      msg <<- conditionMessage(e)
    })
    
    list(filename = fn, ok = ok, msg = msg)
  })
  
  plan(sequential)
  
  s2_res_dt <- rbindlist(s2_results, fill = TRUE)
  s2_ok  <- s2_res_dt[ok == TRUE,  filename]
  s2_bad <- s2_res_dt[ok == FALSE, filename]
  
  if (length(s2_ok))  cat(paste(s2_ok,  collapse = "\n"), "\n", file = processed_log_s2, append = TRUE, sep = "")
  if (length(s2_bad)) cat(paste(s2_bad, collapse = "\n"), "\n", file = failed_log_s2,    append = TRUE, sep = "")
  
  cat(sprintf("\nSECTION 2 staging done. OK: %d | Failed: %d\n\n", length(s2_ok), length(s2_bad)))
  if (length(s2_bad)) {
    cat("First few staging failures:\n")
    print(s2_res_dt[ok == FALSE][1:min(10, .N)])
  }
}

cat("Staging dataset root:\n  ", staging_root, "\n\n", sep = "")

# ==============================================================================
# 2B) ONE-TIME COMPACTION (DuckDB, sequential)
# ==============================================================================
# Reads all staged parquet fragments and rewrites them as ONE final dataset
# partitioned by the same columns. This reduces file counts and improves scan speed.
cat("Starting compaction (DuckDB) -> final dataset:\n  ", final_root, "\n\n", sep = "")

duckdb_temp <- normalizePath(file.path(getwd(), "duckdb_temp"), winslash = "/", mustWork = FALSE)
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
dbExecute(con, "PRAGMA threads=24;")
dbExecute(con, "PRAGMA memory_limit='44GB';")
dbExecute(con, sprintf("PRAGMA temp_directory='%s';", gsub("\\\\", "/", duckdb_temp)))

if (dir.exists(final_root) && length(list.files(final_root, recursive = TRUE)) > 0) {
  cat("WARNING: final_root is not empty:\n  ", final_root, "\n", sep = "")
  cat("Recommendation: delete it before re-running compaction to avoid mixing old/new.\n\n")
}

staging_glob <- paste0(gsub("\\\\", "/", staging_root), "/**/*.parquet")
final_path   <- gsub("\\\\", "/", final_root)

sql_compact <- sprintf("
  COPY (
    SELECT *
    FROM read_parquet('%s')
  )
  TO '%s'
  (FORMAT PARQUET,
   PARTITION_BY (Disease, BSource, BType, Isotype),
   COMPRESSION ZSTD,
   ROW_GROUP_SIZE 1000000);
", staging_glob, final_path)

dbExecute(con, sql_compact)
dbDisconnect(con, shutdown = TRUE)

cat("Compaction done.\nFinal dataset root:\n  ", final_root, "\n\n", sep = "")

# Sanity check: open final dataset with Arrow and inspect schema/preview
ds_final <- arrow::open_dataset(final_root, format = "parquet")
glimpse(ds_final)

ds_stage <- arrow::open_dataset(staging_root, format = "parquet")
glimpse(ds_stage)
# ==============================================================================
