# ==============================================================================
# Script: Part6_PeptideDB_Zenodo.R
# Date: 2026-01
#
# GOAL
# ----
# Create a single, shareable "PeptideDB" folder for Zenodo:
#   D:\OAS\PeptideDB
#
# Contents:
#   1) final_enriched/   (partitioned Parquet, same partitions as Part2 final)
#      - Part2 final dataset + 3 peptide-level annotations:
#          Disease_presence_peptide   (Part3)
#          n_distinct_antibodies      (Part4)
#          peptide_in_cdr3            (Part5)
#
#   2) derived_tables/ (optional copies of the small peptide-level tables)
#      - disease_presence_peptide/
#      - peptide_uniqueness/
#      - peptide_in_cdr3/
#
#   3) README.txt with schema/provenance hints
#
# NOTES
# -----
# - Uses DuckDB COPY to create a new compacted partitioned dataset.
# - Windows-safe paths: forward slashes for DuckDB.
# - Resume/skip: if output already has parquet files, skip unless FORCE_REBUILD=TRUE.
# ==============================================================================

suppressPackageStartupMessages({
  library(DBI)
  library(duckdb)
})

# ==============================================================================
# 0) SETTINGS
# ==============================================================================
FORCE_REBUILD <- FALSE     # set TRUE to rebuild outputs even if exist
COPY_DERIVED  <- TRUE      # copy small derived tables into PeptideDB folder

# DuckDB resources
DUCKDB_THREADS <- 6L
DUCKDB_MEM     <- "44GB"

# Buckets used in your derived tables (must match Parts 3-5)
N_BUCKETS <- 256L

# ==============================================================================
# 1) HELPERS (Windows-safe)
# ==============================================================================
norm  <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)
pjoin <- function(...) norm(paste(..., sep = "/"))

assert_no_backslash <- function(x, name) {
  if (grepl("\\\\", x)) stop("Backslash found in ", name, ":\n", x)
}

dir_has_parquet <- function(d) {
  if (!dir.exists(d)) return(FALSE)
  length(list.files(d, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)) > 0
}

copy_dir_recursive <- function(from, to) {
  if (!dir.exists(from)) stop("Missing folder: ", from)
  dir.create(to, showWarnings = FALSE, recursive = TRUE)
  files <- list.files(from, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
  for (f in files) {
    rel <- substring(f, nchar(from) + 2)   # +2 to drop "/"
    dest <- file.path(to, rel)
    dest_dir <- dirname(dest)
    dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)
    # overwrite=TRUE to make reruns robust
    file.copy(f, dest, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
  }
  invisible(TRUE)
}

# ==============================================================================
# 2) INPUT PATHS (edit if your folders differ)
# ==============================================================================
# Part2 final parquet dataset
part2_final_root <- norm("D:/OAS/OAS_human_heavychain_disease_tryptic/parquet_db_partitioned/final")
part2_glob <- paste0(part2_final_root, "/**/*.parquet")

# Part3 peptide disease presence (output from Part3)
part3_work_root <- norm("OAS_human_heavychain_disease_tryptic/Disease_presence")
part3_pp_pres_dir <- pjoin(part3_work_root, "presence", "peptide")
part3_pp_pres_glob <- paste0(part3_pp_pres_dir, "/**/*.parquet")

# Part4 peptide uniqueness (output from Part4)
part4_work_root <- norm("OAS_human_heavychain_disease_tryptic/Peptide_uniqueness")
part4_counts_dir <- pjoin(part4_work_root, "counts", "peptide")
part4_counts_glob <- paste0(part4_counts_dir, "/**/*.parquet")

# Part5 peptide-in-CDR3 (output from Part5)
part5_work_root <- norm("OAS_human_heavychain_disease_tryptic/CDR3_check")
part5_answer_dir <- pjoin(part5_work_root, "answer", "peptide")
part5_answer_glob <- paste0(part5_answer_dir, "/**/*.parquet")

# Sanity: ensure forward slashes for DuckDB globs
for (x in list(
  part2_glob = part2_glob,
  part3_pp_pres_glob = part3_pp_pres_glob,
  part4_counts_glob = part4_counts_glob,
  part5_answer_glob = part5_answer_glob
)) assert_no_backslash(x[[1]], names(x))

# ==============================================================================
# 3) OUTPUT PATHS (Zenodo share folder)
# ==============================================================================
out_root <- norm("D:/OAS/PeptideDB")
dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

out_final_enriched <- pjoin(out_root, "final_enriched")   # partitioned parquet output
out_derived_root   <- pjoin(out_root, "derived_tables")

out_readme <- pjoin(out_root, "README.txt")

# DuckDB files for this build (kept inside PeptideDB so the whole folder is self-contained)
duckdb_file <- pjoin(out_root, "build_part6.duckdb")
duckdb_temp <- pjoin(out_root, "duckdb_temp_part6")
dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)

assert_no_backslash(out_root, "out_root")
assert_no_backslash(out_final_enriched, "out_final_enriched")
assert_no_backslash(duckdb_file, "duckdb_file")
assert_no_backslash(duckdb_temp, "duckdb_temp")

# ==============================================================================
# 4) OPTIONAL: COPY DERIVED TABLES (small, reusable)
# ==============================================================================
if (COPY_DERIVED) {
  cat("Copying derived tables into:\n  ", out_derived_root, "\n\n", sep = "")
  dir.create(out_derived_root, showWarnings = FALSE, recursive = TRUE)
  
  # Disease presence peptide
  copy_dir_recursive(part3_pp_pres_dir, pjoin(out_derived_root, "disease_presence_peptide"))
  
  # Peptide uniqueness
  copy_dir_recursive(part4_counts_dir, pjoin(out_derived_root, "peptide_uniqueness"))
  
  # Peptide in CDR3
  copy_dir_recursive(part5_answer_dir, pjoin(out_derived_root, "peptide_in_cdr3"))
}

# ==============================================================================
# 5) BUILD ENRICHED FINAL DATASET (DuckDB COPY with joins)
# ==============================================================================
if (!FORCE_REBUILD && dir_has_parquet(out_final_enriched)) {
  cat("\nfinal_enriched already exists -> SKIP (set FORCE_REBUILD=TRUE to rebuild)\n")
} else {
  if (dir.exists(out_final_enriched) && length(list.files(out_final_enriched, recursive = TRUE)) > 0) {
    cat("\nWARNING: out_final_enriched is not empty:\n  ", out_final_enriched, "\n", sep = "")
    cat("Recommendation: delete it before rebuild to avoid mixing old/new.\n\n")
    if (!FORCE_REBUILD) stop("Output not empty and FORCE_REBUILD=FALSE. Aborting for safety.")
  }
  
  dir.create(out_final_enriched, showWarnings = FALSE, recursive = TRUE)
  
  cat("\nBuilding enriched partitioned dataset:\n")
  cat("  Input  (Part2): ", part2_final_root, "\n", sep = "")
  cat("  Output        : ", out_final_enriched, "\n\n", sep = "")
  
  con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  
  dbExecute(con, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
  dbExecute(con, "PRAGMA enable_progress_bar=true;")
  
  # IMPORTANT:
  # - read_parquet(part2_glob) is partitioned by Disease/BSource/BType/Isotype.
  # - the derived tables are peptide-level and bucket-partitioned; we just read them by glob.
  # - we LEFT JOIN so no Part2 rows are lost; COALESCE fills missing values.
  #
  # Output is re-partitioned by Disease/BSource/BType/Isotype and compacted into fewer files.
  
  out_final_path <- gsub("\\\\", "/", out_final_enriched)
  
  sql_copy <- sprintf("
    COPY (
      SELECT
        f.*,

        -- Part3: #distinct diseases per peptide
        CAST(COALESCE(p.Disease_presence_peptide, 0) AS INTEGER) AS Disease_presence_peptide,

        -- Part4: #distinct antibodies per peptide
        CAST(COALESCE(u.n_distinct_antibodies, 0) AS BIGINT)     AS n_distinct_antibodies,

        -- Part5: peptide contained in any observed CDR3 (1/0)
        CAST(COALESCE(c.peptide_in_cdr3, 0) AS INTEGER)          AS peptide_in_cdr3

      FROM read_parquet('%s') f
      LEFT JOIN read_parquet('%s') p
        ON f.Peptide = p.Peptide
      LEFT JOIN read_parquet('%s') u
        ON f.Peptide = u.Peptide
      LEFT JOIN read_parquet('%s') c
        ON f.Peptide = c.Peptide
    )
    TO '%s'
    (FORMAT PARQUET,
     PARTITION_BY (Disease, BSource, BType, Isotype),
     COMPRESSION ZSTD,
     ROW_GROUP_SIZE 1000000);
  ", part2_glob, part3_pp_pres_glob, part4_counts_glob, part5_answer_glob, out_final_path)
  
  assert_no_backslash(sql_copy, "sql_copy")
  
  dbExecute(con, sql_copy)
  
  if (!dir_has_parquet(out_final_enriched)) {
    stop("Build failed: no parquet files found in ", out_final_enriched)
  }
  
  cat("\nEnriched dataset build DONE.\n")
}

# ==============================================================================
# 6) WRITE README (lightweight provenance + schema hints)
# ==============================================================================
readme_lines <- c(
  "PeptideDB (for Zenodo sharing)",
  "",
  "This folder was created by Part6_PeptideDB_Zenodo.R.",
  "",
  "Folders:",
  "  final_enriched/      Partitioned Parquet dataset (query this).",
  "  derived_tables/      (Optional) small peptide-level tables copied from Parts 3–5.",
  "  build_part6.duckdb   DuckDB file used to build final_enriched (optional for recipients).",
  "",
  "final_enriched columns:",
  "  Contains all columns from Part2 final dataset, plus:",
  "    - Disease_presence_peptide   (INTEGER)  #distinct diseases peptide appears in (Part3)",
  "    - n_distinct_antibodies      (BIGINT)   #distinct antibodies peptide appears in (Part4)",
  "    - peptide_in_cdr3            (INTEGER)  1 if peptide found in ANY observed cdr3_aa (Part5)",
  "",
  "Partitioning:",
  "  final_enriched/ is partitioned by: Disease / BSource / BType / Isotype",
  "",
  "Provenance (inputs):",
  paste0("  Part2 final parquet: ", part2_final_root),
  paste0("  Part3 peptide presence: ", part3_pp_pres_dir),
  paste0("  Part4 peptide uniqueness: ", part4_counts_dir),
  paste0("  Part5 peptide_in_cdr3: ", part5_answer_dir),
  "",
  "Notes:",
  "  - Values are LEFT JOINed onto the Part2 rows by Peptide; missing values are filled with 0.",
  "  - This dataset may contain repeated peptides across partitions because Part2 is partitioned by metadata.",
  "  - For purely peptide-level analysis, use derived_tables/* (one row per peptide)."
)

writeLines(readme_lines, con = out_readme)

cat("\nDONE.\n")
cat("PeptideDB root:\n  ", out_root, "\n", sep = "")
cat("Enriched dataset:\n  ", out_final_enriched, "\n", sep = "")
if (COPY_DERIVED) cat("Derived tables:\n  ", out_derived_root, "\n", sep = "")
cat("README:\n  ", out_readme, "\n", sep = "")
