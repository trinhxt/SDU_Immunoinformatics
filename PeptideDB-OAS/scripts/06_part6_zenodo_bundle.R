################################################################################
## scripts/06_part6_zenodo_bundle.R
##
## Part6: Create Zenodo-shareable PeptideDB folder (single root)
##
## Goal
## ----
## Create a single folder (can be zipped) containing:
##   1) final_enriched/  (partitioned Parquet)
##      = Part2 final dataset + peptide-level annotations:
##          - Disease_presence_peptide   (Part3)
##          - n_distinct_antibodies      (Part4)
##          - peptide_in_cdr3            (Part5)
##
##   2) derived_tables/ (optional copies of small peptide-level tables)
##      - disease_presence_peptide/
##      - peptide_uniqueness/
##      - peptide_in_cdr3/
##
##   3) sankey_data_peptide.csv (created here; Part3 already created antibody one)
##
##   4) README.txt (schema + provenance)
##
## Uses config/paths.yml (via helpers_paths.R):
##   - P$final_part2_dir  : Part2 final dataset
##   - P$part3_root       : Part3 outputs
##   - P$part4_root       : Part4 outputs
##   - P$part5_root       : Part5 outputs
##   - P$zenodo_out_root  : output folder for sharing (e.g., D:/OAS/PeptideDB)
##
## Resume/Skip
## ----------
## - If final_enriched/ already has parquet files, skip unless FORCE_REBUILD=TRUE.
## - derived_tables/ copy is idempotent (overwrites files).
##
## Notes
## -----
## - Windows-safe for DuckDB: forward slashes only (helpers)
## - Uses DuckDB COPY for a compact, partitioned output.
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
  library(data.table)
})

main <- function() {
  
  # ---------------------------------------------------------------------------
  # Settings
  # ---------------------------------------------------------------------------
  FORCE_REBUILD <- FALSE     # rebuild final_enriched even if exists
  COPY_DERIVED  <- TRUE      # copy small derived tables into PeptideDB
  BUILD_SANKEY_PEPTIDE <- TRUE
  
  DUCKDB_THREADS <- 6L
  DUCKDB_MEM     <- "44GB"
  ROW_GROUP_SIZE <- 1000000L
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  dir_has_parquet <- function(d) {
    if (!dir.exists(d)) return(FALSE)
    length(list.files(d, pattern = "\\.parquet$", recursive = TRUE, full.names = TRUE)) > 0
  }
  
  copy_dir_recursive <- function(from, to) {
    if (!dir.exists(from)) stop("Missing folder: ", from)
    dir.create(to, showWarnings = FALSE, recursive = TRUE)
    files <- list.files(from, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
    for (f in files) {
      rel <- substring(f, nchar(from) + 2)  # +2 drops "/"
      dest <- file.path(to, rel)
      dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
      file.copy(f, dest, overwrite = TRUE, copy.mode = TRUE, copy.date = TRUE)
    }
    invisible(TRUE)
  }
  
  # ---------------------------------------------------------------------------
  # Inputs
  # ---------------------------------------------------------------------------
  part2_final_root <- norm(P$final_part2_dir)
  if (!dir.exists(part2_final_root)) stop("Part2 final dataset not found: ", part2_final_root)
  part2_glob <- paste0(part2_final_root, "/**/*.parquet")
  
  # Part3 peptide presence
  part3_pp_pres_dir  <- pjoin(norm(P$part3_root), "presence", "peptide")
  if (!dir.exists(part3_pp_pres_dir)) stop("Part3 peptide presence dir not found: ", part3_pp_pres_dir)
  part3_pp_pres_glob <- paste0(part3_pp_pres_dir, "/**/*.parquet")
  
  # Part4 peptide uniqueness counts
  part4_counts_dir  <- pjoin(norm(P$part4_root), "counts", "peptide")
  if (!dir.exists(part4_counts_dir)) stop("Part4 peptide counts dir not found: ", part4_counts_dir)
  part4_counts_glob <- paste0(part4_counts_dir, "/**/*.parquet")
  
  # Part5 peptide_in_cdr3 answers
  part5_answer_dir  <- pjoin(norm(P$part5_root), "answer", "peptide")
  if (!dir.exists(part5_answer_dir)) stop("Part5 peptide answer dir not found: ", part5_answer_dir)
  part5_answer_glob <- paste0(part5_answer_dir, "/**/*.parquet")
  
  # Sanity: forward slashes for DuckDB globs
  assert_no_backslash(part2_glob, "part2_glob")
  assert_no_backslash(part3_pp_pres_glob, "part3_pp_pres_glob")
  assert_no_backslash(part4_counts_glob, "part4_counts_glob")
  assert_no_backslash(part5_answer_glob, "part5_answer_glob")
  
  # ---------------------------------------------------------------------------
  # Outputs
  # ---------------------------------------------------------------------------
  out_root <- norm(P$zenodo_out_root)
  dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
  
  out_final_enriched <- pjoin(out_root, "final_enriched")
  out_derived_root   <- pjoin(out_root, "derived_tables")
  out_readme         <- pjoin(out_root, "README.txt")
  
  out_sankey_peptide_csv <- pjoin(out_root, "sankey_data_peptide.csv")
  
  # DuckDB build artifacts (kept in out_root so the bundle is self-contained)
  duckdb_file <- pjoin(out_root, "build_part6.duckdb")
  duckdb_temp <- pjoin(out_root, "duckdb_temp_part6")
  dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)
  
  assert_no_backslash(out_root, "out_root")
  assert_no_backslash(out_final_enriched, "out_final_enriched")
  assert_no_backslash(duckdb_file, "duckdb_file")
  assert_no_backslash(duckdb_temp, "duckdb_temp")
  
  cat("Part6 output root:\n  ", out_root, "\n\n", sep = "")
  cat("Inputs:\n")
  cat("  Part2 final: ", part2_final_root, "\n", sep = "")
  cat("  Part3 peptide presence: ", part3_pp_pres_dir, "\n", sep = "")
  cat("  Part4 peptide uniqueness: ", part4_counts_dir, "\n", sep = "")
  cat("  Part5 peptide_in_cdr3: ", part5_answer_dir, "\n\n", sep = "")
  
  # ---------------------------------------------------------------------------
  # 4) Optional: copy derived tables
  # ---------------------------------------------------------------------------
  if (isTRUE(COPY_DERIVED)) {
    cat("Copying derived tables into:\n  ", out_derived_root, "\n\n", sep = "")
    dir.create(out_derived_root, showWarnings = FALSE, recursive = TRUE)
    
    copy_dir_recursive(part3_pp_pres_dir, pjoin(out_derived_root, "disease_presence_peptide"))
    copy_dir_recursive(part4_counts_dir,  pjoin(out_derived_root, "peptide_uniqueness"))
    copy_dir_recursive(part5_answer_dir,  pjoin(out_derived_root, "peptide_in_cdr3"))
  }
  
  # ---------------------------------------------------------------------------
  # 5) Build enriched final dataset (DuckDB COPY)
  # ---------------------------------------------------------------------------
  if (!FORCE_REBUILD && dir_has_parquet(out_final_enriched)) {
    cat("final_enriched already exists -> SKIP (set FORCE_REBUILD=TRUE to rebuild)\n\n")
  } else {
    
    # Safety: refuse to rebuild into non-empty dir unless FORCE_REBUILD=TRUE
    if (dir.exists(out_final_enriched) &&
        length(list.files(out_final_enriched, recursive = TRUE)) > 0 &&
        !FORCE_REBUILD) {
      stop(
        "Output not empty and FORCE_REBUILD=FALSE:\n  ", out_final_enriched, "\n",
        "Delete final_enriched/ or set FORCE_REBUILD=TRUE."
      )
    }
    
    dir.create(out_final_enriched, showWarnings = FALSE, recursive = TRUE)
    
    cat("Building enriched partitioned dataset:\n")
    cat("  Output: ", out_final_enriched, "\n\n", sep = "")
    
    con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
    on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
    
    dbExecute(con, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
    dbExecute(con, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
    dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
    dbExecute(con, "PRAGMA enable_progress_bar=true;")
    
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
       ROW_GROUP_SIZE %d);
    ", part2_glob, part3_pp_pres_glob, part4_counts_glob, part5_answer_glob, out_final_path, ROW_GROUP_SIZE)
    
    assert_no_backslash(sql_copy, "sql_copy")
    dbExecute(con, sql_copy)
    
    if (!dir_has_parquet(out_final_enriched)) {
      stop("Build failed: no parquet files found in ", out_final_enriched)
    }
    
    cat("\nEnriched dataset build DONE.\n\n")
  }
  
  # ---------------------------------------------------------------------------
  # 5B) Build sankey_data_peptide.csv (for plotting later)
  # ---------------------------------------------------------------------------
  if (isTRUE(BUILD_SANKEY_PEPTIDE)) {
    
    # Skip if already valid
    if (file.exists(out_sankey_peptide_csv)) {
      ok_cols <- c("Species","Disease_presence_peptide","Disease","Isotype","BSource","BType","n_unique_peptides")
      dt_try <- tryCatch(fread(out_sankey_peptide_csv, nrows = 0), error = function(e) NULL)
      if (!is.null(dt_try) && all(ok_cols %in% names(dt_try))) {
        cat("sankey_data_peptide.csv exists with expected columns -> SKIP\n\n")
      } else {
        cat("sankey_data_peptide.csv exists but schema mismatch -> rebuild\n")
        file.remove(out_sankey_peptide_csv)
      }
    }
    
    if (!file.exists(out_sankey_peptide_csv)) {
      cat("Building sankey_data_peptide.csv...\n")
      
      # Use derived peptide-level disease presence + Part2 partitions for Disease/Isotype/BSource/BType context
      con2 <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
      on.exit(try(dbDisconnect(con2, shutdown = TRUE), silent = TRUE), add = TRUE)
      
      dbExecute(con2, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
      dbExecute(con2, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
      dbExecute(con2, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
      dbExecute(con2, "PRAGMA enable_progress_bar=true;")
      
      sql_pp_stats <- sprintf("
        SELECT
          'Human' AS Species,
          p.Disease_presence_peptide,
          f.Disease,
          f.Isotype,
          f.BSource,
          f.BType,
          COUNT(DISTINCT f.Peptide) AS n_unique_peptides
        FROM read_parquet('%s') f
        LEFT JOIN read_parquet('%s') p
          ON f.Peptide = p.Peptide
        GROUP BY 1,2,3,4,5,6
        ORDER BY 3,4,5,6,2;
      ", part2_glob, part3_pp_pres_glob)
      
      dt_pp <- as.data.table(dbGetQuery(con2, sql_pp_stats))
      
      setcolorder(dt_pp, c(
        "Species","Disease_presence_peptide","Disease","Isotype","BSource","BType","n_unique_peptides"
      ))
      
      fwrite(dt_pp, out_sankey_peptide_csv)
      cat("Wrote:\n  ", out_sankey_peptide_csv, "\n\n", sep = "")
    }
  }
  
  # ---------------------------------------------------------------------------
  # 6) Write README
  # ---------------------------------------------------------------------------
  readme_lines <- c(
    "PeptideDB (for Zenodo sharing)",
    "",
    "This folder was created by scripts/06_part6_zenodo_bundle.R (Part6).",
    "",
    "Folders:",
    "  final_enriched/      Partitioned Parquet dataset (query this).",
    "  derived_tables/      (Optional) small peptide-level tables copied from Parts 3–5.",
    "",
    "Files:",
    "  sankey_data_peptide.csv  Convenience summary for Sankey plotting (created in Part6).",
    "  build_part6.duckdb       DuckDB file used to build final_enriched (optional for recipients).",
    "  README.txt               This file.",
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
  
  cat("DONE (Part6).\n")
  cat("PeptideDB root:\n  ", out_root, "\n", sep = "")
  cat("Enriched dataset:\n  ", out_final_enriched, "\n", sep = "")
  if (isTRUE(COPY_DERIVED)) cat("Derived tables:\n  ", out_derived_root, "\n", sep = "")
  if (isTRUE(BUILD_SANKEY_PEPTIDE)) cat("Sankey peptide CSV:\n  ", out_sankey_peptide_csv, "\n", sep = "")
  cat("README:\n  ", out_readme, "\n", sep = "")
  
  invisible(TRUE)
}

main()
