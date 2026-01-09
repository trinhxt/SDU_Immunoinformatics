################################################################################
## scripts/03_part3_disease_presence.R
##
## Part 3: Disease presence tables + OPTIONAL antibody Sankey stats
##
## Sankey peptide stats are NOT produced here (handled in Part 6).
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
  FORCE_REBUILD      <- FALSE
  WRITE_ANTIBODY_CSV <- TRUE
  
  N_BUCKETS       <- 256L
  DUCKDB_THREADS <- 6L
  DUCKDB_MEM     <- "44GB"
  
  REQUIRED_AB_COLS <- c(
    "Species",
    "Disease_presence_antibody",
    "Disease",
    "Isotype",
    "BSource",
    "BType",
    "n_unique_antibodies"
  )
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  part_done <- function(out_dir) {
    dir.exists(out_dir) &&
      length(list.files(out_dir, pattern = "\\.parquet$", recursive = TRUE)) > 0
  }
  
  sankey_ab_ok <- function(path) {
    if (!file.exists(path)) return(FALSE)
    hdr <- tryCatch(names(fread(path, nrows = 0)), error = function(e) NULL)
    !is.null(hdr) && all(REQUIRED_AB_COLS %in% hdr)
  }
  
  # ---------------------------------------------------------------------------
  # Input (Part2 output)
  # ---------------------------------------------------------------------------
  final_root <- norm(P$final_part2_dir)
  parquet_glob <- paste0(final_root, "/**/*.parquet")
  
  if (!dir.exists(final_root)) {
    stop("Part2 final dataset not found: ", final_root)
  }
  assert_no_backslash(parquet_glob, "parquet_glob")
  
  # ---------------------------------------------------------------------------
  # Output root
  # ---------------------------------------------------------------------------
  work_root <- norm(P$part3_root)
  dir.create(work_root, showWarnings = FALSE, recursive = TRUE)
  
  duckdb_file <- pjoin(work_root, "build.duckdb")
  duckdb_temp <- pjoin(work_root, "duckdb_temp")
  dir.create(duckdb_temp, showWarnings = FALSE, recursive = TRUE)
  
  ab_base_dir <- pjoin(work_root, "base_distinct", "antibody")
  pp_base_dir <- pjoin(work_root, "base_distinct", "peptide")
  ab_pres_dir <- pjoin(work_root, "presence", "antibody")
  pp_pres_dir <- pjoin(work_root, "presence", "peptide")
  
  dir.create(ab_base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(pp_base_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(ab_pres_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(pp_pres_dir, showWarnings = FALSE, recursive = TRUE)
  
  # ---------------------------------------------------------------------------
  # DuckDB connection (core tables)
  # ---------------------------------------------------------------------------
  con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  
  dbExecute(con, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
  dbExecute(con, "PRAGMA enable_progress_bar=true;")
  
  # ---------------------------------------------------------------------------
  # [1/4] Antibody DISTINCT base
  # ---------------------------------------------------------------------------
  if (!FORCE_REBUILD && part_done(ab_base_dir)) {
    cat("[1/4] Antibody DISTINCT base -> SKIP\n")
  } else {
    cat("[1/4] Writing antibody DISTINCT base...\n")
    
    sql <- sprintf("
      COPY (
        SELECT DISTINCT
          Disease, Isotype, BSource, BType, Antibody,
          (hash(Antibody) %% %d)::INTEGER AS bucket
        FROM read_parquet('%s')
        WHERE Disease IS NOT NULL AND Disease <> '' AND Disease <> 'None'
          AND Antibody IS NOT NULL AND Antibody <> ''
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, parquet_glob, ab_base_dir)
    
    dbExecute(con, sql)
  }
  
  # ---------------------------------------------------------------------------
  # [2/4] Peptide DISTINCT base
  # ---------------------------------------------------------------------------
  if (!FORCE_REBUILD && part_done(pp_base_dir)) {
    cat("[2/4] Peptide DISTINCT base -> SKIP\n")
  } else {
    cat("[2/4] Writing peptide DISTINCT base...\n")
    
    sql <- sprintf("
      COPY (
        SELECT DISTINCT
          Disease, Isotype, BSource, BType, Peptide,
          (hash(Peptide) %% %d)::INTEGER AS bucket
        FROM read_parquet('%s')
        WHERE Disease IS NOT NULL AND Disease <> '' AND Disease <> 'None'
          AND Peptide IS NOT NULL AND Peptide <> ''
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, parquet_glob, pp_base_dir)
    
    dbExecute(con, sql)
  }
  
  # ---------------------------------------------------------------------------
  # [3/4] Antibody presence
  # ---------------------------------------------------------------------------
  if (!FORCE_REBUILD && part_done(ab_pres_dir)) {
    cat("[3/4] Antibody presence -> SKIP\n")
  } else {
    cat("[3/4] Writing antibody presence...\n")
    
    sql <- sprintf("
      COPY (
        SELECT
          Antibody,
          (hash(Antibody) %% %d)::INTEGER AS bucket,
          COUNT(DISTINCT Disease) AS Disease_presence_antibody
        FROM read_parquet('%s')
        GROUP BY Antibody, bucket
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, paste0(ab_base_dir, "/**/*.parquet"), ab_pres_dir)
    
    dbExecute(con, sql)
  }
  
  # ---------------------------------------------------------------------------
  # [4/4] Peptide presence
  # ---------------------------------------------------------------------------
  if (!FORCE_REBUILD && part_done(pp_pres_dir)) {
    cat("[4/4] Peptide presence -> SKIP\n")
  } else {
    cat("[4/4] Writing peptide presence...\n")
    
    sql <- sprintf("
      COPY (
        SELECT
          Peptide,
          (hash(Peptide) %% %d)::INTEGER AS bucket,
          COUNT(DISTINCT Disease) AS Disease_presence_peptide
        FROM read_parquet('%s')
        GROUP BY Peptide, bucket
      ) TO '%s'
      (FORMAT PARQUET, PARTITION_BY(bucket), COMPRESSION ZSTD);
    ", N_BUCKETS, paste0(pp_base_dir, "/**/*.parquet"), pp_pres_dir)
    
    dbExecute(con, sql)
  }
  
  cat("\nPart3 core parquet tables DONE.\n\n")
  
  # ---------------------------------------------------------------------------
  # OPTIONAL: sankey_data_antibody.csv (schema-checked)
  # ---------------------------------------------------------------------------
  out_csv <- pjoin(work_root, "sankey_data_antibody.csv")
  
  if (!WRITE_ANTIBODY_CSV) {
    cat("WRITE_ANTIBODY_CSV=FALSE -> skip antibody CSV\n")
  } else if (sankey_ab_ok(out_csv)) {
    cat("sankey_data_antibody.csv exists with correct schema -> SKIP\n")
  } else {
    
    cat("Writing sankey_data_antibody.csv ...\n")
    
    con2 <- dbConnect(duckdb::duckdb(), dbdir = pjoin(work_root, "query.duckdb"))
    on.exit(try(dbDisconnect(con2, shutdown = TRUE), silent = TRUE), add = TRUE)
    
    sql <- sprintf("
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
    ", paste0(ab_base_dir, "/**/*.parquet"),
                   paste0(ab_pres_dir, "/**/*.parquet"))
    
    dt <- as.data.table(dbGetQuery(con2, sql))
    setcolorder(dt, REQUIRED_AB_COLS)
    fwrite(dt, out_csv)
    
    cat("Wrote:\n  ", out_csv, "\n", sep = "")
  }
  
  invisible(TRUE)
}

main()
