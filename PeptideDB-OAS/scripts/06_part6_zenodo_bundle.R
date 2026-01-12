################################################################################
## scripts/06_part6_zenodo_bundle.R
##
## Part6: Build Zenodo-ready PeptideDB
##
## Output structure:
##   D:/OAS/PeptideDB/
##     ├─ final_enriched/      (partitioned parquet, query this)
##     ├─ derived_tables/      (copied once, small peptide-level tables)
##     ├─ README.txt
##     ├─ build_part6.duckdb
##     └─ duckdb_temp_part6/
##
## Design:
##   - Partition-wise enrichment (safe, low-memory)
##   - One final compaction step
##   - Resume-safe
##
## Key fix:
##   - DuckDB COPY ... TO 'path' expects a FILE path (unless PARTITION_BY is used).
##   - We now write stage output into stage_dir/data.parquet (a real file),
##     instead of trying to COPY TO a directory.
################################################################################
################################################################################
## scripts/06_part6_zenodo_bundle.R
##
## Part6: Build Zenodo-ready PeptideDB
##
## SPEED UPGRADE:
##   1) Load global peptide “dimension tables” once into DuckDB:
##        - dim_p3: Disease_presence_peptide
##        - dim_p4: n_distinct_antibodies
##        - dim_p5: peptide_in_cdr3
##      + create indexes on Peptide
##
##   2) For each Part2 partition, first compute DISTINCT peptide keys from that
##      partition, then FILTER dim tables to those keys before joining.
##      This avoids scanning/joining huge global tables for every partition.
##
## Key fix kept:
##   - Per-partition staging writes to stage_dir/data.parquet (file, not directory)
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
  # SETTINGS
  # ---------------------------------------------------------------------------
  FORCE_REBUILD <- FALSE
  COPY_DERIVED  <- TRUE
  
  # Increase threads if your CPU allows (e.g., 8/12/16)
  DUCKDB_THREADS <- 8L
  DUCKDB_MEM     <- "44GB"
  
  # Larger row groups often write faster (try 2e6–5e6 if you have RAM)
  ROW_GROUP_SIZE <- 2e6
  
  # Whether to build/rebuild the DuckDB dimension tables
  # (Set TRUE if your part3/4/5 inputs changed and you want to refresh dims)
  REBUILD_DIMS <- FALSE
  
  # ---------------------------------------------------------------------------
  # INPUT PATHS
  # ---------------------------------------------------------------------------
  part2_final_root <- norm(P$final_part2_dir)
  
  part3_pp_pres_dir <- pjoin(norm(P$part3_root), "presence", "peptide")
  part4_counts_dir  <- pjoin(norm(P$part4_root), "counts", "peptide")
  part5_answer_dir  <- pjoin(norm(P$part5_root), "answer", "peptide")
  
  # DuckDB glob paths (always forward slashes)
  part3_glob <- paste0(part3_pp_pres_dir, "/**/*.parquet")
  part4_glob <- paste0(part4_counts_dir,  "/**/*.parquet")
  part5_glob <- paste0(part5_answer_dir,  "/**/*.parquet")
  
  # ---------------------------------------------------------------------------
  # OUTPUT PATHS (Zenodo bundle)
  # ---------------------------------------------------------------------------
  out_root   <- norm(P$zenodo_out_root)
  out_final  <- pjoin(out_root, "final_enriched")
  out_stage  <- pjoin(out_root, "_staging_enriched")
  out_readme <- pjoin(out_root, "README.txt")
  
  duckdb_file <- pjoin(out_root, "build_part6.duckdb")
  duckdb_temp <- pjoin(out_root, "duckdb_temp_part6")
  dir.create(duckdb_temp, recursive = TRUE, showWarnings = FALSE)
  
  cat("Part6 output root:\n  ", out_root, "\n\n", sep = "")
  
  # ---------------------------------------------------------------------------
  # 1) COPY DERIVED TABLES (ONCE)
  # ---------------------------------------------------------------------------
  if (COPY_DERIVED) {
    derived_root <- pjoin(out_root, "derived_tables")
    
    if (!FORCE_REBUILD &&
        dir.exists(derived_root) &&
        length(list.files(derived_root, recursive = TRUE, all.files = TRUE)) > 0) {
      
      cat("derived_tables already exists with content -> SKIP copy\n\n")
      
    } else {
      
      cat("Copying derived tables into:\n  ", derived_root, "\n\n", sep = "")
      unlink(derived_root, recursive = TRUE, force = TRUE)
      dir.create(derived_root, recursive = TRUE, showWarnings = FALSE)
      
      file.copy(part3_pp_pres_dir,
                pjoin(derived_root, "disease_presence_peptide"),
                recursive = TRUE)
      
      file.copy(part4_counts_dir,
                pjoin(derived_root, "peptide_uniqueness"),
                recursive = TRUE)
      
      file.copy(part5_answer_dir,
                pjoin(derived_root, "peptide_in_cdr3"),
                recursive = TRUE)
    }
  }
  
  # ---------------------------------------------------------------------------
  # 2) PREPARE STAGING + FINAL OUTPUT
  # ---------------------------------------------------------------------------
  if (FORCE_REBUILD) {
    unlink(out_stage, recursive = TRUE, force = TRUE)
    unlink(out_final, recursive = TRUE, force = TRUE)
  }
  
  if (dir.exists(out_final) &&
      length(list.files(out_final, recursive = TRUE, all.files = TRUE)) > 0) {
    
    cat("final_enriched already exists -> SKIP build\n")
    return(invisible(TRUE))
  }
  
  dir.create(out_stage, recursive = TRUE, showWarnings = FALSE)
  dir.create(out_final, recursive = TRUE, showWarnings = FALSE)
  
  # ---------------------------------------------------------------------------
  # 3) DISCOVER PART2 PARTITIONS
  # ---------------------------------------------------------------------------
  part_dirs <- list.dirs(part2_final_root, recursive = TRUE, full.names = TRUE)
  part_dirs <- part_dirs[
    grepl("Disease=", part_dirs) &
      grepl("BSource=", part_dirs) &
      grepl("BType=", part_dirs) &
      grepl("Isotype=", part_dirs)
  ]
  
  cat("Found Part2 partitions:", length(part_dirs), "\n\n")
  
  # ---------------------------------------------------------------------------
  # 4) DUCKDB SETUP
  # ---------------------------------------------------------------------------
  con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file)
  on.exit(try(dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
  
  dbExecute(con, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
  dbExecute(con, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
  dbExecute(con, sprintf("PRAGMA temp_directory='%s';", duckdb_temp))
  dbExecute(con, "PRAGMA preserve_insertion_order=false;")
  dbExecute(con, "PRAGMA enable_progress_bar=true;")
  
  # ---------------------------------------------------------------------------
  # 4b) BUILD/LOAD GLOBAL DIM TABLES ONCE (major speedup)
  # ---------------------------------------------------------------------------
  # These tables turn "read_parquet(glob) join ..." repeated 333x
  # into fast indexed joins inside DuckDB.
  #
  # We keep only needed columns to reduce size.
  #
  # NOTE: If these dims already exist, we reuse them (unless REBUILD_DIMS=TRUE).
  # ---------------------------------------------------------------------------
  table_exists <- function(con, tbl) {
    q <- sprintf("
      SELECT COUNT(*) AS n
      FROM information_schema.tables
      WHERE table_schema='main' AND table_name='%s'
    ", tbl)
    as.integer(dbGetQuery(con, q)$n) > 0
  }
  
  if (REBUILD_DIMS) {
    dbExecute(con, "DROP TABLE IF EXISTS dim_p3;")
    dbExecute(con, "DROP TABLE IF EXISTS dim_p4;")
    dbExecute(con, "DROP TABLE IF EXISTS dim_p5;")
  }
  
  if (!table_exists(con, "dim_p3")) {
    cat("Building dim_p3 from Part3 parquet...\n")
    dbExecute(con, sprintf("
      CREATE TABLE dim_p3 AS
      SELECT
        Peptide,
        Disease_presence_peptide::INTEGER AS Disease_presence_peptide
      FROM read_parquet('%s');
    ", part3_glob))
  } else {
    cat("dim_p3 exists -> reuse\n")
  }
  
  if (!table_exists(con, "dim_p4")) {
    cat("Building dim_p4 from Part4 parquet...\n")
    dbExecute(con, sprintf("
      CREATE TABLE dim_p4 AS
      SELECT
        Peptide,
        n_distinct_antibodies::BIGINT AS n_distinct_antibodies
      FROM read_parquet('%s');
    ", part4_glob))
  } else {
    cat("dim_p4 exists -> reuse\n")
  }
  
  if (!table_exists(con, "dim_p5")) {
    cat("Building dim_p5 from Part5 parquet...\n")
    dbExecute(con, sprintf("
      CREATE TABLE dim_p5 AS
      SELECT
        Peptide,
        peptide_in_cdr3::INTEGER AS peptide_in_cdr3
      FROM read_parquet('%s');
    ", part5_glob))
  } else {
    cat("dim_p5 exists -> reuse\n")
  }
  
  # Indexes help when filtering/joining dims by Peptide
  # (DuckDB may or may not use them in all cases, but generally helps)
  cat("Creating indexes on dim tables...\n")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_dim_p3_pep ON dim_p3(Peptide);")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_dim_p4_pep ON dim_p4(Peptide);")
  dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_dim_p5_pep ON dim_p5(Peptide);")
  
  cat("\n")
  
  # ---------------------------------------------------------------------------
  # 5) PARTITION-WISE ENRICHMENT (SAFE LOOP) - FAST VERSION
  #    - Extract DISTINCT peptide keys for this partition
  #    - Filter dim tables to just those keys
  #    - Join filtered dims to the partition rows
  #    - Write stage_dir/data.parquet
  # ---------------------------------------------------------------------------
  part2_root_fwd <- gsub("\\\\", "/", part2_final_root)
  out_stage_fwd  <- gsub("\\\\", "/", out_stage)
  
  for (i in seq_along(part_dirs)) {
    
    src_part <- gsub("\\\\", "/", part_dirs[i])
    
    # Build relative partition path (starts with "/Disease=.../BSource=.../...")
    rel_part <- sub(part2_root_fwd, "", src_part, fixed = TRUE)
    
    # Staging directory for this partition
    stage_dir  <- paste0(out_stage_fwd, rel_part)
    stage_file <- paste0(stage_dir, "/data.parquet")
    
    # Resume-safe: skip only if output file exists
    if (file.exists(stage_file)) next
    
    dir.create(stage_dir, recursive = TRUE, showWarnings = FALSE)
    
    cat(sprintf("[%d/%d] Enriching partition:\n  %s\n",
                i, length(part_dirs), src_part))
    
    # Fast SQL:
    # - keys: distinct peptides present in this partition
    # - p3/p4/p5: restrict dims to those peptides only
    # - join those small filtered dims back to the partition rows
    sql_part <- sprintf("
      COPY (
        WITH keys AS (
          SELECT DISTINCT Peptide
          FROM read_parquet('%s/**/*.parquet')
          WHERE Peptide IS NOT NULL
        ),
        p3 AS (
          SELECT d.Peptide, d.Disease_presence_peptide
          FROM dim_p3 d
          INNER JOIN keys k USING (Peptide)
        ),
        p4 AS (
          SELECT d.Peptide, d.n_distinct_antibodies
          FROM dim_p4 d
          INNER JOIN keys k USING (Peptide)
        ),
        p5 AS (
          SELECT d.Peptide, d.peptide_in_cdr3
          FROM dim_p5 d
          INNER JOIN keys k USING (Peptide)
        )
        SELECT
          f.*,
          COALESCE(p3.Disease_presence_peptide, 0)::INTEGER AS Disease_presence_peptide,
          COALESCE(p4.n_distinct_antibodies, 0)::BIGINT     AS n_distinct_antibodies,
          COALESCE(p5.peptide_in_cdr3, 0)::INTEGER          AS peptide_in_cdr3
        FROM read_parquet('%s/**/*.parquet') f
        LEFT JOIN p3 USING (Peptide)
        LEFT JOIN p4 USING (Peptide)
        LEFT JOIN p5 USING (Peptide)
      )
      TO '%s'
      (FORMAT PARQUET,
       COMPRESSION ZSTD,
       ROW_GROUP_SIZE %d);
    ",
                        src_part,
                        src_part,
                        stage_file,
                        as.integer(ROW_GROUP_SIZE)
    )
    
    dbExecute(con, sql_part)
  }
  
  # ---------------------------------------------------------------------------
  # 6) FINAL COMPACTION (partitioned output)
  # ---------------------------------------------------------------------------
  cat("\nFinal compaction -> final_enriched\n")
  
  out_final_fwd <- gsub("\\\\", "/", out_final)
  
  sql_final <- sprintf("
    COPY (
      SELECT * FROM read_parquet('%s/**/*.parquet')
    )
    TO '%s'
    (FORMAT PARQUET,
     PARTITION_BY (Disease, BSource, BType, Isotype),
     COMPRESSION ZSTD,
     ROW_GROUP_SIZE %d);
  ",
                       out_stage_fwd,
                       out_final_fwd,
                       as.integer(ROW_GROUP_SIZE)
  )
  
  dbExecute(con, sql_final)
  
  # Cleanup staging
  unlink(out_stage, recursive = TRUE, force = TRUE)
  
  # ---------------------------------------------------------------------------
  # 7) README
  # ---------------------------------------------------------------------------
  writeLines(c(
    "PeptideDB (Zenodo-ready)",
    "",
    "final_enriched/",
    "  Partitioned Parquet dataset (query this).",
    "",
    "derived_tables/",
    "  disease_presence_peptide/",
    "  peptide_uniqueness/",
    "  peptide_in_cdr3/",
    "",
    "Enriched columns:",
    "  - Disease_presence_peptide",
    "  - n_distinct_antibodies",
    "  - peptide_in_cdr3",
    "",
    "Partitioning:",
    "  Disease / BSource / BType / Isotype",
    "",
    "Build notes:",
    "  - Dimension tables (dim_p3/dim_p4/dim_p5) are loaded once into DuckDB",
    "  - Each partition filters dims by its peptide keys before joining (faster)"
  ), out_readme)
  
  cat("\nDONE (Part6)\n")
  cat("Final dataset:\n  ", out_final, "\n", sep = "")
  cat("README:\n  ", out_readme, "\n", sep = "")
  
  invisible(TRUE)
}

main()
