################################################################################
## scripts/02_part2_digestion.R
##
## Part2: Build peptide–antibody Parquet DB from processed OAS files
##
## Uses config/paths.yml (via helpers_paths.R):
##   - P$oas_processed_dir : input processed OAS *.csv.gz (from Part1)
##   - P$metadata_csv      : metadata CSV (from Part1)
##   - P$work_root         : working root for Part2 outputs
##   - P$final_part2_dir   : final compacted Parquet dataset directory
##
## Outputs (under P$work_root):
##   - digest_csv_gz/                   (Section 1 outputs; one per input file)
##   - parquet_db_partitioned/_staging/ (Section 2A per-file staging folders)
##   - parquet_db_partitioned/final/    (Section 2B final compacted dataset)
##   - Log_files/                       (resume logs)
##
## Update (2026-01):
## - After successful compaction, delete _staging to free disk space
##   (only if final dataset is verified to contain Parquet files).
################################################################################

# ==============================================================================
# Load config + helpers
# ==============================================================================
source("R/helpers_paths.R")
cfg <- load_config()
P   <- get_paths(cfg)

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Biostrings)
  library(cleaver)
  library(DBI)
  library(duckdb)
  library(future)
  library(future.apply)
  library(dplyr)
})

main <- function() {
  
  cat("Arrow version:", as.character(packageVersion("arrow")), "\n\n")
  
  # ---------------------------------------------------------------------------
  # Settings
  # ---------------------------------------------------------------------------
  FORCE_REBUILD_S1  <- FALSE  # set TRUE to redo digestion CSV.GZ
  FORCE_REBUILD_S2  <- FALSE  # set TRUE to redo staging (per-file parquet)
  FORCE_COMPACT     <- FALSE  # set TRUE to redo final compaction (delete final first is best)
  
  # Workers (Section 1 + Section 2A)
  N_WORKERS <- 6L
  
  # DuckDB compaction settings
  DUCKDB_THREADS <- 24L
  DUCKDB_MEM     <- "44GB"
  ROW_GROUP_SIZE <- 1000000L
  
  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------
  final_dataset_ok <- function(final_root) {
    dir.exists(final_root) &&
      length(list.files(final_root, pattern = "\\.parquet$", recursive = TRUE)) > 0
  }
  
  # ---------------------------------------------------------------------------
  # Paths from config
  # ---------------------------------------------------------------------------
  in_dir    <- P$oas_processed_dir
  work_root <- P$work_root
  
  if (!dir.exists(in_dir)) stop("Input processed dir not found: ", in_dir)
  if (!file.exists(P$metadata_csv)) stop("Metadata CSV not found: ", P$metadata_csv)
  
  # Section 1 output
  digest_csv_dir <- pjoin(work_root, "digest_csv_gz")
  dir.create(digest_csv_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Section 2 output roots
  parquet_db_dir <- pjoin(work_root, "parquet_db_partitioned")
  staging_root   <- pjoin(parquet_db_dir, "_staging")
  final_root     <- norm(P$final_part2_dir)  # allow config override
  dir.create(parquet_db_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(staging_root,   showWarnings = FALSE, recursive = TRUE)
  dir.create(final_root,     showWarnings = FALSE, recursive = TRUE)
  
  # Logs
  log_dir <- pjoin(work_root, "Log_files")
  dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
  
  processed_log_s1 <- pjoin(log_dir, "Section1_digestion_processed_files.txt")
  failed_log_s1    <- pjoin(log_dir, "Section1_digestion_failed_files.txt")
  processed_log_s2 <- pjoin(log_dir, "Section2_parquet_processed_files.txt")
  failed_log_s2    <- pjoin(log_dir, "Section2_parquet_failed_files.txt")
  
  for (f in c(processed_log_s1, failed_log_s1, processed_log_s2, failed_log_s2)) {
    if (!file.exists(f)) file.create(f)
  }
  
  processed_s1 <- readLines(processed_log_s1, warn = FALSE)
  processed_s2 <- readLines(processed_log_s2, warn = FALSE)
  
  cat("Input processed OAS dir:\n  ", in_dir, "\n", sep = "")
  cat("Work root:\n  ", work_root, "\n", sep = "")
  cat("Digest CSV dir:\n  ", digest_csv_dir, "\n", sep = "")
  cat("Staging root:\n  ", staging_root, "\n", sep = "")
  cat("Final Parquet dataset:\n  ", final_root, "\n\n", sep = "")
  
  # ---------------------------------------------------------------------------
  # 1) Reference peptide filter set
  # ---------------------------------------------------------------------------
  ref_rdata <- P$reference_tryptic_rdata
  if (!file.exists(ref_rdata)) {
    stop(
      "Reference peptide file not found:\n  ", ref_rdata, "\n\n",
      "See data/reference/README.md for instructions to build it."
    )
  }
  
  cat("Loading reference peptides:\n  ",
      normalizePath(ref_rdata, winslash = "/", mustWork = FALSE), "\n\n", sep = "")
  
  load(ref_rdata) # creates UniProtNCBI_Tryptic
  ref_dt <- unique(data.table(Peptide = as.character(UniProtNCBI_Tryptic)))
  setkey(ref_dt, Peptide)
  rm(UniProtNCBI_Tryptic)
  gc()
  
  # ---------------------------------------------------------------------------
  # 2) Metadata (select disease-mapped files + partition constants)
  # ---------------------------------------------------------------------------
  meta_cols <- c("Filename", "Disease", "BSource", "BType", "Isotype")
  metadata  <- fread(P$metadata_csv, select = meta_cols)
  
  metadata[, (meta_cols) := lapply(.SD, as.character), .SDcols = meta_cols]
  
  metadata <- metadata[
    !is.na(Disease) & Disease != "None" & Disease != "" &
      !is.na(BSource) & BSource != "" &
      !is.na(BType)   & BType   != "" &
      !is.na(Isotype) & Isotype != ""
  ]
  setkey(metadata, Filename)
  
  # ---------------------------------------------------------------------------
  # 3) Input discovery (only files present in metadata)
  # ---------------------------------------------------------------------------
  all_files <- list.files(in_dir, full.names = TRUE, recursive = FALSE, pattern = "\\.csv\\.gz$")
  if (!length(all_files)) stop("No *.csv.gz files found in: ", in_dir)
  
  task_dt <- data.table(file_path = all_files, filename = basename(all_files))
  task_dt <- task_dt[filename %in% metadata$Filename]
  setorder(task_dt, filename)
  
  cat("Total disease-mapped input files:", nrow(task_dt), "\n\n")
  
  # ---------------------------------------------------------------------------
  # SECTION 1: DIGESTION -> digest_csv_gz/*.csv.gz  (parallel)
  # ---------------------------------------------------------------------------
  KEEP_INPUT_COLS_S1 <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr3_aa")
  S1_OUT_COLS <- c("Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa")
  
  digest_one_to_csv_gz <- function(file_path, fn, ref_dt, out_dir, keep_cols) {
    
    tab <- arrow::read_csv_arrow(file_path, col_select = keep_cols)
    if (!("sequence_alignment_aa" %in% names(tab))) stop("Missing sequence_alignment_aa in: ", fn)
    
    dt <- as.data.table(tab)
    dt[, sequence_alignment_aa := as.character(sequence_alignment_aa)]
    
    for (cc in setdiff(keep_cols, "sequence_alignment_aa")) {
      if (!cc %in% names(dt)) dt[, (cc) := NA_character_]
    }
    dt[, `:=`(
      v_call  = as.character(v_call),
      d_call  = as.character(d_call),
      j_call  = as.character(j_call),
      cdr3_aa = as.character(cdr3_aa)
    )]
    
    seq_vec <- unique(dt$sequence_alignment_aa[!is.na(dt$sequence_alignment_aa) & nzchar(dt$sequence_alignment_aa)])
    
    out_file <- file.path(out_dir, fn)
    tmp_file <- paste0(out_file, ".tmp_", Sys.getpid())
    
    if (!length(seq_vec)) {
      out_dt <- data.table(Peptide=character(), Antibody=character(), v_call=character(),
                           d_call=character(), j_call=character(), cdr3_aa=character())
      fwrite(out_dt, tmp_file, sep = ",", compress = "gzip")
      file.rename(tmp_file, out_file)
      return(invisible(TRUE))
    }
    
    aa  <- Biostrings::AAStringSet(stats::setNames(seq_vec, seq_vec))
    dig <- cleaver::cleave(aa, enzym = "trypsin", missedCleavages = 0:1, unique = TRUE)
    
    lens    <- lengths(dig)
    pep_seq <- as.character(unlist(dig, use.names = FALSE))
    pep_ab  <- rep.int(names(dig), lens)
    out_dt  <- data.table(Peptide = pep_seq, Antibody = pep_ab)
    
    out_dt <- unique(out_dt[nchar(Peptide) > 4], by = c("Peptide", "Antibody"))
    setkey(out_dt, Peptide)
    out_dt <- out_dt[!ref_dt, on = "Peptide"]
    setkey(out_dt, NULL)
    
    ann_dt <- unique(dt[
      !is.na(sequence_alignment_aa) & nzchar(sequence_alignment_aa),
      .(sequence_alignment_aa, v_call, d_call, j_call, cdr3_aa)
    ], by = "sequence_alignment_aa")
    setkey(ann_dt, sequence_alignment_aa)
    
    out_dt <- ann_dt[out_dt, on = c(sequence_alignment_aa = "Antibody")]
    out_dt[, `:=`(Antibody = sequence_alignment_aa, sequence_alignment_aa = NULL)]
    
    for (cc in S1_OUT_COLS) if (!cc %in% names(out_dt)) out_dt[, (cc) := NA_character_]
    out_dt <- out_dt[, ..S1_OUT_COLS]
    
    fwrite(out_dt, tmp_file, sep = ",", compress = "gzip")
    file.rename(tmp_file, out_file)
    
    invisible(TRUE)
  }
  
  cat("========== PART2 / SECTION 1: DIGESTION -> CSV.GZ ==========\n")
  cat("Output dir:", digest_csv_dir, "\n\n")
  
  s1_tasks <- task_dt[, .(file_path, filename)]
  s1_tasks[, out_path := file.path(digest_csv_dir, filename)]
  
  if (FORCE_REBUILD_S1) {
    cat("FORCE_REBUILD_S1=TRUE: re-digest all files (will overwrite outputs).\n\n")
  } else {
    s1_tasks <- s1_tasks[!(filename %in% processed_s1 & file.exists(out_path))]
  }
  
  cat("Section 1 files to digest now:", nrow(s1_tasks), "\n\n")
  
  if (nrow(s1_tasks) == 0) {
    cat("Section 1: nothing to do.\n\n")
  } else {
    if (.Platform$OS.type == "windows") plan(multisession, workers = N_WORKERS) else plan(multicore, workers = N_WORKERS)
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
  
  # ---------------------------------------------------------------------------
  # SECTION 2: STAGING + COMPACTION
  # ---------------------------------------------------------------------------
  cat("========== PART2 / SECTION 2: STAGING + COMPACTION ==========\n")
  cat("Input digest CSV dir:", digest_csv_dir, "\n")
  cat("Staging root        :", staging_root, "\n")
  cat("Final root          :", final_root, "\n\n")
  
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
  
  safe_folder_name <- function(x) gsub("[\\\\/:*?\"<>|]", "_", x)
  
  digest_files <- list.files(digest_csv_dir, pattern = "\\.csv\\.gz$", full.names = TRUE)
  digest_dt <- data.table(file_path = digest_files, filename = basename(digest_files))
  setorder(digest_dt, filename)
  
  cat("Digest CSV.GZ files found:", nrow(digest_dt), "\n\n")
  
  s2_tasks <- digest_dt
  if (!FORCE_REBUILD_S2) {
    s2_tasks <- s2_tasks[!(filename %in% processed_s2)]
  } else {
    cat("FORCE_REBUILD_S2=TRUE: will attempt to stage all digest files.\n\n")
  }
  
  cat("Section 2 files to STAGE now:", nrow(s2_tasks), "\n\n")
  
  # 2A) Parallel staging
  if (nrow(s2_tasks) == 0) {
    cat("Section 2 staging: nothing to do.\n\n")
  } else {
    if (.Platform$OS.type == "windows") plan(multisession, workers = N_WORKERS) else plan(multicore, workers = N_WORKERS)
    options(future.globals.maxSize = 10 * 1024^3)
    
    s2_results <- future_lapply(seq_len(nrow(s2_tasks)), function(i) {
      fp <- s2_tasks$file_path[i]
      fn <- s2_tasks$filename[i]
      
      m <- metadata[J(fn)]
      if (nrow(m) == 0) {
        return(list(filename = fn, ok = FALSE, msg = paste0("No metadata row found (cannot partition): ", fn)))
      }
      m <- m[1]
      
      stage_dir <- file.path(staging_root, paste0("file_", safe_folder_name(fn)))
      dir.create(stage_dir, showWarnings = FALSE, recursive = TRUE)
      
      if (!FORCE_REBUILD_S2 &&
          length(list.files(stage_dir, pattern = "\\.parquet$", recursive = TRUE)) > 0) {
        return(list(filename = fn, ok = TRUE, msg = "staged_exists"))
      }
      
      ok <- TRUE
      msg <- "ok"
      tryCatch({
        tab <- arrow::read_csv_arrow(fp)
        dt  <- as.data.table(tab)
        
        dt[, filename := fn]
        dt[, `:=`(
          Disease = as.character(m$Disease),
          BSource = as.character(m$BSource),
          BType   = as.character(m$BType),
          Isotype = as.character(m$Isotype)
        )]
        
        dt <- coerce_schema_s2(dt)
        
        do.call(arrow::write_dataset, c(
          list(
            dataset = dt,
            path = stage_dir,
            format = "parquet",
            partitioning = PARTITION_COLS,
            existing_data_behavior = "overwrite"
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
  
  # 2B) Compaction (DuckDB, sequential)
  final_has_any <- dir.exists(final_root) && length(list.files(final_root, recursive = TRUE)) > 0
  if (final_has_any && !FORCE_COMPACT) {
    cat("Compaction SKIP: final_root already has data.\n")
    cat("  ", final_root, "\n\n", sep = "")
  } else {
    
    if (final_has_any && FORCE_COMPACT) {
      cat("FORCE_COMPACT=TRUE but final_root is not empty.\n")
      cat("Recommendation: delete final_root before compaction to avoid mixing old/new:\n  ",
          final_root, "\n\n", sep = "")
    }
    
    cat("Starting compaction (DuckDB) -> final dataset:\n  ", final_root, "\n\n", sep = "")
    
    duckdb_temp_compact <- pjoin(work_root, "duckdb_temp_part2_compact")
    dir.create(duckdb_temp_compact, showWarnings = FALSE, recursive = TRUE)
    
    conC <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
    on.exit(try(dbDisconnect(conC, shutdown = TRUE), silent = TRUE), add = TRUE)
    
    dbExecute(conC, sprintf("PRAGMA threads=%d;", DUCKDB_THREADS))
    dbExecute(conC, sprintf("PRAGMA memory_limit='%s';", DUCKDB_MEM))
    dbExecute(conC, sprintf("PRAGMA temp_directory='%s';", gsub("\\\\", "/", duckdb_temp_compact)))
    dbExecute(conC, "PRAGMA enable_progress_bar=true;")
    
    staging_glob <- paste0(gsub("\\\\", "/", staging_root), "/**/*.parquet")
    final_path   <- gsub("\\\\", "/", final_root)
    
    assert_no_backslash(staging_glob, "staging_glob")
    assert_no_backslash(final_path, "final_path")
    
    sql_compact <- sprintf("
      COPY (
        SELECT *
        FROM read_parquet('%s')
      )
      TO '%s'
      (FORMAT PARQUET,
       PARTITION_BY (Disease, BSource, BType, Isotype),
       COMPRESSION ZSTD,
       ROW_GROUP_SIZE %d);
    ", staging_glob, final_path, ROW_GROUP_SIZE)
    
    dbExecute(conC, sql_compact)
    
    cat("Compaction finished.\nChecking final dataset...\n")
    
    if (final_dataset_ok(final_root)) {
      cat("Final dataset verified (Parquet files found).\n")
      cat("Final dataset root:\n  ", final_root, "\n\n", sep = "")
      
      if (dir.exists(staging_root)) {
        cat("Removing staging directory to free disk space:\n  ",
            staging_root, "\n", sep = "")
        unlink(staging_root, recursive = TRUE, force = TRUE)
      }
      
    } else {
      cat("WARNING: Final dataset check FAILED. Staging preserved for debugging:\n  ",
          staging_root, "\n\n", sep = "")
    }
  }
  
  # Optional sanity check (light)
  cat("Sanity check: open final dataset (Arrow) and show schema preview\n")
  ds_final <- arrow::open_dataset(final_root, format = "parquet")
  glimpse(ds_final)
  cat("\n")
  
  invisible(TRUE)
}

main()
