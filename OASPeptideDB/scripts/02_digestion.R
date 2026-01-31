################################################################################
## scripts/02_digestion.R
##
## Part2: Build peptide–antibody Parquet DB from processed OAS files
##
## Logic:
## 1. Digestion (CSV -> CSV):
##    - Reads processed OAS files (sequence_alignment_aa)
##    - In-silico digestion (Trypsin, 0/1 missed cleavages)
##    - Calculates Peptide Start/End indices
##    - Flags "CDR3-Spanning" peptides using Coordinate Overlap
##    - Filters against UniProt reference background
##
## 2. Staging (CSV -> Parquet):
##    - Converts digested CSVs to Parquet partitions (Disease/Isotype/etc)
##
## 3. Compaction (Parquet -> Parquet):
##    - Merges many small staging files into optimized final dataset
################################################################################

# ==============================================================================
# Load Config
# ==============================================================================
source("scripts/00_config.R")

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
  
  cat("Arrow version:", as.character(packageVersion("arrow")), "\n")
  
  # ---------------------------------------------------------------------------
  # Settings
  # ---------------------------------------------------------------------------
  FORCE_REBUILD_S1  <- FALSE  # Redo digestion?
  FORCE_REBUILD_S2  <- FALSE  # Redo staging?
  FORCE_COMPACT     <- FALSE  # Redo final compaction?
  
  # Ensure sub-directories exist
  dir.create(P$digest_csv_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(P$staging_root,   showWarnings = FALSE, recursive = TRUE)
  dir.create(P$db_state1,      showWarnings = FALSE, recursive = TRUE)
  dir.create(P$log_dir,        showWarnings = FALSE, recursive = TRUE)
  
  # Log files
  logs <- list(
    s1_proc = file.path(P$log_dir, "Section1_digestion_processed_files.txt"),
    s1_fail = file.path(P$log_dir, "Section1_digestion_failed_files.txt"),
    s2_proc = file.path(P$log_dir, "Section2_parquet_processed_files.txt"),
    s2_fail = file.path(P$log_dir, "Section2_parquet_failed_files.txt")
  )
  for (f in logs) if (!file.exists(f)) file.create(f)
  
  # ---------------------------------------------------------------------------
  # 1) Load Reference Filter
  # ---------------------------------------------------------------------------
  if (!file.exists(P$ref_rdata)) stop("Ref file missing: ", P$ref_rdata)
  
  cat("Loading reference peptides...\n")
  load(P$ref_rdata) # Expects object 'UniProtNCBI_Tryptic'
  ref_dt <- unique(data.table(Peptide = as.character(UniProtNCBI_Tryptic)))
  setkey(ref_dt, Peptide)
  rm(UniProtNCBI_Tryptic); gc()
  
  # ---------------------------------------------------------------------------
  # 2) Load Metadata
  # ---------------------------------------------------------------------------
  meta_cols <- c("Filename", "Disease", "BSource", "BType", "Isotype")
  if (!file.exists(P$metadata_csv)) stop("Metadata missing: ", P$metadata_csv)
  
  metadata <- fread(P$metadata_csv, select = meta_cols)
  metadata[, (meta_cols) := lapply(.SD, as.character), .SDcols = meta_cols]
  
  # Filter for valid metadata
  metadata <- metadata[
    !is.na(Disease) & Disease != "" &
      !is.na(BSource) & BSource != "" &
      !is.na(BType)   & BType   != "" &
      !is.na(Isotype) & Isotype != ""
  ]
  setkey(metadata, Filename)
  
  # ---------------------------------------------------------------------------
  # 3) Identify Input Files
  # ---------------------------------------------------------------------------
  all_files <- list.files(P$processed_dir, full.names = TRUE, pattern = "\\.csv\\.gz$")
  task_dt <- data.table(file_path = all_files, filename = basename(all_files))
  
  # Only process files that are in our cleaned metadata
  task_dt <- task_dt[filename %in% metadata$Filename]
  setorder(task_dt, filename)
  
  cat("Files to process (Matched Metadata):", nrow(task_dt), "\n\n")
  
  # ===========================================================================
  # SECTION 1: DIGESTION
  # ===========================================================================
  cat("=== SECTION 1: DIGESTION ===\n")
  
  # Columns to read from input
  KEEP_INPUT_COLS <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr3_aa")
  # Columns to write to output
  S1_OUT_COLS <- c("Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa", 
                   "Start", "End", "Is_CDR3_Spanning")
  
  # Define the worker function
  digest_worker <- function(file_path, fn, ref_dt, out_dir) {
    
    # Read Arrow -> DT
    tab <- arrow::read_csv_arrow(file_path, col_select = KEEP_INPUT_COLS)
    dt  <- as.data.table(tab)
    
    # Basic Cleaning
    dt <- dt[!is.na(sequence_alignment_aa) & nzchar(sequence_alignment_aa)]
    if (nrow(dt) == 0) return(invisible(TRUE)) # Empty file handled silently
    
    # Prepare unique sequences for digestion
    # We map back later using 'sequence_alignment_aa' as the key
    seq_vec <- unique(dt$sequence_alignment_aa)
    
    # Biostrings / Cleaver Digestion
    aa  <- Biostrings::AAStringSet(stats::setNames(seq_vec, seq_vec))
    dig <- cleaver::cleave(aa, enzym = "trypsin", missedCleavages = 0:1, unique = TRUE)
    
    # Expand Digestion Results
    lens    <- lengths(dig)
    pep_seq <- as.character(unlist(dig, use.names = FALSE))
    pep_ab  <- rep.int(names(dig), lens) # The parent sequence
    
    # ----------------------------------------------------------
    # CRITICAL UPDATE: Calculate Indices & CDR3 Overlap
    # ----------------------------------------------------------
    
    # 1. Find Start Index (Vectorized)
    # regexpr(fixed=TRUE) is fast and strictly literal
    match_indices <- regexpr(pep_seq, pep_ab, fixed = TRUE)
    starts <- as.integer(match_indices)
    ends   <- starts + nchar(pep_seq) - 1
    
    pep_dt <- data.table(
      Peptide  = pep_seq,
      Antibody = pep_ab,
      Start    = starts,
      End      = ends
    )
    
    # 2. Filter Filter
    # Remove short peptides (< 6 AA)
    pep_dt <- pep_dt[nchar(Peptide) >= 6]
    # Remove UniProt Reference Peptides (Background filter)
    setkey(pep_dt, Peptide)
    pep_dt <- pep_dt[!ref_dt] # Anti-join
    
    # 3. Join back Metadata (V/D/J/CDR3)
    # 'Antibody' col in pep_dt is the full sequence, which matches 'sequence_alignment_aa' in dt
    dt_meta <- unique(dt[, .(sequence_alignment_aa, v_call, d_call, j_call, cdr3_aa)])
    pep_dt <- dt_meta[pep_dt, on = c(sequence_alignment_aa = "Antibody")]
    
    # Rename back for clarity
    setnames(pep_dt, "sequence_alignment_aa", "Antibody")
    
    # 4. CDR3 Spanning Logic (Coordinate Overlap)
    # Find where CDR3 is inside the Antibody
    # We use regexpr again on the parent antibody to find the CDR3 location
    cdr3_loc <- regexpr(pep_dt$cdr3_aa, pep_dt$Antibody, fixed = TRUE)
    
    pep_dt[, C_Start := as.integer(cdr3_loc)]
    pep_dt[, C_End   := C_Start + nchar(cdr3_aa) - 1]
    
    # Overlap Logic: (PepStart <= CDR3End) AND (PepEnd >= CDR3Start)
    # Also ensure CDR3 was actually found (C_Start > 0)
    pep_dt[, Is_CDR3_Spanning := (Start <= C_End) & (End >= C_Start) & (C_Start > 0)]
    
    # Cleanup
    pep_dt[is.na(Is_CDR3_Spanning), Is_CDR3_Spanning := FALSE]
    
    # Select Final Columns
    pep_dt <- pep_dt[, ..S1_OUT_COLS]
    
    # Write Output
    out_file <- file.path(out_dir, fn)
    fwrite(pep_dt, out_file, sep = ",", compress = "gzip")
    
    return(invisible(TRUE))
  }
  
  # Filter tasks for Section 1
  processed_s1 <- readLines(logs$s1_proc, warn = FALSE)
  s1_tasks <- task_dt
  s1_tasks[, out_path := file.path(P$digest_csv_dir, filename)]
  
  if (!FORCE_REBUILD_S1) {
    s1_tasks <- s1_tasks[!(filename %in% processed_s1 & file.exists(out_path))]
  }
  
  if (nrow(s1_tasks) > 0) {
    cat("Digesting", nrow(s1_tasks), "files...\n")
    
    # Parallel Execution
    if (.Platform$OS.type == "windows") plan(multisession, workers = N_CORES) else plan(multicore, workers = N_CORES)
    
    res <- future_lapply(seq_len(nrow(s1_tasks)), function(i) {
      fn <- s1_tasks$filename[i]
      fp <- s1_tasks$file_path[i]
      tryCatch({
        digest_worker(fp, fn, ref_dt, P$digest_csv_dir)
        return(list(fn = fn, ok = TRUE))
      }, error = function(e) {
        return(list(fn = fn, ok = FALSE, msg = conditionMessage(e)))
      })
    })
    plan(sequential)
    
    # Logging
    ok_files <- sapply(res, function(x) if(x$ok) x$fn else NA)
    ok_files <- na.omit(ok_files)
    if(length(ok_files)) write(ok_files, logs$s1_proc, append = TRUE)
    cat("S1 Done. OK:", length(ok_files), "\n")
  } else {
    cat("S1 Skipped (All done).\n")
  }
  
  # ===========================================================================
  # SECTION 2: STAGING (CSV -> PARQUET)
  # ===========================================================================
  cat("\n=== SECTION 2: STAGING ===\n")
  
  # Identify Digested CSVs
  dig_files <- list.files(P$digest_csv_dir, full.names = TRUE, pattern = "\\.csv\\.gz$")
  s2_tasks  <- data.table(file_path = dig_files, filename = basename(dig_files))
  
  # Filter tasks
  processed_s2 <- readLines(logs$s2_proc, warn = FALSE)
  if (!FORCE_REBUILD_S2) {
    s2_tasks <- s2_tasks[!(filename %in% processed_s2)]
  }
  
  if (nrow(s2_tasks) > 0) {
    cat("Staging", nrow(s2_tasks), "files to Parquet...\n")
    
    # Define Worker for Staging
    stage_worker <- function(fp, fn, meta_row, stage_root) {
      # Read CSV
      dt <- fread(fp)
      
      # Add Metadata for Partitioning
      dt[, `:=`(
        filename = fn,
        Disease  = as.character(meta_row$Disease),
        BSource  = as.character(meta_row$BSource),
        BType    = as.character(meta_row$BType),
        Isotype  = as.character(meta_row$Isotype)
      )]
      
      # Write Partitioned Parquet
      # Creates structure: /stage_root/file_XXX/Disease=Y/BSource=Z/...
      # We partition by file first to avoid write conflicts
      this_stage_dir <- file.path(stage_root, paste0("file_", fn))
      
      arrow::write_dataset(
        dt, 
        path = this_stage_dir, 
        format = "parquet",
        partitioning = c("Disease", "BSource", "BType", "Isotype"),
        existing_data_behavior = "overwrite"
      )
    }
    
    # Parallel Staging
    if (.Platform$OS.type == "windows") plan(multisession, workers = N_CORES) else plan(multicore, workers = N_CORES)
    
    res2 <- future_lapply(seq_len(nrow(s2_tasks)), function(i) {
      fn <- s2_tasks$filename[i]
      fp <- s2_tasks$file_path[i]
      
      # Get Metadata
      m <- metadata[J(fn)]
      if (nrow(m) == 0) return(list(fn = fn, ok = FALSE, msg = "No Metadata"))
      
      tryCatch({
        stage_worker(fp, fn, m, P$staging_root)
        return(list(fn = fn, ok = TRUE))
      }, error = function(e) {
        return(list(fn = fn, ok = FALSE, msg = conditionMessage(e)))
      })
    })
    plan(sequential)
    
    # Log
    ok_files2 <- sapply(res2, function(x) if(x$ok) x$fn else NA)
    ok_files2 <- na.omit(ok_files2)
    if(length(ok_files2)) write(ok_files2, logs$s2_proc, append = TRUE)
    cat("S2 Done. OK:", length(ok_files2), "\n")
    
  } else {
    cat("S2 Skipped (All done).\n")
  }
  
  # ===========================================================================
  # SECTION 3: COMPACTION (DUCKDB)
  # ===========================================================================
  cat("\n=== SECTION 3: COMPACTION ===\n")
  
  # Check if final exists
  final_exists <- length(list.files(P$db_state1, pattern = "parquet", recursive = TRUE)) > 0
  
  if (final_exists && !FORCE_COMPACT) {
    cat("Final DB exists. Skipping compaction.\n")
  } else {
    cat("Compacting Staging -> Final DB...\n")
    
    con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:")
    on.exit(dbDisconnect(con, shutdown = TRUE))
    
    dbExecute(con, paste0("PRAGMA memory_limit='", DUCKDB_MEMORY_LIMIT, "'"))
    dbExecute(con, paste0("PRAGMA threads=", N_CORES))
    
    # Source pattern: all parquet files in staging
    src_glob <- paste0(P$staging_root, "/**/*.parquet")
    dest_path <- P$db_state1
    
    # SQL Copy
    sql <- sprintf("
      COPY (SELECT * FROM read_parquet('%s'))
      TO '%s'
      (FORMAT PARQUET, PARTITION_BY (Disease, BSource, BType, Isotype), 
       COMPRESSION ZSTD, ROW_GROUP_SIZE %d)
    ", src_glob, dest_path, ROW_GROUP_SIZE)
    
    tryCatch({
      dbExecute(con, sql)
      cat("Compaction Complete.\n")
      
      # Optional: Clean Staging
      # unlink(P$staging_root, recursive = TRUE)
      
    }, error = function(e) {
      cat("Compaction Failed:", conditionMessage(e), "\n")
    })
  }
  
  cat("\nPipeline Finished.\n")
}

main()