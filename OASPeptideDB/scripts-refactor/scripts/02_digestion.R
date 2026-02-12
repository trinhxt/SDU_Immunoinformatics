################################################################################
## scripts/02_digestion.R
##
## Part 2: Tryptic Digestion, Filtering, and Feature Engineering
##
## Workflow:
## 1. Reads processed Antibody sequences.
## 2. In-Silico Tryptic Digestion (Cleaver).
## 3. Filters:
##    - Removes short peptides (< 6 AA).
##    - Removes peptides found in UniProt Reference (Global Background).
##
## 4. BRANCHED LOGIC:
##    A. IF HEALTHY (Disease == "None"):
##       - Goal: Build a background exclusion list.
##       - Filter: KEEP ALL peptides (Ignore CDR3 overlap). 
##         *Reason: Disease CDR3s can mimic Healthy Framework regions.*
##       - Output: Saves ONLY the 'Peptide' column (distinct).
##
##    B. IF DISEASE:
##       - Goal: Find candidate biomarkers.
##       - Feature Eng: Calculate CDR3 Overlap.
##       - Filter: KEEP ONLY peptides overlapping CDR3 (CDR3_spanning_count > 0).
##       - Output: Saves FULL metadata (Patient, Disease, Counts, etc.).
##
## 5. Saves result as Flat PARQUET Files (One per input).
##
## Updates:
## - Logic: Split processing for Healthy vs. Disease.
## - Format: Parquet with ZSTD compression.
## - UI: Added Progress Bar (progressr).
################################################################################

if (!exists("P")) stop("Configuration 'P' not found. Please run this via DBbuild.R")

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Biostrings)
  library(cleaver)
  library(future)
  library(future.apply)
  library(dplyr)
  library(progressr) # Added for progress bar
  library(stringi)   # Added for fast vectorized string operations
})

# Enable global progress handlers
handlers(global = TRUE)
handlers("txtprogressbar") 

main <- function() {
  
  cat("================================================================\n")
  cat(" PART 2: DIGESTION (Branch: Healthy vs Disease) \n")
  cat("================================================================\n")
  cat("Input Dir:     ", P$processed_dir, "\n")
  cat("Output Dir:    ", P$digest_csv_dir, "\n")
  cat("Cores:         ", N_CORES, "\n")
  cat("Strategy:      Healthy -> Save Peptide Only\n")
  cat("               Disease -> Filter CDR3 > 0 -> Save Full Meta\n")
  cat("================================================================\n\n")
  
  # ---------------------------------------------------------------------------
  # 0. Setup & Validation
  # ---------------------------------------------------------------------------
  
  # Ensure output directory exists
  dir.create(P$digest_csv_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(P$log_dir,        showWarnings = FALSE, recursive = TRUE)
  
  # Log file for failures
  log_fail <- file.path(P$log_dir, "digestion_failures.txt")
  
  # Load Reference Data (UniProt Tryptic Peptides)
  if (!file.exists(P$ref_rdata)) stop("Reference RData missing: ", P$ref_rdata)
  cat("Loading Reference Proteome Filter...\n")
  load(P$ref_rdata) # Expects 'UniProtNCBI_Tryptic'
  ref_dt <- unique(data.table(Peptide = as.character(UniProtNCBI_Tryptic)))
  setkey(ref_dt, Peptide)
  rm(UniProtNCBI_Tryptic); gc()
  
  # Load Metadata 
  if (!file.exists(P$metadata_csv)) stop("Metadata CSV missing: ", P$metadata_csv)
  cat("Loading Metadata...\n")
  
  meta_cols <- c("Filename", "Patient", "Disease", "BSource", "BType", "Isotype")
  metadata  <- fread(P$metadata_csv, select = meta_cols)
  setkey(metadata, Filename)
  
  # ---------------------------------------------------------------------------
  # 1. Identify and Filter Tasks (Resumability)
  # ---------------------------------------------------------------------------
  
  # List all possible inputs 
  all_files <- list.files(P$processed_dir, full.names = TRUE, pattern = "\\.csv\\.gz$")
  task_dt   <- data.table(file_path = all_files, filename = basename(all_files))
  
  # Filter 1: Must be in our metadata
  task_dt <- task_dt[filename %in% metadata$Filename]
  
  # Filter 2: RESUMABILITY CHECK
  task_dt[, target_name := gsub("\\.csv\\.gz$", ".parquet", filename)]
  task_dt[, target_path := file.path(P$digest_csv_dir, target_name)]
  
  # Exclude existing files from the task list
  tasks_to_run <- task_dt[!file.exists(target_path)]
  setorder(tasks_to_run, filename)
  
  cat("Total Valid Files:             ", nrow(task_dt), "\n")
  cat("Already Finished:              ", nrow(task_dt) - nrow(tasks_to_run), "\n")
  cat("Remaining Tasks to Run:        ", nrow(tasks_to_run), "\n\n")
  
  if (nrow(tasks_to_run) == 0) {
    cat("All valid files digested. Exiting.\n")
    return(invisible(NULL))
  }
  
  # ===========================================================================
  # SECTION 1: DIGESTION WORKER
  # ===========================================================================
  
  digest_worker <- function(file_path, fn, ref_dt, meta_row, out_dir) {
    
    # Check Logic Branch
    is_healthy <- (meta_row$Disease == "None")
    
    # 1. Read Data
    #    We need Sequence Alignment for everyone.
    #    We only need CDR3 and V/D/J for Disease analysis. 
    #    But reading them all is usually fast enough.
    cols_to_read <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr3_aa")
    tab <- arrow::read_csv_arrow(file_path, col_select = all_of(cols_to_read))
    dt  <- as.data.table(tab)
    
    # Basic Clean
    dt <- dt[!is.na(sequence_alignment_aa) & nzchar(sequence_alignment_aa)]
    if (nrow(dt) == 0) return(invisible(TRUE))
    
    # 2. Digestion (Cleaver)
    #    Digest unique sequences only
    seq_vec <- unique(dt$sequence_alignment_aa)
    aa      <- Biostrings::AAStringSet(stats::setNames(seq_vec, seq_vec))
    dig     <- cleaver::cleave(aa, enzym = "trypsin", missedCleavages = 0:2, unique = TRUE)
    
    # Expand results
    lens    <- lengths(dig)
    pep_seq <- as.character(unlist(dig, use.names = FALSE))
    pep_ab  <- rep.int(names(dig), lens) # Parent Sequence
    
    pep_dt <- data.table(
      Peptide  = pep_seq,
      Antibody = pep_ab
    )
    
    # 3. Common Filters
    #    Length >= 6 (Mass Spec visibility)
    pep_dt <- pep_dt[nchar(Peptide) >= 6]
    
    #    Reference Filter (UniProt Background)
    setkey(pep_dt, Peptide)
    pep_dt <- pep_dt[!ref_dt] 
    
    if (nrow(pep_dt) == 0) return(invisible(TRUE))
    
    # -------------------------------------------------------------------------
    # 4. BRANCH: HEALTHY PROCESSING
    # -------------------------------------------------------------------------
    if (is_healthy) {
      
      # Logic: 
      # - We keep ALL peptides (CDR3 + Framework) to serve as a comprehensive background.
      # - We DO NOT calculate overlap (saves compute).
      # - We save ONLY the Peptide column (saves huge disk space).
      
      # Unique Peptides only for this file (we don't need counts for background)
      final_dt <- unique(pep_dt[, .(Peptide)])
      
      # Write result immediately
      out_fn   <- gsub("\\.csv\\.gz$", ".parquet", fn)
      out_file <- file.path(out_dir, out_fn)
      
      arrow::write_parquet(final_dt, out_file, compression = "zstd", compression_level = 15)
      return(invisible(TRUE))
    }
    
    # -------------------------------------------------------------------------
    # 5. BRANCH: DISEASE PROCESSING
    # -------------------------------------------------------------------------
    
    # A. Join Metadata (V/D/J/CDR3)
    dt_meta <- unique(dt[, .(sequence_alignment_aa, v_call, d_call, j_call, cdr3_aa)])
    pep_dt  <- dt_meta[pep_dt, on = c(sequence_alignment_aa = "Antibody")]
    setnames(pep_dt, "sequence_alignment_aa", "Antibody")
    
    # B. Calculate CDR3 Overlap (Vectorized)
    p_loc <- stringi::stri_locate_first_fixed(pep_dt$Antibody, pep_dt$Peptide)
    c_loc <- stringi::stri_locate_first_fixed(pep_dt$Antibody, pep_dt$cdr3_aa)
    
    p_starts <- p_loc[,1]; p_ends <- p_loc[,2]
    c_starts <- c_loc[,1]; c_ends <- c_loc[,2]
    
    invalid <- is.na(p_starts) | is.na(c_starts)
    overlap <- pmax(0L, pmin(p_ends, c_ends) - pmax(p_starts, c_starts) + 1L)
    overlap[invalid] <- 0L
    
    pep_dt[, CDR3_spanning_count := overlap]
    
    # CDR3 Percentage
    cdr3_len <- nchar(pep_dt$cdr3_aa)
    pep_dt[, CDR3_spanning_pct := ifelse(cdr3_len > 0, overlap / cdr3_len, 0)]
    
    # C. FILTER: Strict CDR3 Overlap for Disease
    #    We only want peptides that actually represent the hypervariable region.
    pep_dt <- pep_dt[CDR3_spanning_count > 0]
    
    if (nrow(pep_dt) == 0) return(invisible(TRUE))
    
    # D. Attach Clinical Metadata
    pep_dt[, `:=`(
      Patient = as.character(meta_row$Patient),
      Disease = as.character(meta_row$Disease),
      BSource = as.character(meta_row$BSource),
      BType   = as.character(meta_row$BType),
      Isotype = as.character(meta_row$Isotype),
      filename = fn 
    )]
    
    # E. Select Final Columns (Full Schema)
    final_cols <- c(
      "Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa",
      "CDR3_spanning_count", "CDR3_spanning_pct",
      "filename", "Patient", "Disease", "BSource", "BType", "Isotype"
    )
    final_dt <- pep_dt[, ..final_cols]
    
    # F. Write Result
    out_fn   <- gsub("\\.csv\\.gz$", ".parquet", fn)
    out_file <- file.path(out_dir, out_fn)
    
    arrow::write_parquet(final_dt, out_file, compression = "zstd", compression_level = 15)
    
    return(invisible(TRUE))
  }
  
  # ===========================================================================
  # EXECUTION
  # ===========================================================================
  
  cat("Starting digestion on", nrow(tasks_to_run), "files...\n")
  
  if (.Platform$OS.type == "windows") plan(multisession, workers = N_CORES) else plan(multicore, workers = N_CORES)
  
  with_progress({
    p <- progressor(along = 1:nrow(tasks_to_run))
    
    res <- future_lapply(seq_len(nrow(tasks_to_run)), function(i) {
      fn <- tasks_to_run$filename[i]
      fp <- tasks_to_run$file_path[i]
      
      # Get Metadata
      m_row <- metadata[J(fn)]
      
      if (nrow(m_row) == 0) {
        return(list(fn = fn, ok = FALSE, msg = "Metadata not found"))
      }
      
      tryCatch({
        digest_worker(fp, fn, ref_dt, m_row, P$digest_csv_dir)
        p(sprintf("Done: %s", fn)) 
        return(list(fn = fn, ok = TRUE))
      }, error = function(e) {
        p(sprintf("Failed: %s", fn)) 
        return(list(fn = fn, ok = FALSE, msg = conditionMessage(e)))
      })
    })
  })
  
  plan(sequential)
  
  # Reporting
  ok_count <- sum(sapply(res, function(x) x$ok))
  fail_count <- sum(sapply(res, function(x) !x$ok))
  
  cat("\nDigestion Complete.\n")
  cat("Successful Files: ", ok_count, "\n")
  cat("Failed Files:     ", fail_count, "\n")
  
  failures <- Filter(function(x) !x$ok, res)
  if (length(failures) > 0) {
    fail_dt <- rbindlist(lapply(failures, as.data.table))
    fwrite(fail_dt, log_fail)
    cat("Failures written to:", log_fail, "\n")
  }
}

# Run Main
main()
