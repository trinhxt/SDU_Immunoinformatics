################################################################################
## scripts/02_digestion.R
##
## Part 2: Tryptic Digestion Only (Simplified)
##
## Workflow:
## 1. Reads processed Antibody sequences (from Part 1).
## 2. In-Silico Tryptic Digestion (Cleaver).
## 3. Filters:
##    - Removes short peptides (< 6 AA).
##    - Removes peptides found in UniProt Reference (Background).
## 4. Saves result as CSV.GZ.
##
## Updates:
## - Resumable: Checks output folder and skips existing files.
## - Filtered: Only processes files where Disease != "None".
################################################################################

# ==============================================================================
# Load Configuration
# ==============================================================================
source("scripts/00_config.R")

suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Biostrings)
  library(cleaver)
  library(future)
  library(future.apply)
  library(dplyr)
})

main <- function() {
  
  cat("================================================================\n")
  cat(" PART 2: DIGESTION (Simplified + Disease Filter)\n")
  cat("================================================================\n")
  cat("Input Dir:     ", P$processed_dir, "\n")
  cat("Output Dir:    ", P$digest_csv_dir, "\n")
  cat("Cores:         ", N_CORES, "\n")
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
  
  # Load Metadata (to filter for valid files)
  if (!file.exists(P$metadata_csv)) stop("Metadata CSV missing: ", P$metadata_csv)
  cat("Loading Metadata...\n")
  meta_cols <- c("Filename", "Disease", "BSource", "BType", "Isotype")
  metadata  <- fread(P$metadata_csv, select = meta_cols)
  
  # ---------------------------------------------------------------------------
  # CRITICAL UPDATE: Filter for Disease != "None"
  # ---------------------------------------------------------------------------
  initial_count <- nrow(metadata)
  metadata <- metadata[Disease != "None"]
  filtered_count <- nrow(metadata)
  
  cat(sprintf("Metadata Filtering:\n - Total entries: %d\n - Disease != 'None': %d\n - Removed: %d\n\n", 
              initial_count, filtered_count, initial_count - filtered_count))
  
  if (filtered_count == 0) stop("No files found with Disease != 'None'. Check metadata.")
  
  setkey(metadata, Filename)
  
  # ---------------------------------------------------------------------------
  # 1. Identify and Filter Tasks (Resumability)
  # ---------------------------------------------------------------------------
  
  # List all possible inputs
  all_files <- list.files(P$processed_dir, full.names = TRUE, pattern = "\\.csv\\.gz$")
  task_dt   <- data.table(file_path = all_files, filename = basename(all_files))
  
  # Filter 1: Must be in our filtered metadata (Disease != None)
  task_dt <- task_dt[filename %in% metadata$Filename]
  
  # Filter 2: RESUMABILITY CHECK
  # Check which files already exist in the output directory
  existing_outputs <- list.files(P$digest_csv_dir, pattern = "\\.csv\\.gz$")
  
  # Exclude existing files from the task list
  tasks_to_run <- task_dt[!(filename %in% existing_outputs)]
  setorder(tasks_to_run, filename)
  
  cat("Valid Disease Files found on disk: ", nrow(task_dt), "\n")
  cat("Already Finished:                  ", length(existing_outputs), "\n")
  cat("Remaining Tasks to Run:            ", nrow(tasks_to_run), "\n\n")
  
  if (nrow(tasks_to_run) == 0) {
    cat("All valid files digested. Exiting.\n")
    return(invisible(NULL))
  }
  
  # ===========================================================================
  # SECTION 1: DIGESTION WORKER
  # ===========================================================================
  
  digest_worker <- function(file_path, fn, ref_dt, out_dir) {
    
    # 1. Read Data
    cols_to_read <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr3_aa")
    tab <- arrow::read_csv_arrow(file_path, col_select = cols_to_read)
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
    
    # 3. Filters
    #    Length >= 6 (Mass Spec visibility)
    pep_dt <- pep_dt[nchar(Peptide) >= 6]
    
    #    Reference Filter (Remove common host proteins)
    setkey(pep_dt, Peptide)
    pep_dt <- pep_dt[!ref_dt] 
    
    # 4. Join Metadata (V/D/J/CDR3)
    #    Map metadata back to peptides using the Parent Antibody sequence
    dt_meta <- unique(dt[, .(sequence_alignment_aa, v_call, d_call, j_call, cdr3_aa)])
    pep_dt  <- dt_meta[pep_dt, on = c(sequence_alignment_aa = "Antibody")]
    
    #    Rename back for clarity
    setnames(pep_dt, "sequence_alignment_aa", "Antibody")
    
    # 5. Write Result
    #    Columns: Peptide, Antibody, V, D, J, CDR3 (No Start/End)
    final_cols <- c("Peptide", "Antibody", "v_call", "d_call", "j_call", "cdr3_aa")
    
    # Check if cols exist (safety)
    found_cols <- intersect(final_cols, names(pep_dt))
    pep_dt <- pep_dt[, ..found_cols]
    
    out_file <- file.path(out_dir, fn)
    fwrite(pep_dt, out_file, sep = ",", compress = "gzip")
    
    return(invisible(TRUE))
  }
  
  # ===========================================================================
  # EXECUTION
  # ===========================================================================
  
  cat("Starting digestion on", nrow(tasks_to_run), "files...\n")
  
  # Parallel Execution
  if (.Platform$OS.type == "windows") plan(multisession, workers = N_CORES) else plan(multicore, workers = N_CORES)
  
  # Future Apply Loop
  res <- future_lapply(seq_len(nrow(tasks_to_run)), function(i) {
    fn <- tasks_to_run$filename[i]
    fp <- tasks_to_run$file_path[i]
    
    tryCatch({
      digest_worker(fp, fn, ref_dt, P$digest_csv_dir)
      return(list(fn = fn, ok = TRUE))
    }, error = function(e) {
      return(list(fn = fn, ok = FALSE, msg = conditionMessage(e)))
    })
  })
  
  plan(sequential)
  
  # Reporting
  ok_count <- sum(sapply(res, function(x) x$ok))
  fail_count <- sum(sapply(res, function(x) !x$ok))
  
  cat("\nDigestion Complete.\n")
  cat("Successful: ", ok_count, "\n")
  cat("Failed:     ", fail_count, "\n")
  
  # Log Failures
  failures <- Filter(function(x) !x$ok, res)
  if (length(failures) > 0) {
    fail_dt <- rbindlist(lapply(failures, as.data.table))
    fwrite(fail_dt, log_fail)
    cat("Failures written to:", log_fail, "\n")
  }
}

# Run Main
main()
