################################################################################
## scripts/01_part1_download.R
##
## Part1: Download + preprocess OAS CSV.GZ files
##
## Outputs
## -------
## 1) Processed OAS files (same filenames as OAS download):
##      P$oas_processed_dir/*.csv.gz
##    Columns kept:
##      sequence_alignment_aa, v_call, d_call, j_call, cdr1_aa, cdr2_aa, cdr3_aa
##
## 2) Metadata table:
##      P$metadata_csv  (CSV)
##
## Fast SKIP rule (set-based, no per-file checking)
## ------------------------------------------------
## If ALL are true:
##   - P$metadata_csv exists and has Filename
##   - Filename set matches ALL filenames in bulk_download.sh
##   - processed_dir contains ALL expected *.csv.gz and no extras
## Then we SKIP the whole script (return) and do NOT download/parse anything.
##
## Notes
## -----
## - Uses Windows-safe download.file(mode="wb").
## - Raw staging folder uses P$oas_raw_dir if provided, else tempdir().
## - No quit() is used (safe for source()).
################################################################################


# ==============================================================================
# Load config + path helpers
# ==============================================================================
source("R/helpers_paths.R")
cfg <- load_config()
P   <- get_paths(cfg)

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
  library(dplyr)
  library(jsonlite)
})

main <- function() {
  
  # ----------------------------------------------------------------------------
  # Helpers
  # ----------------------------------------------------------------------------
  read_bulk_urls <- function(sh_path) {
    if (!file.exists(sh_path)) stop("bulk download script not found: ", sh_path)
    lines <- readLines(sh_path, warn = FALSE)
    lines <- trimws(lines)
    lines <- lines[nzchar(lines)]
    lines <- lines[grepl("^wget\\b", lines)]
    urls <- sub("^wget\\s+", "", lines)
    urls <- trimws(urls)
    urls <- urls[nzchar(urls)]
    urls
  }
  
  read_metadata_first_line <- function(gz_path) {
    first_line <- readLines(gzfile(gz_path), n = 1, warn = FALSE)
    first_line <- gsub("\\bNaN\\b", "null", first_line)
    clean_json <- gsub('^"|"$', "", first_line)
    clean_json <- gsub('""', '"', clean_json)
    jsonlite::fromJSON(clean_json)
  }
  
  fast_skip_check <- function(bulk_sh, metadata_csv, processed_dir) {
    
    urls <- read_bulk_urls(bulk_sh)
    expected <- sort(unique(basename(urls)))
    
    if (!file.exists(metadata_csv)) {
      return(list(skip = FALSE, reason = "metadata_missing", expected = expected))
    }
    md <- tryCatch(fread(metadata_csv, showProgress = FALSE), error = function(e) NULL)
    if (is.null(md) || !("Filename" %in% names(md))) {
      return(list(skip = FALSE, reason = "metadata_unreadable_or_no_Filename", expected = expected))
    }
    meta_files <- sort(unique(as.character(md$Filename)))
    
    if (!dir.exists(processed_dir)) {
      return(list(skip = FALSE, reason = "processed_dir_missing", expected = expected, meta_files = meta_files))
    }
    proc_files <- sort(unique(basename(list.files(processed_dir, pattern = "\\.csv\\.gz$", full.names = TRUE))))
    
    miss_in_meta <- setdiff(expected, meta_files)
    miss_in_dir  <- setdiff(expected, proc_files)
    extra_meta   <- setdiff(meta_files, expected)
    extra_dir    <- setdiff(proc_files, expected)
    
    ok <- (length(miss_in_meta) == 0 &&
             length(miss_in_dir)  == 0 &&
             length(extra_meta)   == 0 &&
             length(extra_dir)    == 0)
    
    list(
      skip = ok,
      reason = if (ok) "complete" else "incomplete",
      expected_n = length(expected),
      meta_n     = length(meta_files),
      dir_n      = length(proc_files),
      expected   = expected,
      meta_files = meta_files,
      proc_files = proc_files,
      miss_in_meta = miss_in_meta,
      miss_in_dir  = miss_in_dir,
      extra_meta   = extra_meta,
      extra_dir    = extra_dir
    )
  }
  
  process_one_url <- function(url, i, n, raw_dir, out_dir) {
    file_name <- basename(url)
    cat(sprintf("[Part1] (%d/%d) %s\n", i, n, file_name))
    
    dest_raw <- file.path(raw_dir, file_name)
    
    # Download
    if (!file.exists(dest_raw)) {
      ok_dl <- TRUE
      tryCatch({
        download.file(url, destfile = dest_raw, mode = "wb", quiet = TRUE)
      }, error = function(e) {
        ok_dl <<- FALSE
        cat("  - Download failed: ", conditionMessage(e), "\n", sep = "")
      })
      if (!ok_dl || !file.exists(dest_raw)) {
        return(list(ok = FALSE, filename = file_name, meta = NULL))
      }
    } else {
      cat("  - Raw exists, skip download\n")
    }
    
    out_file <- file.path(out_dir, file_name)
    tmp_out  <- paste0(out_file, ".tmp_", Sys.getpid())
    
    # Already processed
    if (file.exists(out_file)) {
      cat("  - Processed exists, skip processing\n")
      meta <- tryCatch(read_metadata_first_line(dest_raw), error = function(e) NULL)
      meta_dt <- if (is.null(meta)) data.table(Filename = file_name) else as.data.table(t(meta))
      meta_dt[, Filename := file_name]
      return(list(ok = TRUE, filename = file_name, meta = meta_dt))
    }
    
    # Metadata
    meta <- tryCatch(read_metadata_first_line(dest_raw), error = function(e) {
      cat("  - Metadata parse failed: ", conditionMessage(e), "\n", sep = "")
      NULL
    })
    meta_dt <- if (is.null(meta)) data.table() else as.data.table(t(meta))
    meta_dt[, Filename := file_name]
    
    # Data extraction
    ok_proc <- TRUE
    tryCatch({
      tab <- arrow::read_csv_arrow(dest_raw, skip = 1)
      
      keep <- c("sequence_alignment_aa", "v_call", "d_call", "j_call", "cdr1_aa", "cdr2_aa", "cdr3_aa")
      present <- intersect(keep, names(tab))
      
      if (!("sequence_alignment_aa" %in% present)) {
        stop("Missing required column sequence_alignment_aa in: ", file_name)
      }
      
      dt <- as.data.table(tab[, present, drop = FALSE])
      for (cc in setdiff(keep, names(dt))) dt[, (cc) := NA_character_]
      
      dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
      fwrite(dt, file = tmp_out, compress = "gzip")
      file.rename(tmp_out, out_file)
      
    }, error = function(e) {
      ok_proc <<- FALSE
      cat("  - Processing failed: ", conditionMessage(e), "\n", sep = "")
      if (file.exists(tmp_out)) file.remove(tmp_out)
    })
    
    if (ok_proc) {
      try(file.remove(dest_raw), silent = TRUE)
    }
    
    list(ok = ok_proc, filename = file_name, meta = meta_dt)
  }
  
  # ----------------------------------------------------------------------------
  # Paths
  # ----------------------------------------------------------------------------
  shell_script_file <- P$oas_download_script
  processed_dir     <- P$oas_processed_dir
  metadata_csv      <- P$metadata_csv
  
  raw_dir <- if (nzchar(P$oas_raw_dir)) P$oas_raw_dir else tempdir()
  dir.create(raw_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(metadata_csv), showWarnings = FALSE, recursive = TRUE)
  
  cat("Bulk download script :", shell_script_file, "\n")
  cat("Raw staging dir      :", raw_dir, "\n")
  cat("Processed dir        :", processed_dir, "\n")
  cat("Metadata CSV         :", metadata_csv, "\n\n")
  
  # ----------------------------------------------------------------------------
  # Fast skip check (set compare)
  # ----------------------------------------------------------------------------
  chk <- fast_skip_check(shell_script_file, metadata_csv, processed_dir)
  
  if (isTRUE(chk$skip)) {
    cat("Part1 SKIP: metadata + processed dir match bulk_download.sh (", chk$expected_n, " files).\n", sep = "")
    return(invisible(TRUE))
  }
  
  cat("Part1 NOT complete -> will process.\n")
  cat("Expected:", chk$expected_n,
      " | Metadata:", if (!is.null(chk$meta_n)) chk$meta_n else NA,
      " | Processed dir:", if (!is.null(chk$dir_n)) chk$dir_n else NA,
      "\n\n")
  
  # ----------------------------------------------------------------------------
  # Process all URLs (sequential)
  # ----------------------------------------------------------------------------
  urls <- read_bulk_urls(shell_script_file)
  if (!length(urls)) stop("No URLs found in: ", shell_script_file)
  
  n <- length(urls)
  results <- vector("list", n)
  for (i in seq_along(urls)) {
    results[[i]] <- process_one_url(urls[[i]], i, n, raw_dir = raw_dir, out_dir = processed_dir)
  }
  
  res_dt <- rbindlist(lapply(results, function(x) data.table(
    Filename = x$filename,
    ok = isTRUE(x$ok)
  )), fill = TRUE)
  
  cat(sprintf("\nFinished Part1. OK: %d | Failed: %d\n", sum(res_dt$ok), sum(!res_dt$ok)))
  
  # ----------------------------------------------------------------------------
  # Build + write metadata (fresh)
  # ----------------------------------------------------------------------------
  meta_list <- lapply(results, function(x) x$meta)
  meta_list <- meta_list[!vapply(meta_list, is.null, logical(1))]
  metadata_dt <- rbindlist(meta_list, fill = TRUE)
  
  if (!("Filename" %in% names(metadata_dt))) metadata_dt[, Filename := NA_character_]
  
  if (all(c("Author", "Subject") %in% names(metadata_dt))) {
    metadata_dt[, Patient := paste(Author, Subject, sep = "_")]
    metadata_dt[, Patient := paste0("P", as.numeric(factor(Patient)))]
  }
  
  metadata_dt[] <- lapply(metadata_dt, function(col) {
    if (is.list(col)) sapply(col, toString) else col
  })
  
  fwrite(metadata_dt, metadata_csv)
  cat("Wrote metadata CSV:\n  ", metadata_csv, "\n", sep = "")
  
  invisible(TRUE)
}

main()
