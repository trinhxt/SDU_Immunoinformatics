################################################################################
## 00_config.R
##
## Single configuration file for the OAS Proteomics Workflow.
## Source this at the top of your scripts: source("scripts/00_config.R")
################################################################################

# ==============================================================================
# 1. Project Directories (Customize these)
# ==============================================================================
# Where the raw OAS shell script is located (Download this from OAS website)
DOWNLOAD_SHELL_SCRIPT <- "data/bulk_download_human_unpaired_heavychain.sh"

# Where to store the raw downloaded files (can be deleted after processing)
RAW_DATA_DIR <- "D:/OAS/unpaired/raw"

# Where to save the cleaned, processed CSV.GZ files
PROCESSED_DATA_DIR <- "D:/OAS/unpaired/antibody"

# Where to save the tryptic peptide data (CSV.GZ file)
TRYPTIC_DATA_DIR <- "D:/OAS/unpaired/tryptic"

# Where to save the metadata CSV
METADATA_FILE <- "D:/OAS/unpaired/OAS_metadata.csv"

# Where to build the final Parquet database
PARQUET_DB_DIR <- "D:/OAS/unpaired/tryptic/parquet_db"

# Where to store temporary files
WORK_DIR <- "E:/OAS_temp"

# Reference file for filtering (Uniprot Tryptic Peptides)
REF_PEPTIDE_RDATA <- "data/UniProtNCBI_Tryptic.RData"


# ==============================================================================
# 2. Computational Settings
# ==============================================================================
# Number of cores for parallel processing (downloading/digesting)
N_CORES <- 6L

# DuckDB settings for database compaction
DUCKDB_MEMORY_LIMIT <- "48GB"
ROW_GROUP_SIZE <- 1000000L


# ==============================================================================
# 3. Helper Functions (Path Normalization)
# ==============================================================================
# Ensure directories exist
for (d in c(RAW_DATA_DIR, PROCESSED_DATA_DIR, PARQUET_DB_DIR, WORK_DIR, dirname(METADATA_FILE))) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Normalize paths for Windows/Linux compatibility
norm_path <- function(p) normalizePath(p, winslash = "/", mustWork = FALSE)
pjoin <- function(...) norm_path(paste(..., sep = "/"))

P <- list(
  shell_script    = norm_path(DOWNLOAD_SHELL_SCRIPT),
  raw_dir         = norm_path(RAW_DATA_DIR),
  tryptic_dir     = norm_path(TRYPTIC_DATA_DIR),
  processed_dir   = norm_path(PROCESSED_DATA_DIR),
  metadata_csv    = norm_path(METADATA_FILE),
  parquet_dir     = norm_path(PARQUET_DB_DIR),
  work_dir        = norm_path(WORK_DIR),
  
  # Part 2 Specific Paths (Derived from WORK_DIR and PARQUET_DB_DIR)
  work_root       = norm_path(WORK_DIR),
  digest_csv_dir  = pjoin(TRYPTIC_DATA_DIR, "digest_csv_gz"),
  log_dir         = pjoin(TRYPTIC_DATA_DIR, "log_files"),
  staging_root    = pjoin(PARQUET_DB_DIR, "_staging"),
  db_stage1       = pjoin(PARQUET_DB_DIR, "db_stage1"),
  
  ref_rdata       = norm_path(REF_PEPTIDE_RDATA)
)

message("Configuration loaded.")
message(" - Processed Input: ", P$processed_dir)
message(" - Final Output DB: ", P$db_stage1)