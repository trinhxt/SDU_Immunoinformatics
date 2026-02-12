# Computational Settings (Adjust below numbers based on your CPU and RAM)
N_CORES <- 6                   # Number of CPU Cores 
DUCKDB_MEMORY_LIMIT <- "48GB"  # DuckDB memory limit (reduce if you have less than 64GB RAM)

################################################################################
## DBbuild.R
##
## Goal: Master Controller for the OAS Proteomics Pipeline (Cross-Platform)
##       - Installs ALL required dependencies automatically.
##       - Interactive folder selection via svDialogs.
##       - Maps to EXISTING folder structure to allow resuming without rebuilds.
##
## Usage:
##   1. Ensure scripts 01-05 have 'source("scripts/00_config.R")' REMOVED.
##   2. Run: source("DBbuild.R")
##   3. Select the folder containing 'antibody', 'raw', and 'tryptic'.
################################################################################

# ==============================================================================
# 1. Dependency Management (Install Missing Packages)
# ==============================================================================
cat("=== Checking Dependencies ===\n")

# List of all packages used across scripts 01-05
required_packages <- c(
  "svDialogs",    # For interactive dialogs
  "parallel",     # Core R package (usually pre-installed)
  "data.table",   # Fast data manipulation
  "jsonlite",     # JSON parsing (01_download.R)
  "R.utils",      # Utilities (01_download.R)
  "arrow",        # Parquet file handling
  "Biostrings",   # Biological sequences (Bioconductor)
  "cleaver",      # In-silico digestion (Bioconductor)
  "future",       # Parallel processing
  "future.apply", # Parallel apply functions
  "dplyr",        # Data manipulation
  "progressr",    # Progress bars
  "stringi",      # Fast string operations
  "duckdb",       # Database engine
  'DT', 'shiny', 'shinyBS', 'shinyFiles', 'shinyjs', 'shinythemes', 'shinyWidgets' #Shinyapp
)

# Function to check and install standard CRAN packages
install_if_missing <- function(pkgs) {
  missing <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(missing) > 0) {
    cat("Installing missing CRAN packages:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
}

# Function to check and install Bioconductor packages
install_bioc_if_missing <- function(pkgs) {
  missing <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
  if (length(missing) > 0) {
    cat("Installing missing Bioconductor packages:", paste(missing, collapse = ", "), "\n")
    if (!require("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(missing, update = FALSE)
  }
}

# Separate lists
bioc_pkgs <- c("Biostrings", "cleaver")
cran_pkgs <- setdiff(required_packages, bioc_pkgs)

# Install
install_if_missing(cran_pkgs)
install_bioc_if_missing(bioc_pkgs)

# Load essential libraries for this script
suppressPackageStartupMessages({
  library(svDialogs)
  library(parallel)
})

cat("All dependencies are satisfied.\n\n")


# ==============================================================================
# 2. Setup & Interactive Selection
# ==============================================================================
cat("=== OAS Pipeline Setup ===\n")

# Select Data Directory
# We ask the user to select the ROOT folder of their data.
user_selected_dir <- dlg_dir(default = getwd(), 
                             title = "Select your Data Root Folder (at least 200GB free space)")$res

if (length(user_selected_dir) == 0) {
  stop("Operation cancelled: No directory selected.")
}

cat("\nSelected Data Root: ", user_selected_dir, "\n")


# ==============================================================================
# 3. Path Configuration (Mapped to Your Tree)
# ==============================================================================

# Current Project Root (where this script and the 'scripts/' folder are)
PROJECT_ROOT <- getwd()

# Data Root (The folder you selected)
DATA_ROOT <- user_selected_dir

# Helper for OS-agnostic path joining
pjoin <- function(...) file.path(..., fsep = .Platform$file.sep)

# GLOBAL CONFIGURATION LIST 'P'
P <- list(
  # --- Inputs (Code Repository) ---
  shell_script    = pjoin(PROJECT_ROOT, "data", "bulk_download_human_unpaired_heavychain.sh"),
  ref_rdata       = pjoin(PROJECT_ROOT, "data", "UniProtNCBI_Tryptic.RData"),
  
  # --- Data Directories (Mapped to User Selection) ---
  root_dir        = DATA_ROOT,
  
  # 1. Raw Downloads
  raw_dir         = pjoin(DATA_ROOT, "raw"),
  
  # 2. Processed Antibody CSVs
  processed_dir   = pjoin(DATA_ROOT, "antibody"),
  
  # 3. Metadata (Assumed at the root of your data folder)
  metadata_csv    = pjoin(DATA_ROOT, "OAS_metadata.csv"),
  
  # 4. Tryptic Data Structure
  tryptic_dir     = pjoin(DATA_ROOT, "tryptic"),
  digest_csv_dir  = pjoin(DATA_ROOT, "tryptic", "digest_csv_gz"),
  log_dir         = pjoin(DATA_ROOT, "tryptic", "log_files"),
  parquet_dir     = pjoin(DATA_ROOT, "tryptic", "parquet_db"),
  
  # 5. Database Stages
  db_stage1       = pjoin(DATA_ROOT, "tryptic", "parquet_db", "db_stage1"),
  
  # 6. Temporary Work Dir (Safe to create/delete)
  work_dir        = pjoin(DATA_ROOT, "temp_work")
)

# ==============================================================================
# 4. Validation & Initialization
# ==============================================================================

cat("\nVerifying Directory Structure...\n")
dirs_to_check <- c(
  P$raw_dir, 
  P$processed_dir, 
  P$tryptic_dir,
  P$digest_csv_dir,
  P$log_dir,
  P$parquet_dir,
  P$work_dir
)

# Auto-create folders if they don't exist
for (d in dirs_to_check) {
  if (!dir.exists(d)) {
    cat("  [Creating] ", d, "\n")
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  } else {
    cat("  [Found]    ", d, "\n")
  }
}

cat("\nConfiguration Loaded.\n")
cat("  - Cores:     ", N_CORES, "\n")
cat("  - Memory:    ", DUCKDB_MEMORY_LIMIT, "\n\n")

# ==============================================================================
# 5. Pipeline Execution
# ==============================================================================

pipeline_scripts <- list(
  "Step 1: Download & Preprocess OAS"           = "scripts/01_download.R",
  "Step 2: Tryptic Digestion"                   = "scripts/02_digestion.R",
  "Step 3: Build Fact Table (db_stage1)"        = "scripts/03_build_fact_table.R",
  "Step 4: Build Dimension Table (db_stage2)"   = "scripts/04_build_dimension_table.R",
  "Step 5: Combine & Finalize (db_stage3)"      = "scripts/05_combine_db.R"
)

run_step <- function(step_name, script_path) {
  cat("============================================================\n")
  cat(">>> Running:", step_name, "\n")
  cat("Script     :", script_path, "\n")
  cat("============================================================\n")
  
  if (!file.exists(script_path)) {
    stop(paste("CRITICAL ERROR: Script not found at:", script_path))
  }
  
  t_start <- Sys.time()
  
  # Source with local=FALSE so scripts can see the global 'P', 'N_CORES', etc.
  tryCatch({
    source(script_path, local = FALSE)
  }, error = function(e) {
    cat("\n❌ FAILURE in", step_name, "\n")
    stop(conditionMessage(e))
  })
  
  t_end <- Sys.time()
  cat("\n✅ Finished:", step_name, "\n")
  cat("Time:", round(difftime(t_end, t_start, units="mins"), 2), "mins.\n\n")
}

# Execute Sequential Steps
for (nm in names(pipeline_scripts)) {
  run_step(nm, pipeline_scripts[[nm]])
}

cat("\n============================================================\n")
cat("🎉 PIPELINE COMPLETED SUCCESSFULLY 🎉\n")
cat("Data Location: ", DATA_ROOT, "\n")
cat("============================================================\n")