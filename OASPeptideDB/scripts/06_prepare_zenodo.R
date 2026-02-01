################################################################################
## scripts/06_prepare_zenodo.R
##
## Goal: Package the databases for Zenodo publication.
##
## Workflow:
## 1. Asks user for an Output Directory (Interactive).
## 2. Creates a release folder: 'OAS_Tryptic_Peptide_Database' inside that dir.
## 3. Copies & Renames databases:
##    - db_stage1 -> 'data/observations' (Fact Table)
##    - db_stage2 -> 'data/peptide_metrics' (Dimension Table)
## 4. Generates documentation (README.txt, quick_start.R).
################################################################################

source("scripts/00_config.R")

# fs is faster/safer for directory copying than base R
if (!require("fs")) install.packages("fs")
library(fs) 

# Helper to ask for directory interactively
get_user_dir <- function(default_path) {
  out_path <- NA
  
  if (interactive()) {
    cat("\n--------------------------------------------------------\n")
    cat(" ACTION REQUIRED: Select destination folder for Release.\n")
    cat(" If cancelled, will default to: ", default_path, "\n")
    cat("--------------------------------------------------------\n")
    flush.console()
    
    if (.Platform$OS.type == "windows") {
      out_path <- utils::choose.dir(default = default_path, caption = "Select Folder for Zenodo Release")
    } else if (requireNamespace("tcltk", quietly = TRUE)) {
      out_path <- tcltk::tk_choose.dir(default = default_path, caption = "Select Folder for Zenodo Release")
    }
  }
  
  if (is.na(out_path) || !nzchar(out_path)) {
    cat("No folder selected (or non-interactive). Using default: ", default_path, "\n")
    return(default_path)
  }
  return(out_path)
}

main <- function() {
  
  # ---------------------------------------------------------------------------
  # 1. Setup Input/Output
  # ---------------------------------------------------------------------------
  INPUT_FACT  <- file.path(P$parquet_dir, "db_stage1")
  INPUT_DIM   <- file.path(P$parquet_dir, "db_stage2")
  
  # Validate Inputs
  if (!dir.exists(INPUT_FACT)) stop("db_stage1 missing! Run scripts/03 first.")
  if (!dir.exists(INPUT_DIM))  stop("db_stage2 missing! Run scripts/04 first.")
  
  # Select Output Directory
  base_output_dir <- get_user_dir(P$work_dir)
  RELEASE_DIR     <- file.path(base_output_dir, "OAS_Tryptic_Peptide_Database")
  
  cat("=== STEP 6: PREPARING ZENODO RELEASE ===\n")
  cat("Source Fact DB:      ", INPUT_FACT, "\n")
  cat("Source Dimension DB: ", INPUT_DIM, "\n")
  cat("Target Release Dir:  ", RELEASE_DIR, "\n\n")
  
  # ---------------------------------------------------------------------------
  # 2. Create Directory Structure
  # ---------------------------------------------------------------------------
  if (dir.exists(RELEASE_DIR)) {
    cat("Removing previous release build at target...\n")
    dir_delete(RELEASE_DIR)
  }
  dir_create(RELEASE_DIR)
  dir_create(file.path(RELEASE_DIR, "data"))
  
  # ---------------------------------------------------------------------------
  # 3. Copy & Rename Databases
  # ---------------------------------------------------------------------------
  cat("Copying Observations DB (Fact Table)... (This may take time)\n")
  # Rename to 'observations' for clarity
  dir_copy(INPUT_FACT, file.path(RELEASE_DIR, "data", "observations"))
  
  cat("Copying Metrics DB (Dimension Table)... (This may take time)\n")
  # Rename to 'peptide_metrics' for clarity
  dir_copy(INPUT_DIM, file.path(RELEASE_DIR, "data", "peptide_metrics"))
  
  # ---------------------------------------------------------------------------
  # 4. Generate User Documentation (README)
  # ---------------------------------------------------------------------------
  cat("Generating Documentation...\n")
  
  readme_text <- "
OAS Human Unpaired Heavy Chain Tryptic Peptide Database
=======================================================

This dataset contains tryptic peptides derived from the OAS (Observed Antibody Space) database.
It is structured as a 'Star Schema' with two main components: Observations (Fact) and Metrics (Dimension).

CONTENTS
--------

1. data/observations/ (Parquet Dataset)
   ------------------------------------
   * Description: Contains every observation of a tryptic peptide in the source data.
   * Structure: Partitioned by Disease, BSource, BType, Isotype.
   * Columns:
     - Peptide: Amino acid sequence.
     - Antibody: Full antibody variable region sequence.
     - cdr3_aa: CDR3 region sequence.
     - CDR3_spanning_count: Integer (Number of overlapping amino acids between Peptide and CDR3).
     - Patient: Patient ID (derived from source metadata).
     - Disease, BSource, BType, Isotype: Clinical metadata.
     - v_call, d_call, j_call: Germline gene annotations.
     - filename: Source file name.

2. data/peptide_metrics/ (Parquet Dataset)
   ---------------------------------------
   * Description: Contains aggregated statistics for every unique peptide.
   * Structure: Partitioned by Hash(Peptide) for optimized joining.
   * Columns:
     - Peptide: Amino acid sequence (Primary Key).
     - N_Patients: Number of unique patients carrying this peptide.
     - N_Diseases: Number of distinct diseases this peptide appears in.
     - N_Antibodies: Number of unique antibody sequences containing this peptide.
     - partition_id: Optimized integer ID for joins.

HOW TO USE (Query Strategy)
---------------------------
To analyze this data efficiently, perform a JOIN between the two tables on the 'Peptide' column.
- Use 'observations' to filter by clinical criteria (e.g., Disease = 'Lupus').
- Use 'peptide_metrics' to filter by statistical criteria (e.g., N_Patients > 10).

See 'quick_start.R' for a working code example.
"
  writeLines(readme_text, file.path(RELEASE_DIR, "README.txt"))
  
  # ---------------------------------------------------------------------------
  # 5. Generate Quick Start Script
  # ---------------------------------------------------------------------------
  cat("Generating Quick Start Script...\n")
  
  quick_start_code <- '
# Quick Start: Querying the OAS Peptide Database
# Required: install.packages("duckdb")

library(duckdb)

# 1. Connect (In-Memory is fine for querying)
con <- dbConnect(duckdb())

# 2. Register the Parquet Files as Views (Virtual Tables)
#    This is "Zero-Copy" - it reads metadata instantly without loading data.
dbExecute(con, "CREATE VIEW observations AS SELECT * FROM read_parquet(\'data/observations/**/*.parquet\')")
dbExecute(con, "CREATE VIEW metrics      AS SELECT * FROM read_parquet(\'data/peptide_metrics/**/*.parquet\')")

# 3. Example Query: The "Single Database" Experience
#    Goal: Find peptides in SARS-CoV-2 patients that are shared by at least 5 patients.

sql <- "
  SELECT 
    obs.Peptide,
    obs.Antibody,
    obs.Patient,
    obs.Disease,
    met.N_Patients,
    met.N_Diseases
  FROM observations obs
  JOIN metrics met ON obs.Peptide = met.Peptide
  WHERE 
    obs.Disease = \'SARS-CoV-2\'      -- Filter Raw Data (Uses Partition Pruning)
    AND met.N_Patients >= 5           -- Filter Stats (Uses Hash Join)
  LIMIT 20
"

cat("Running Query...\\n")
results <- dbGetQuery(con, sql)
print(results)
'
  writeLines(quick_start_code, file.path(RELEASE_DIR, "quick_start.R"))
  
  cat("================================================================\n")
  cat("RELEASE PREPARED SUCCESSFULLY!\n")
  cat("Location: ", RELEASE_DIR, "\n")
  cat("Action:   Zip this folder and upload to Zenodo.\n")
  cat("================================================================\n")
}

main()