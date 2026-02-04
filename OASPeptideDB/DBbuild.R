################################################################################
## DBbuild.R
##
## Goal: Run the full DAT-DB build pipeline sequentially.
##
## Prerequisite:
##   - Ensure 'scripts/00_config.R' is set up correctly.
##   - Ensure renv is active (renv::restore() if needed).
##
## How to use:
##   1. Comment/Uncomment lines in `pipeline_scripts` to select steps.
##   2. Run: source("DBbuild.R")
################################################################################

# ==============================================================================
# 1) Pipeline Definition
#    (Comment out lines with '#' to skip specific steps)
# ==============================================================================
pipeline_scripts <- list(
  
  # --- Setup & Acquisition ---
  "Step 1: Download & Preprocess OAS"      = "scripts/01_download.R",
  
  # --- Processing ---
  "Step 2: Tryptic Digestion"              = "scripts/02_digestion.R",
  
  # --- Database Construction ---
  "Step 3: Build Fact Table (db_stage1)"   = "scripts/03_build_fact_table.R",
  "Step 4: Build Dimension/Metrics (db_stage2)" = "scripts/04_build_dimension_table.R",
  "Step 5: Combine & Finalize (db_stage3)" = "scripts/05_combine_db.R"
)

# ==============================================================================
# 2) Runner Function (Fail-Fast)
# ==============================================================================
run_step <- function(step_name, script_path) {
  cat("\n============================================================\n")
  cat(">>> Running:", step_name, "\n")
  cat("Script     :", script_path, "\n")
  cat("Start Time :", format(Sys.time(), "%H:%M:%S"), "\n")
  cat("============================================================\n\n")
  
  if (!file.exists(script_path)) {
    stop(sprintf("❌ FAILED: Script file not found: %s", script_path))
  }
  
  # Run the script in the global environment
  ok <- tryCatch(
    {
      source(script_path, local = FALSE)
      TRUE
    },
    error = function(e) {
      cat("\n------------------------------------------------------------\n")
      cat("❌ FAILED:", step_name, "\n")
      cat("Error message:\n")
      cat(conditionMessage(e), "\n")
      cat("------------------------------------------------------------\n\n")
      FALSE
    }
  )
  
  if (!isTRUE(ok)) {
    stop(sprintf("Pipeline stopped due to error in: %s", step_name))
  }
  
  cat(sprintf("\n✅ Completed: %s\n", step_name))
  invisible(TRUE)
}

# ==============================================================================
# 3) Execution Loop
# ==============================================================================
if (length(pipeline_scripts) == 0) {
  stop("pipeline_scripts list is empty. Please uncomment at least one step.")
}

cat("\n============================================================\n")
cat("🚀 STARTING PIPELINE \n")
cat("The following steps will be executed:\n")
for (nm in names(pipeline_scripts)) {
  cat(" - ", nm, "  [", pipeline_scripts[[nm]], "]\n", sep = "")
}
cat("============================================================\n\n")

# Iterate through selected scripts
for (nm in names(pipeline_scripts)) {
  run_step(nm, pipeline_scripts[[nm]])
}

cat("\n============================================================\n")
cat("🎉 PIPELINE COMPLETED SUCCESSFULLY 🎉\n")
cat("============================================================\n\n")

invisible(TRUE)