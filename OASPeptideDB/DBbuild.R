################################################################################
## DBbuild.R
##
## Run selected pipeline scripts sequentially (Part1 → Part7).
##
## How to use
## ----------
## 1) Edit the `pipeline_scripts` list below.
## 2) Comment out any lines you DON'T want to run.
## 3) Run:
##      source("DBbuild.R")
##
################################################################################

# ==============================================================================
# 1) Choose which scripts to run (comment out any you want to skip)
# ==============================================================================
pipeline_scripts <- list(
  "Part1 – Download & preprocess OAS"          = "scripts/01_download.R"#,
  #"Part2 – Digestion & Parquet DB build"       = "scripts/02_part2_digestion.R",
  #"Part3 – Disease presence"                   = "scripts/03_part3_disease_presence.R",
  #"Part4 – Peptide uniqueness"                 = "scripts/04_part4_peptide_uniqueness.R",
  #"Part5 – Peptide in CDR3"                    = "scripts/05_part5_cdr3.R",
  #"Part6 – Zenodo PeptideDB bundle"            = "scripts/06_part6_zenodo_bundle.R",
  #"Part7 – QC check of part 6"                 = "scripts/07_part7_QC_final_enriched.R"
)

# ==============================================================================
# 2) Runner (fail-fast)
# ==============================================================================
run_step <- function(step_name, script_path) {
  cat("\n============================================================\n")
  cat(">>> Running:", step_name, "\n")
  cat("Script     :", script_path, "\n")
  cat("============================================================\n\n")
  
  if (!file.exists(script_path)) {
    stop(sprintf("❌ FAILED: script not found: %s", script_path))
  }
  
  ok <- tryCatch(
    {
      source(script_path, local = TRUE)
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
    stop(sprintf("Pipeline stopped at: %s", step_name))
  }
  
  cat(sprintf("\n✅ Completed: %s\n\n", step_name))
  invisible(TRUE)
}

# ==============================================================================
# 3) Execute
# ==============================================================================
if (length(pipeline_scripts) == 0) {
  stop("pipeline_scripts is empty. Add at least one script to run.")
}

cat("\n============================================================\n")
cat("Pipeline start\n")
cat("Scripts to run:\n")
for (nm in names(pipeline_scripts)) {
  cat(" - ", nm, "  (", pipeline_scripts[[nm]], ")\n", sep = "")
}
cat("============================================================\n\n")

for (nm in names(pipeline_scripts)) {
  run_step(nm, pipeline_scripts[[nm]])
}

cat("\n============================================================\n")
cat("🎉 PIPELINE COMPLETED SUCCESSFULLY 🎉\n")
cat("============================================================\n\n")

invisible(TRUE)
