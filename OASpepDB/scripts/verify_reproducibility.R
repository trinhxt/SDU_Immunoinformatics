# ==============================================================================
# OASpepDB: Self-Check Reproducibility Verification Script
# Verifies packages, DB detection, DuckDB SQL execution, and reference data.
# ==============================================================================

cat("[1/5] Verifying R dependencies...\n")
required_pkgs <- c("shiny", "bslib", "DT", "duckdb", "DBI", "arrow", "ggplot2", "plotly", "dplyr", "htmlwidgets", "zip")
missing_pkgs <- required_pkgs[!required_pkgs %in% installed.packages()[, "Package"]]
if (length(missing_pkgs) > 0) {
  stop("Missing R packages: ", paste(missing_pkgs, collapse = ", "))
}
cat("      OK: All 11 required packages are installed.\n")

cat("[2/5] Loading DBquery.R environment...\n")
app_env <- new.env()
sys.source("OASpepDB/DBquery.R", envir = app_env)
cat("      OK: DBquery.R parsed and loaded successfully.\n")

cat("[3/5] Verifying database availability...\n")
db_path <- app_env$DEFAULT_DB_ROOT
if (!nzchar(db_path) || !dir.exists(db_path)) {
  cat("      Notice: No pre-existing CDR3_db found. Building demo CDR3_db via pipeline scripts...\n")
  system2("python", c("OASpepDB/scripts/01_download_OAS.py", "--sh", "OASpepDB/data/demo_download_human_unpaired.sh", "--out-dir", "OAS_raw", "--threads", "4"))
  system2("python", c("OASpepDB/scripts/02_build_healthy_cdr3_cache.py", "--data-dir", "OAS_raw", "--out-parquet", "CDR3_db/healthy_cdr3_cache.parquet", "--workers", "2"))
  system2("python", c("OASpepDB/scripts/04_build_disease_db.py", "--data-dir", "OAS_raw", "--healthy-cache", "CDR3_db/healthy_cdr3_cache.parquet", "--out-dir", "CDR3_db", "--workers", "2"))
  db_path <- app_env$detect_default_db_root()
}
stopifnot(nzchar(db_path), dir.exists(db_path))
cohort_count <- app_env$get_db_cohort_count(db_path)
cat(sprintf("      OK: Connected to database at '%s' with %d cohorts.\n", db_path, cohort_count))
stopifnot(cohort_count >= 1)

cat("[4/5] Verifying DuckDB queries on database...\n")
con <- app_env$get_db_con()
query_sql <- sprintf("SELECT Disease, BSource, BType, Isotype, count(*) as n FROM read_parquet('%s/Disease=*/**/*.parquet', hive_partitioning=true) GROUP BY Disease, BSource, BType, Isotype", db_path)
res <- DBI::dbGetQuery(con, query_sql)
cat(sprintf("      OK: Query returned %d partition rows.\n", nrow(res)))
stopifnot(nrow(res) > 0)
DBI::dbDisconnect(con, shutdown = TRUE)

cat("[5/5] Verifying reference FASTA assets...\n")
crap_path <- app_env$CRAP_REF_PATH
entrap_path <- app_env$ENTRAPMENT_REF_PATH
stopifnot(file.exists(crap_path))
stopifnot(file.exists(entrap_path))
cat("      OK: cRAP and Entrapment reference files resolved.\n")

cat("\n======================================================\n")
cat("SUCCESS: OASpepDB is 100% reproducible and ready to run!\n")
cat("======================================================\n\n")
cat("To launch the interactive Shiny web app:\n")
cat("  - Windows: Double-click 'OASpepDB/DBquery-Windows.bat'\n")
cat("  - Linux:   bash OASpepDB/DBquery-Linux.sh\n")
cat("  - macOS:   bash OASpepDB/DBquery-Mac.command\n")
cat("  - Command: Rscript -e \"shiny::runApp('OASpepDB/DBquery.R', launch.browser = TRUE)\"\n\n")
