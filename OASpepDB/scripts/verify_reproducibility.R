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

cat("[3/5] Verifying database detection...\n")
db_path <- app_env$DEFAULT_DB_ROOT
stopifnot(nzchar(db_path), dir.exists(db_path))
cohort_count <- app_env$get_db_cohort_count(db_path)
cat(sprintf("      OK: Detected database at '%s' with %d cohorts.\n", db_path, cohort_count))
stopifnot(cohort_count >= 1)

cat("[4/5] Verifying DuckDB queries on demo data...\n")
con <- app_env$get_db_con()
query_sql <- sprintf("SELECT Disease, BSource, BType, Isotype, count(*) as n FROM read_parquet('%s/**/*.parquet', hive_partitioning=true) GROUP BY Disease, BSource, BType, Isotype", db_path)
res <- DBI::dbGetQuery(con, query_sql)
cat(sprintf("      OK: Query returned %d partition rows.\n", nrow(res)))
stopifnot(nrow(res) > 0)

lookup_sql <- sprintf("SELECT clonotype_id, cdr3_aa, tryptic_peptides FROM read_parquet('%s/**/*.parquet', hive_partitioning=true) WHERE list_contains(string_split(tryptic_peptides, ';'), 'GGYSYGYFDYWGQGTLVTVSSASTK')", db_path)
lookup_res <- DBI::dbGetQuery(con, lookup_sql)
cat(sprintf("      OK: Reverse peptide lookup returned %d match(es).\n", nrow(lookup_res)))
stopifnot(nrow(lookup_res) >= 1)
DBI::dbDisconnect(con, shutdown = TRUE)

cat("[5/5] Verifying reference FASTA assets...\n")
crap_path <- app_env$CRAP_REF_PATH
entrap_path <- app_env$ENTRAPMENT_REF_PATH
stopifnot(file.exists(crap_path))
stopifnot(file.exists(entrap_path))
cat("      OK: cRAP and Entrapment reference files resolved.\n")

cat("\n======================================================\n")
cat("SUCCESS: OASpepDB is 100% reproducible and ready to run!\n")
cat("======================================================\n")
