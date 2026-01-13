################################################################################
# scripts/07_part7_QC_final_vs_final_enriched_minimal.R
#
# Minimal Part7 QC:
#   - glimpse() both datasets
#   - compare: partitions (folder count), rows, columns
#
# SUCCESS criteria:
#   1) Same number of partition folders (Disease/BSource/BType/Isotype)
#   2) Same row count (Arrow metadata)
#   3) final_enriched contains ALL columns of final
#   4) final_enriched has additional enrichment columns
################################################################################

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
})

# ------------------------------------------------------------------------------
# CONFIG
# ------------------------------------------------------------------------------
parquet_dir_final   <- "D:/OAS/OAS_human_heavychain_disease_tryptic/parquet_db_partitioned/final"
parquet_dir_enrich  <- "D:/OAS/OASpeptideDB"

# Expected enrichment columns (must exist in final_enriched but not required in final)
enrich_added_cols <- c(
  "Disease_presence_peptide",
  "n_distinct_antibodies",
  "peptide_in_cdr3"
)

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
cat_line <- function(...) cat(..., "\n", sep = "")

count_partition_folders <- function(root) {
  dirs <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  part_dirs <- dirs[
    grepl("Disease=", dirs) &
      grepl("BSource=", dirs) &
      grepl("BType=", dirs) &
      grepl("Isotype=", dirs)
  ]
  length(part_dirs)
}

ds_info <- function(ds) {
  list(
    n_files = length(ds$files),
    n_rows  = ds$num_rows,   # Arrow metadata (fast)
    n_cols  = ncol(ds),
    cols    = names(ds)
  )
}

# ------------------------------------------------------------------------------
# 1) Open datasets (LAZY)
# ------------------------------------------------------------------------------
cat_line("Opening final dataset:\n  ", parquet_dir_final)
ds_final <- open_dataset(parquet_dir_final, format = "parquet")

cat_line("Opening final_enriched dataset:\n  ", parquet_dir_enrich)
ds_enrich <- open_dataset(parquet_dir_enrich, format = "parquet")

# ------------------------------------------------------------------------------
# 2) glimpse() overview
# ------------------------------------------------------------------------------
cat_line("\n--- final (Part2) ---")
glimpse(ds_final)

cat_line("\n--- final_enriched (Part6) ---")
glimpse(ds_enrich)

# ------------------------------------------------------------------------------
# 3) Compare rows + columns + partitions
# ------------------------------------------------------------------------------
info_final  <- ds_info(ds_final)
info_enrich <- ds_info(ds_enrich)

part_final  <- count_partition_folders(parquet_dir_final)
part_enrich <- count_partition_folders(parquet_dir_enrich)

cat_line("\n==================== QC SUMMARY ====================")
cat_line("Files:")
cat_line("  final        : ", info_final$n_files)
cat_line("  final_enriched: ", info_enrich$n_files)

cat_line("\nRows (Arrow metadata):")
cat_line("  final        : ", format(info_final$n_rows, big.mark = ","))
cat_line("  final_enriched: ", format(info_enrich$n_rows, big.mark = ","))

cat_line("\nColumns:")
cat_line("  final        : ", info_final$n_cols)
cat_line("  final_enriched: ", info_enrich$n_cols)

cat_line("\nPartition folders (Disease/BSource/BType/Isotype):")
cat_line("  final        : ", part_final)
cat_line("  final_enriched: ", part_enrich)
cat_line("====================================================\n")

# ------------------------------------------------------------------------------
# 4) SUCCESS / FAIL checks
# ------------------------------------------------------------------------------
ok <- TRUE

# (1) partitions
if (part_final != part_enrich) {
  ok <- FALSE
  cat_line("❌ FAIL: Partition folder count differs (final=", part_final,
           ", final_enriched=", part_enrich, ")")
} else {
  cat_line("✅ OK: Partition folder count matches (", part_final, ")")
}

# (2) rows (if Arrow reports them)
if (!is.na(info_final$n_rows) && !is.na(info_enrich$n_rows)) {
  if (info_final$n_rows != info_enrich$n_rows) {
    ok <- FALSE
    cat_line("❌ FAIL: Row count mismatch (final=", info_final$n_rows,
             ", final_enriched=", info_enrich$n_rows, ")")
  } else {
    cat_line("✅ OK: Row count matches")
  }
} else {
  cat_line("⚠️ NOTE: Arrow num_rows is NA for at least one dataset; row parity not checked.")
}

# (3) columns: final_enriched must contain all final columns
missing_from_enrich <- setdiff(info_final$cols, info_enrich$cols)
if (length(missing_from_enrich) > 0) {
  ok <- FALSE
  cat_line("❌ FAIL: final_enriched is missing columns from final:")
  print(missing_from_enrich)
} else {
  cat_line("✅ OK: final_enriched contains all final columns")
}

# (4) added enrichment columns exist
missing_added <- setdiff(enrich_added_cols, info_enrich$cols)
if (length(missing_added) > 0) {
  ok <- FALSE
  cat_line("❌ FAIL: final_enriched is missing enrichment columns:")
  print(missing_added)
} else {
  cat_line("✅ OK: enrichment columns present in final_enriched")
}

# (Optional) show what extra columns final_enriched has (besides final)
extra_in_enrich <- setdiff(info_enrich$cols, info_final$cols)
cat_line("\nExtra columns in final_enriched (vs final):")
print(extra_in_enrich)

# ------------------------------------------------------------------------------
# DONE
# ------------------------------------------------------------------------------
if (ok) {
  cat_line("\n🎉 SUCCESS: Part6 output looks correct (partitions match, rows match, columns preserved + added).")
} else {
  cat_line("\n🛑 QC FAILED: One or more checks did not pass. See messages above.")
}
