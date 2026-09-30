#!/usr/bin/env Rscript
################################################################################
## OASpepDB: Step 03 - In Silico Digestion of Human Reference Proteomes
## (03_insilico_digest.R)
## Purpose:
##   Perform in silico tryptic digestion (Trypsin/P, missed cleavages 0-2,
##   length 6-45 aa, preserving native I and L residues) on two human reference
##   proteome datasets:
##     1. UniProt Swiss-Prot (Canonical + Reviewed Isoforms)
##     2. NCBI RefSeq Human Proteome (GRCh38.p14: NP_ + XP_)
##   Outputs:
##     - Parquet: Data/negative_human_reference_peptides.parquet
##     - DuckDB:  Data/negative_reference.duckdb (indexed B-Tree table)
##
## Design notes:
##   - Native amino acid sequences are strictly preserved (no I->L conversion;
##     isobaric normalization is deferred to the FASTA export stage).
##   - Chunked streaming (5,000 proteins/chunk) for minimal memory footprint.
##   - 100% English codebase and logging.
################################################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(data.table)
  library(arrow)
  library(duckdb)
})

# ==============================================================================
# 1. CONFIGURATION & FILE PATHS
# ==============================================================================
cat("\n")
cat("================================================================================\n")
cat("  OASpepDB: IN SILICO DIGESTION & NEGATIVE REFERENCE GENERATOR                  \n")
cat("================================================================================\n")

base_dir <- if (dir.exists(file.path(getwd(), "Data"))) getwd() else if (dir.exists(file.path(getwd(), "OASpepDB", "Data"))) file.path(getwd(), "OASpepDB") else "C:/Users/TXT/Documents/GitHub/SDU_Immunoinformatics/OASpepDB"
data_dir <- file.path(base_dir, "Data")

fasta_uniprot <- file.path(data_dir, "UniProt_SP_canonical_isoform_2026_09_26.fasta")
fasta_refseq  <- file.path(data_dir, "GCF_000001405.40_GRCh38.p14_protein.faa")

out_parquet <- file.path(data_dir, "negative_human_reference_peptides.parquet")
out_duckdb  <- file.path(data_dir, "negative_reference.duckdb")

# Digestion parameters
MIN_LEN <- 6L
MAX_LEN <- 45L
MAX_MISSED <- 2L
CHUNK_SIZE <- 5000L

cat("Data Directory:     ", data_dir, "\n")
cat("UniProt FASTA:      ", basename(fasta_uniprot), "\n")
cat("RefSeq FASTA:       ", basename(fasta_refseq), "\n")
cat("Output Parquet:     ", basename(out_parquet), "\n")
cat("Output DuckDB:      ", basename(out_duckdb), "\n")
cat(sprintf("Digestion Rules:    Trypsin/P | Missed 0-%d | Length %d-%d aa | Native I/L Preserved\n",
            MAX_MISSED, MIN_LEN, MAX_LEN))
cat("================================================================================\n\n")

# Verify input files exist
if (!file.exists(fasta_uniprot)) stop("UniProt FASTA file not found: ", fasta_uniprot)
if (!file.exists(fasta_refseq))  stop("RefSeq FASTA file not found: ", fasta_refseq)

# ==============================================================================
# 2. HIGH-PERFORMANCE IN SILICO TRYPTIC DIGESTION FUNCTION
# ==============================================================================
digest_protein_chunk <- function(seq_chunk, min_len = 6L, max_len = 45L, missed = 2L) {
  out <- vector("list", length(seq_chunk))
  
  for (k in seq_along(seq_chunk)) {
    seq <- seq_chunk[[k]]
    # Cleave after Lysine (K) or Arginine (R)
    cuts <- gregexpr("[KR]", seq)[[1]]
    cuts <- cuts[cuts > 0]
    n <- nchar(seq)
    
    if (length(cuts) == 0) {
      cuts <- n
    } else if (cuts[length(cuts)] != n) {
      cuts <- c(cuts, n)
    }
    cuts <- c(0L, cuts)
    m <- length(cuts)
    
    p_list <- character()
    for (mc in 0:missed) {
      if (m - 1L - mc <= 0L) next
      starts <- cuts[1:(m - 1L - mc)] + 1L
      ends   <- cuts[(2L + mc):m]
      lens   <- ends - starts + 1L
      
      valid <- (lens >= min_len) & (lens <= max_len)
      if (any(valid)) {
        p_list <- c(p_list, substring(seq, starts[valid], ends[valid]))
      }
    }
    
    # Store native peptides without I->L substitution
    out[[k]] <- p_list
  }
  
  # Deduplicate within chunk to minimize memory usage
  unique(unlist(out, use.names = FALSE))
}

process_fasta_file <- function(file_path, label) {
  cat(sprintf("[1/2] Loading %s...", label))
  t0 <- Sys.time()
  aa_set <- Biostrings::readAAStringSet(file_path)
  seqs <- as.character(aa_set)
  rm(aa_set); invisible(gc())
  n_total <- length(seqs)
  t_load <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf(" Done! (%d proteins, %.2f seconds)\n", n_total, t_load))
  
  cat(sprintf("[2/2] In silico digesting %d proteins (chunk size = %d)...\n", n_total, CHUNK_SIZE))
  t1 <- Sys.time()
  
  chunks <- split(seqs, ceiling(seq_along(seqs) / CHUNK_SIZE))
  chunk_results <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    chunk_results[[i]] <- digest_protein_chunk(chunks[[i]], MIN_LEN, MAX_LEN, MAX_MISSED)
    if (i %% 5 == 0 || i == length(chunks)) {
      cat(sprintf("   - Processed %d / %d chunks (%.1f%%)...\n",
                  i, length(chunks), (i / length(chunks)) * 100))
    }
  }
  
  rm(seqs, chunks); invisible(gc())
  unique_peps <- unique(unlist(chunk_results, use.names = FALSE))
  rm(chunk_results); invisible(gc())
  
  t_dig <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  cat(sprintf("   => Finished %s: %d unique native peptides (%.2f seconds)\n\n",
              label, length(unique_peps), t_dig))
  
  return(unique_peps)
}

# ==============================================================================
# 3. DIGEST UNIPROT AND REFSEQ DATASETS
# ==============================================================================
time_start_all <- Sys.time()

cat(">>> STEP 1: DIGESTING UNIPROT CANONICAL + ISOFORMS <<<\n")
uniprot_peps <- process_fasta_file(fasta_uniprot, "UniProt Swiss-Prot")

cat(">>> STEP 2: DIGESTING NCBI REFSEQ (GRCh38.p14) <<<\n")
refseq_peps <- process_fasta_file(fasta_refseq, "NCBI RefSeq Human")

# ==============================================================================
# 4. GLOBAL MERGE AND DEDUPLICATION
# ==============================================================================
cat(">>> STEP 3: GLOBAL MERGE AND DEDUPLICATION <<<\n")
t_merge <- Sys.time()

all_peptides <- unique(c(uniprot_peps, refseq_peps))
rm(uniprot_peps, refseq_peps); invisible(gc())

total_unique_peptides <- length(all_peptides)
cat(sprintf("Total unique native human reference peptides: %s\n",
            format(total_unique_peptides, big.mark = ",")))
cat(sprintf("Merge and deduplication time: %.2f seconds\n\n",
            as.numeric(difftime(Sys.time(), t_merge, units = "secs"))))

# Build compact data.table
ref_dt <- data.table(
  peptide = all_peptides,
  len     = as.integer(nchar(all_peptides))
)
rm(all_peptides); invisible(gc())

# ==============================================================================
# 5. EXPORT TO PARQUET & DUCKDB INDEX
# ==============================================================================
cat(">>> STEP 4: EXPORTING TO PARQUET & DUCKDB <<<\n")

# Write Parquet with ZSTD compression
t_io <- Sys.time()
cat("Writing Parquet file: ", out_parquet, "...\n")
arrow::write_parquet(ref_dt, out_parquet, compression = "zstd", compression_level = 7)
parquet_size_mb <- file.size(out_parquet) / (1024^2)
cat(sprintf("   => Parquet export completed: %.2f MB\n", parquet_size_mb))

# Ingest into indexed DuckDB database
cat("Creating indexed DuckDB database: ", out_duckdb, "...\n")
if (file.exists(out_duckdb)) unlink(out_duckdb)
con <- dbConnect(duckdb::duckdb(), out_duckdb)

dbExecute(con, "SET threads TO 4;")
dbExecute(con, sprintf("
  CREATE TABLE negative_peptides AS 
  SELECT peptide, len 
  FROM read_parquet('%s');
", gsub("\\\\", "/", out_parquet)))

cat("Building UNIQUE B-Tree index on 'peptide' column...\n")
dbExecute(con, "CREATE UNIQUE INDEX idx_peptide ON negative_peptides(peptide);")

row_count <- dbGetQuery(con, "SELECT COUNT(*) AS total FROM negative_peptides;")$total
dbDisconnect(con, shutdown = TRUE)

cat(sprintf("   => Successfully populated DuckDB: %s indexed records!\n\n",
            format(row_count, big.mark = ",")))

# ==============================================================================
# 6. QUALITY CONTROL & SANITY BENCHMARKS
# ==============================================================================
cat(">>> STEP 5: SANITY QUALITY CHECKS <<<\n")

# Standard housekeeping and germline antibody peptides (Native sequence with original I/L)
test_peptides <- c(
  "LVNEVTEFAK",         # Serum Albumin (HSA)
  "FKDLGEENFK",         # Serum Albumin (HSA)
  "SYELPDGQVITIGNER",   # Actin cytoplasmic 1 (ACTB, native with I)
  "GALQNIIPASTGAAK",    # GAPDH (native with I)
  "AEDTAVYYCAK",        # Germline IGHV3-23 FR3 tryptic peptide
  "ASTKGPSVFPLAPSSK"    # Germline IGHG1 CH1 tryptic peptide
)

qc_hits <- test_peptides %in% ref_dt$peptide
for (idx in seq_along(test_peptides)) {
  status <- if (qc_hits[idx]) "[PASS]" else "[FAIL]"
  cat(sprintf("   %s Benchmark peptide: %-22s -> %s\n", status, test_peptides[idx],
              ifelse(qc_hits[idx], "Found in Negative Index", "MISSING")))
}

if (!all(qc_hits)) {
  warning("Warning: Some benchmark peptides were not found in the negative index!")
} else {
  cat("\n==> ALL QUALITY CONTROL BENCHMARKS PASSED (100%)!\n")
}

time_total <- as.numeric(difftime(Sys.time(), time_start_all, units = "secs"))
cat("\n================================================================================\n")
cat(sprintf("IN SILICO DIGESTION COMPLETE (Total elapsed: %.2f seconds)!\n", time_total))
cat(sprintf("Parquet Artifact: %s (%.2f MB)\n", out_parquet, parquet_size_mb))
cat(sprintf("DuckDB Artifact:  %s (%.2f MB)\n", out_duckdb, file.size(out_duckdb) / (1024^2)))
cat("================================================================================\n")
