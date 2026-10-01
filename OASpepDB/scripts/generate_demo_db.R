# ==============================================================================
# Generate Lightweight Demo CDR3 Database for Reviewer Reproducibility
# Creates a minimal Hive-partitioned Parquet dataset for 5 sample cohorts.
# ==============================================================================

suppressPackageStartupMessages({
  library(arrow)
})

script_dir <- tryCatch({
  normalizePath(dirname(sys.frame(1)$ofile))
}, error = function(e) {
  getwd()
})

base_dir <- if (basename(script_dir) %in% c("scripts", "Scripts")) dirname(script_dir) else file.path(getwd(), "OASpepDB")
out_root <- file.path(base_dir, "demo_CDR3_db")

demo_records <- list(
  # 1. COVID-19
  list(
    disease = "COVID-19",
    disease_long = "Severe Acute Respiratory Syndrome Coronavirus 2 (COVID-19)",
    bsource = "PBMC",
    btype = "Memory-B-Cells",
    isotype = "IGHG1",
    clonotype_id = paste0("CLON_COV_", sprintf("%04d", 1:8)),
    cdr3_aa = c(
      "ARGGYSYGYFDY", "ARDLYYYYYGMDV", "ARHGGNYWFDY", "ARDLVRGVIYYYYGMDV",
      "ARGPQYYDFWSGYFDY", "ARVIGSGSRWFDP", "ARDTGYYDSSGYYY", "AKDFGSSWYFDY"
    ),
    v_call = c("IGHV3-23*01", "IGHV1-69*01", "IGHV4-39*01", "IGHV3-30*01", "IGHV1-18*01", "IGHV3-7*01", "IGHV5-51*01", "IGHV3-21*01"),
    d_call = c("IGHD3-10*01", "IGHD2-2*01", "IGHD6-19*01", "IGHD3-3*01", "IGHD1-26*01", "IGHD2-15*01", "IGHD5-12*01", "IGHD3-9*01"),
    j_call = c("IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ5*02", "IGHJ4*02", "IGHJ4*02"),
    v_id = c(97.2, 95.8, 98.1, 94.5, 96.0, 93.8, 97.5, 96.4),
    patients = as.integer(c(5, 4, 3, 3, 2, 2, 1, 1)),
    redundancy = as.integer(c(120, 85, 45, 30, 22, 18, 9, 5)),
    file = "Kim_et_al_2020_COV_PBMC.csv.gz"
  ),
  # 2. HIV
  list(
    disease = "HIV",
    disease_long = "Human Immunodeficiency Virus Infection",
    bsource = "PBMC",
    btype = "Memory-B-Cells",
    isotype = "IGHG1",
    clonotype_id = paste0("CLON_HIV_", sprintf("%04d", 1:6)),
    cdr3_aa = c(
      "ARGKYYDFWSGYYTYFDY", "ARDRGQLLYYYYYGMDV", "AREGTTTVTTFDY",
      "ARHGYYYGSGSYFDY", "ARVLPFDFWSGYPDY", "ARDLGIAAPYYYYFDY"
    ),
    v_call = c("IGHV1-2*02", "IGHV1-69*01", "IGHV3-23*01", "IGHV4-34*01", "IGHV1-46*01", "IGHV3-30*01"),
    d_call = c("IGHD3-10*01", "IGHD2-2*01", "IGHD6-13*01", "IGHD3-3*01", "IGHD1-1*01", "IGHD2-21*01"),
    j_call = c("IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02"),
    v_id = c(88.4, 91.2, 93.0, 89.5, 92.1, 90.7),
    patients = as.integer(c(4, 3, 2, 2, 1, 1)),
    redundancy = as.integer(c(95, 60, 35, 20, 12, 8)),
    file = "Schoofs_et_al_2016_HIV_PBMC.csv.gz"
  ),
  # 3. SLE
  list(
    disease = "SLE",
    disease_long = "Systemic Lupus Erythematosus",
    bsource = "PBMC",
    btype = "Plasma-Cells",
    isotype = "IGHG1",
    clonotype_id = paste0("CLON_SLE_", sprintf("%04d", 1:5)),
    cdr3_aa = c(
      "ARRRYSYGYFDY", "ARDLRYYYYGMDV", "ARHRGNYWFDY", "ARDLRRGVIYYYYGMDV", "ARRGPQYYDFWSGYFDY"
    ),
    v_call = c("IGHV3-23*01", "IGHV4-34*01", "IGHV1-69*01", "IGHV3-30*01", "IGHV1-18*01"),
    d_call = c("IGHD3-10*01", "IGHD2-2*01", "IGHD6-19*01", "IGHD3-3*01", "IGHD1-26*01"),
    j_call = c("IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ6*02", "IGHJ4*02"),
    v_id = c(96.1, 95.0, 97.4, 93.9, 95.8),
    patients = as.integer(c(3, 2, 2, 1, 1)),
    redundancy = as.integer(c(55, 40, 28, 15, 10)),
    file = "Tipton_et_al_2015_SLE_PBMC.csv.gz"
  ),
  # 4. Asthma
  list(
    disease = "Asthma",
    disease_long = "Bronchial Asthma",
    bsource = "PBMC",
    btype = "Memory-B-Cells",
    isotype = "IGHA1",
    clonotype_id = paste0("CLON_AST_", sprintf("%04d", 1:4)),
    cdr3_aa = c(
      "ARGWTYFDY", "ARDPVTGYYYYMDV", "ARSSSWYFDY", "AKDVGSYFDY"
    ),
    v_call = c("IGHV3-23*01", "IGHV1-69*01", "IGHV3-7*01", "IGHV3-21*01"),
    d_call = c("IGHD3-10*01", "IGHD2-15*01", "IGHD5-12*01", "IGHD3-9*01"),
    j_call = c("IGHJ4*02", "IGHJ6*02", "IGHJ4*02", "IGHJ4*02"),
    v_id = c(97.0, 96.2, 95.5, 98.0),
    patients = as.integer(c(2, 2, 1, 1)),
    redundancy = as.integer(c(32, 24, 14, 6)),
    file = "Berkers_et_al_2019_Asthma_PBMC.csv.gz"
  ),
  # 5. CLL
  list(
    disease = "CLL",
    disease_long = "Chronic Lymphocytic Leukemia",
    bsource = "PBMC",
    btype = "Unsorted-B-Cells",
    isotype = "IGHM",
    clonotype_id = paste0("CLON_CLL_", sprintf("%04d", 1:4)),
    cdr3_aa = c(
      "ARDLRYYYYGMDV", "AREGTTTVTTFDY", "ARHGYYYGSGSYFDY", "ARDLGIAAPYYYYFDY"
    ),
    v_call = c("IGHV1-69*01", "IGHV3-23*01", "IGHV4-34*01", "IGHV3-30*01"),
    d_call = c("IGHD2-2*01", "IGHD6-13*01", "IGHD3-3*01", "IGHD2-21*01"),
    j_call = c("IGHJ6*02", "IGHJ4*02", "IGHJ4*02", "IGHJ4*02"),
    v_id = c(99.2, 98.5, 99.0, 98.8),
    patients = as.integer(c(3, 2, 1, 1)),
    redundancy = as.integer(c(400, 150, 80, 45)),
    file = "Stamatopoulos_et_al_2017_CLL_PBMC.csv.gz"
  )
)

fr3_tail <- "YYCAR"
fwr4_aa <- "WGQGTLVTVSS"
isotype_tail <- "ASTK"

for (rec in demo_records) {
  n <- length(rec$clonotype_id)
  cassettes <- paste0(fr3_tail, rec$cdr3_aa, fwr4_aa, isotype_tail)
  cdr3_lens <- as.integer(nchar(rec$cdr3_aa))
  
  # Generate tryptic peptides: [KR](?!P) digestion simulation
  tryptics <- sapply(cassettes, function(cas) {
    # Splits at K or R
    peps <- unlist(strsplit(gsub("([KR])", "\\1\n", cas), "\n"))
    peps <- peps[nchar(peps) >= 6]
    paste(peps, collapse = ";")
  }, USE.NAMES = FALSE)

  df <- data.frame(
    clonotype_id = rec$clonotype_id,
    cdr3_cassette = cassettes,
    cdr3_aa = rec$cdr3_aa,
    cdr3_len = cdr3_lens,
    fr3_tail = rep(fr3_tail, n),
    fwr4_aa = rep(fwr4_aa, n),
    isotype_tail = rep(isotype_tail, n),
    tryptic_peptides = tryptics,
    v_call = rec$v_call,
    d_call = rec$d_call,
    j_call = rec$j_call,
    v_identity = rec$v_id,
    sequence_alignment_aa = paste0("EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAV", cassettes),
    Redundancy = rec$redundancy,
    N_Patients = rec$patients,
    Data_file = rep(rec$file, n),
    Disease_long = rep(rec$disease_long, n),
    stringsAsFactors = FALSE
  )

  part_dir <- file.path(
    out_root,
    paste0("Disease=", rec$disease),
    paste0("BSource=", rec$bsource),
    paste0("BType=", rec$btype),
    paste0("Isotype=", rec$isotype)
  )
  dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(part_dir, "data_0.parquet")
  arrow::write_parquet(df, out_file, compression = "zstd", compression_level = 9)
  message(sprintf("Wrote %d demo records to %s", nrow(df), out_file))
}

message("Demo CDR3 database generation complete at: ", out_root)
