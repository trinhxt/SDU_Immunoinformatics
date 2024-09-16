#-------------------------------------------------------------------------------
# Ensure required packages are installed and loaded with auto-update option
# Set a custom library path
setwd("/work/Immunoinformatics/Antibody-DB/")
custom_lib <- "/work/Immunoinformatics/R-packages-RStudio/"
#custom_lib <- "/work/Immunoinformatics/R-packages-Ubuntu/"
.libPaths(custom_lib)

# Define the required packages
cran_packages <- c("arrow", "dplyr", "jsonlite", "seqinr", "data.table", "fs", "stringr", "stringi", "duckdb", "DBI", "RSQLite")
bioc_packages <- c("cleaver", "Biostrings")

# Load all packages
lapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

#-------------------------------------------------------------------------------

# Read metadata file
metadata <- fread("OAS_metadata.csv")
colnames(metadata)

# Summarize the data: group by Disease and sum the Unique sequences
disease_summary <- metadata[, .(Total_unique_sequences = sum(`Unique sequences`)), by = Disease]

# Order by Total_unique_sequences in descending order
disease_summary <- disease_summary[order(-Total_unique_sequences)]

# Add the Disease_no column for diseases other than "None"
disease_summary[Disease != "None", Disease_no := paste0("D", seq_len(.N))]

# Keep Disease_no as "None" for the "None" Disease
disease_summary[Disease == "None", Disease_no := "None"]

#-------------------------------------------------------------------------------
# Define the working folder and the target folder with the new naming convention
antibody_folder <- "OAS_full"
peptide_folder1 <- "OAS_tryptic/Disease_index_ab/"
peptide_folder2 <- "OAS_tryptic/Disease_index_filename/"
duckdb_folder_pep <- "OAS_pepDB/SQL"
csv_folder_pep <- "OAS_pepDB/CSV"
duckdb_folder_ab <- "OAS_abDB/SQL"
progress_file1 <- file.path("Log_files", "Part4_pepDB_prep.json")

# Create the target folders if they don't exist
if (!dir_exists(duckdb_folder_pep)) {
  dir_create(duckdb_folder_pep)
}

if (!dir_exists(duckdb_folder_ab)) {
  dir_create(duckdb_folder_ab)
}

if (!dir_exists(csv_folder_pep)) {
  dir_create(csv_folder_pep)
}

# Load progress from previous run if exists
if (file_exists(progress_file1)) {
  completed_diseases_redun <- fromJSON(progress_file1)
} else {
  completed_diseases_redun <- c()
}

# Filter the disease_summary to only include diseases that haven't been processed yet
if (length(completed_diseases_redun) > 0) {
  disease_summary <- disease_summary[!(Disease_no %in% completed_diseases_redun)]
}

# Get the list of all files in the working folder (antibody data)
antibody_files <- list.files(path = antibody_folder, pattern = "\\.csv\\.gz$", full.names = TRUE)
file_base_names <- path_file(antibody_files)

# Get list of all peptide files in OAS_tryptic/Disease_index_ab
peptide_files1 <- list.files(path = peptide_folder1, pattern = "\\.csv\\.gz$", full.names = TRUE)
# Get list of all peptide files in OAS_tryptic/Disease_index_filename
peptide_files2 <- list.files(path = peptide_folder2, pattern = "\\.csv\\.gz$", full.names = TRUE)

# Initialize a list to store error logs
error_log <- list()

# Initialize the output CSV file path for metadata2
metadata2_file <- file.path("OAS_metadata2.csv")

#-------------------------------------------------------------------------------
# Loop through each disease in disease_summary
for (disease_no in rev(disease_summary$Disease_no)) {
  
  # Skip if the disease_no is "None"
  if (disease_no == "None") next
  
  # Get the corresponding disease name
  disease_name <- disease_summary[Disease_no == disease_no, Disease]
  
  # Display a message indicating which disease is being processed
  message(sprintf("Processing disease: %s (Disease No: %s)", disease_name, disease_no))
  
  # Get the file list for the current disease
  files_for_disease <- metadata[Disease == disease_name, Filename]
  
  # Filter the files to process by matching the base names to files_for_disease
  files_to_process <- antibody_files[file_base_names %in% files_for_disease]
  
  # Get the correct tryptic files for the current disease_no
  tryptic_file1 <- peptide_files1[grepl(paste0("^", disease_no, "_"), basename(peptide_files1))]
  tryptic_file2 <- peptide_files2[grepl(paste0("^", disease_no, "_"), basename(peptide_files2))]
  
  #-----------------------------------------------------------------------------
  # Step 1: Read tryptic_file1 and prepare table1 with unique sequences and N_antibody, and Sequence_length
  
  # Read the tryptic data from the file using read_csv_arrow from the arrow package
  tryptic_data1 <- read_csv_arrow(tryptic_file1)
  
  # Create table1 by extracting unique sequences, counting antibodies, and calculating sequence lengths
  table1 <- tryptic_data1 %>%
    group_by(Sequence) %>%
    reframe(
      N_antibody = n(), # Total number of antibodies related to each sequence
      Length_aa = nchar(Sequence) # Directly calculate length of each unique sequence
    ) %>%
    distinct(Sequence, .keep_all = TRUE) %>%  # Ensure only unique sequences are retained
    as.data.table()
  
  # clean memory
  rm(tryptic_data1)
  #-----------------------------------------------------------------------------
  # Step 2: Read all antibody files for the related disease to extract CDR3 column
  CDR3 <- c()
  CDR3 <- unique(unlist(lapply(files_to_process, function(file) {
    tryCatch({
      # Read the file using read_csv_arrow
      data <- read_csv_arrow(file)
      # Extract unique CDR3 sequences from the file
      unique_cdr3 <- unique(data$cdr3_aa)
      return(unique_cdr3)
    }, error = function(e) {
      # Log the error and skip the file
      error_log <<- append(error_log, list(list(file = file, error = conditionMessage(e))))
      return(NULL)
    })
  })))
  
  # Digest CDR3 peptides to tryptic peptides
  CDR3_tryptic_all <- character()
  # Process in chunks of 5,000,000
  total_chunks <- ceiling(length(CDR3) / 5e6) # 5,000,000 per chunk
  for (chunk in seq_len(total_chunks)) {
    cat(sprintf("     Processing chunk %d of %d\n", chunk, total_chunks))
    
    # Define chunk start and end
    start_index <- (chunk - 1) * 5e6 + 1
    end_index <- min(chunk * 5e6, length(CDR3))
    
    # Extract chunk of sequences
    seq_ab_chunk <- CDR3[start_index:end_index]
    names_ab_chunk <- paste0("CDR3_", seq(start_index, end_index))
    
    tryptic_combined <- AAStringSet(setNames(seq_ab_chunk, names_ab_chunk))
    tryptic_combined <- cleave(tryptic_combined, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
    CDR3_tryptic <- unique(as.character(unlist(tryptic_combined)))
    
    # Accumulate results in the combined data table
    CDR3_tryptic_all <- unique(c(CDR3_tryptic_all, CDR3_tryptic))
    
    # Clean memory
    rm(list = c("CDR3_tryptic", "names_ab_chunk", "tryptic_combined"))
    gc()
  }
  
  table1$CDR3 <- ifelse(table1$Sequence %in% CDR3_tryptic_all, 1, 0)
  
  max_antibody <- max(table1$N_antibody)
  N_peptide_unique <- uniqueN(table1$Sequence[table1$N_antibody==1])
  N_peptide_nonunique <- uniqueN(table1$Sequence[table1$N_antibody>1])
  N_peptide_IGHA <- uniqueN(table1$Sequence[table1$Isotype == "IGHA"])
  N_peptide_IGHD <- uniqueN(table1$Sequence[table1$Isotype == "IGHD"])
  N_peptide_IGHE <- uniqueN(table1$Sequence[table1$Isotype == "IGHE"])
  N_peptide_IGHG <- uniqueN(table1$Sequence[table1$Isotype == "IGHG"])
  N_peptide_IGHM <- uniqueN(table1$Sequence[table1$Isotype == "IGHM"])
  N_peptide_Bulk <- uniqueN(table1$Sequence[table1$Isotype == "Bulk"])
  
  #-----------------------------------------------------------------------------
  # Step 3: Read tryptic_file2 and extract Peptides, Filename, Isotype, Btype, Bsource, Patient
  tryptic_data2 <- read_csv_arrow(tryptic_file2)
  
  # Subset tryptic_data2 based on table1$Sequence
  table2 <- tryptic_data2 %>%
    filter(Sequence %in% table1$Sequence)
  
  # Clear tryptic_data2 from memory
  rm(tryptic_data2)
  
  # Join table2 with specific columns from metadata based on Filename
  table2 <- table2 %>%
    left_join(metadata %>% select(Filename, Patient, BSource, BType, Isotype), by = "Filename")
  
  
  # Add a column "N_patient" to count the number of unique Patients per Sequence
  table2 <- table2 %>%
    group_by(Sequence) %>%
    mutate(N_patient = n_distinct(Patient)) %>%
    ungroup()
  
  # Join table2 with table1 based on Sequence
  table1 <- table2 %>%
    left_join(table1, by = "Sequence")

  N_peptide <- uniqueN(table1$Sequence)  # Number of unique peptides in table1
  max_patient <- max(table1$N_patient)
  
  
  #-----------------------------------------------------------------------------
  # Step 4: Save table1 to duckDB database
  
  # Get the tryptic file name without the .csv.gz extension
  tryptic_file1name <- basename(tryptic_file1)  # Extracts the filename from the path
  db_filename <- tools::file_path_sans_ext(tools::file_path_sans_ext(tryptic_file1name))  # Removes both .csv and .gz extensions
  
  # Specify the databas and CSV output folder paths
  db_file <- file.path(duckdb_folder_pep, paste0(db_filename, ".duckdb"))
  csv_file <- file.path(csv_folder_pep, paste0(db_filename, ".csv.gz"))
  
  # Create a SQLite connection
  #con <- dbConnect(RSQLite::SQLite(), dbname = db_file)
  # Create a connection to a DuckDB database file 
  con <- dbConnect(duckdb::duckdb(), dbdir = db_file)
  
  # Save table1 to the SQLite database
  dbWriteTable(con, "table1", table1, overwrite = TRUE)
  
  # Close the database connection
  dbDisconnect(con)
  
  # Save table1 to a compressed CSV file
  fwrite(table1, csv_file, compress = "gzip")
  
  # Display a message to confirm the successful saving of data
  message(sprintf("     Saved peptides (table1) to duckDB and CSV."))
  
  # Clear table1 from memory
  rm(table1, table2)
  
  
  #-----------------------------------------------------------------------------
  # Step 5: make table2 to store antibody data and save to duckDB database
  # Process files and combine them into one data table using rbindlist for speed
  table2 <- rbindlist(lapply(files_to_process, function(file) {
    tryCatch({
      # Read the file using read_csv_arrow
      data <- read_csv_arrow(file, col_select = c("sequence_alignment_aa",	"v_call",	"j_call",	"cdr3_aa"))
      data$Filename <- basename(file)
      return(data)
    }, error = function(e) {
      # Log the error and skip the file
      error_log <<- append(error_log, list(list(file = file, error = conditionMessage(e))))
      return(NULL)
    })
  }), use.names = TRUE, fill = TRUE)  # Ensure consistent column names and fill missing columns
  table2 <- unique(table2)
  
  N_antibody <- uniqueN(table2$sequence_alignment_aa) # Number of unique antibody sequences in table2
  
  # Save table2 to duckDB database
  # Get the tryptic file name without the .csv.gz extension
  tryptic_file1name <- basename(tryptic_file1)  # Extracts the filename from the path
  db_filename <- tools::file_path_sans_ext(tools::file_path_sans_ext(tryptic_file1name))  # Removes both .csv and .gz extensions
  
  # Specify the database output folder paths
  db_file <- file.path(duckdb_folder_ab, paste0(db_filename, ".duckdb"))
  
  # Create a database connection
  #con <- dbConnect(RSQLite::SQLite(), dbname = db_file)
  con <- dbConnect(duckdb::duckdb(), dbdir = db_file)
  
  # Save table2 to the SQLite database
  dbWriteTable(con, "table2", table2, overwrite = TRUE)
  
  # Close the database connection
  dbDisconnect(con)
  
  # Display a message to confirm the successful saving of data
  message(sprintf("     Saved antibodies (table2) to duckDB."))
  
  # Clear table2 from memory
  rm(table2)       
  
  
  
  #-----------------------------------------------------------------------------
  # Update metadata2
  N_patient <- uniqueN(metadata[metadata$Disease == disease_name, Patient])  # Number of unique patients
  
  # Store the results in a data.table for metadata2
  metadata2_entry <- data.table(
    Disease_no = disease_no,
    Disease = disease_name,
    N_antibody = N_antibody,
    N_peptide = N_peptide,
    N_patient = N_patient,
    max_antibody = max_antibody,
    max_patient = patient,
    N_peptide_unique = N_peptide_unique,
    N_peptide_nonunique = N_peptide_nonunique,
    N_peptide_IGHA = N_peptide_IGHA ,
    N_peptide_IGHD = N_peptide_IGHD ,
    N_peptide_IGHE = N_peptide_IGHE ,
    N_peptide_IGHG = N_peptide_IGHG ,
    N_peptide_IGHM = N_peptide_IGHM ,
    N_peptide_Bulk = N_peptide_Bulk 
  )
  
  # Save metadata2 incrementally to the CSV file
  fwrite(metadata2_entry, metadata2_file, append = TRUE, col.names = !file.exists(metadata2_file))
  
  
  
  #-----------------------------------------------------------------------------
  # Update the list of completed diseases and save progress
  completed_diseases_redun <- c(completed_diseases_redun, disease_no)
  write_json(completed_diseases_redun, progress_file1)
  
  # Trigger garbage collection to free up memory
  gc()             
  
}


# Save the error log if there were any errors
if (length(error_log) > 0) {
  write_json(error_log, file.path(duckdb_folder_pep, "error_log.json"))
}

# At this point, all diseases have been processed and the corresponding files have been saved in the target folder