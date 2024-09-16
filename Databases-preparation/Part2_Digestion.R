#-------------------------------------------------------------------------------
# Set a custom library path (optional)
custom_lib <- "/work/Immunoinformatics/R-packages-RStudio/"
#custom_lib <- "/work/Immunoinformatics/R-packages-Ubuntu/"
.libPaths(custom_lib)

# Load all packages
cran_packages <- c("arrow", "dplyr", "jsonlite", "seqinr", "data.table", "fs", "stringr", "stringi", "parallel", "DBI", "RSQLite", "tools")
bioc_packages <- c("cleaver", "Biostrings", "GenomeInfoDb")
lapply(c(cran_packages, bioc_packages), require, character.only = TRUE)



################################################################################
## In silico digestion of proteins in UniProt and NCBI-RefSeq
## Read UniProt data
fasta_file <- ("UniProt_TR_SP_Human_2024_07_25.fasta")
# Read the FASTA file
sequences <- readAAStringSet(fasta_file) # use readAAStringSet for protein sequences
# Extract protein names and sequences
protein_names <- names(sequences)
protein_sequences <- as.character(sequences)
# Create a data table with protein names and sequences
UniProt_protein <- data.table(
  ProteinName = protein_names,
  Sequence = protein_sequences
)
# Extract the accession numbers using regular expressions
UniProt_protein[, Accession := sub(".*\\|(.*?)\\|.*", "\\1", ProteinName)]

## In silico digestion of proteins in UniProt
## create AAStringSet object
p <- AAStringSet(setNames(UniProt_protein$Sequence, UniProt_protein$Accession))
## cleavage miss zero cleavage position
UniProt_Tryptic <- cleave(p, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
UniProt_Tryptic <- unique(as.character(unlist(UniProt_Tryptic)))
UniProt_Tryptic <- UniProt_Tryptic[nchar(UniProt_Tryptic)>1]

## Read NCBI_RefSeq data
fasta_file <- ("NCBI_RefSeq_Human_2024_07_25.fasta")
# Read the FASTA file
sequences <- readAAStringSet(fasta_file) # use readAAStringSet for protein sequences
# Extract protein names and sequences
protein_names <- names(sequences)
protein_sequences <- as.character(sequences)
# Create a data table with protein names and sequences
NCBI_protein <- data.table(
  ProteinName = protein_names,
  Sequence = protein_sequences
)
# Extract the accession numbers using regular expressions
NCBI_protein[, Accession := sub("^(\\S+).*", "\\1", ProteinName)]

## In silico digestion of proteins in NCBI
## create AAStringSet object
p <- AAStringSet(setNames(NCBI_protein$Sequence, NCBI_protein$Accession))
## cleavage miss zero cleavage position
NCBI_Tryptic <- cleave(p, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
NCBI_Tryptic <- unique(as.character(unlist(NCBI_Tryptic)))
NCBI_Tryptic <- NCBI_Tryptic[nchar(NCBI_Tryptic)>1]

# save peptides of union of NCBI and UniProt
UniProtNCBI_Tryptic <- union(NCBI_Tryptic, UniProt_Tryptic)
save(UniProtNCBI_Tryptic, file = "UniProtNCBI_Tryptic.RData")





################################################################################
## In silico digestion of antibodies in OAS_full folder
# Create directories if they don't exist
dir.create("OAS_tryptic", showWarnings = FALSE)
dir.create("OAS_tryptic/None", showWarnings = FALSE)
dir.create("OAS_tryptic/Disease_index_ab", showWarnings = FALSE)
dir.create("OAS_tryptic/Disease_index_filename", showWarnings = FALSE)

# Load UniProtNCBI_Tryptic
load("UniProtNCBI_Tryptic.RData")

# Read metadata file
metadata <- fread("OAS_metadata.csv")

# Summarize the data: group by Disease and sum the Unique sequences
disease_summary <- metadata[, .(Total_unique_sequences = sum(`Unique sequences`)), by = Disease]

# Order by Total_unique_sequences in descending order
disease_summary <- disease_summary[order(-Total_unique_sequences)]

# Add the Disease_no column for diseases other than "None"
disease_summary[Disease != "None", Disease_no := paste0("D", seq_len(.N))]

# Keep Disease_no as "None" for the "None" Disease
disease_summary[Disease == "None", Disease_no := "None"]

# Filter metadata for files with Disease = "None" and with specific diseases
metadata_none <- metadata[Disease == "None"]
metadata_disease <- metadata[Disease != "None"]

# Get list of file paths for files with Disease = "None"
file_paths_none <- file.path("OAS_full", metadata_none$Filename)

# Get list of file paths for files with specific diseases
file_paths_disease <- file.path("OAS_full", metadata_disease$Filename)






################################################################################
# Section 1: Digest antibodies in Disease-specific files
# Each disease is one data with 2 columns: Sequence, Antibody
processed_diseases_file <- "Log_files/Part2_Digestion_processed_diseases1.txt"

# Load the list of already processed diseases, if the file exists
if (file.exists(processed_diseases_file)) {
  processed_diseases <- readLines(processed_diseases_file)
} else {
  processed_diseases <- character(0)
}

# Process files with specific diseases in groups corresponding to each disease
if (nrow(disease_summary) > 0) {
  for (disease in rev(disease_summary$Disease[disease_summary$Disease != "None"])) {
    
    # Check if this disease has already been processed
    if (disease %in% processed_diseases) {
      cat(sprintf("Skipping already processed disease: %s\n", disease))
      next
    }
    
    cat(sprintf("Processing files for disease: %s\n", disease))
    
    # Filter metadata for current disease
    disease_metadata <- metadata[Disease == disease]
    
    # Get file paths for the current disease
    disease_file_paths <- file.path("OAS_full", disease_metadata$Filename)
    
    if (length(disease_file_paths) > 0) {
      # Use rbindlist instead of do.call(rbind, ...)
      df_combined <- rbindlist(lapply(disease_file_paths, function(fp) {
        # Read the CSV file
        df <- read_csv_arrow(fp, skip = 0, col_select = "sequence_alignment_aa")
        # Convert all columns to character
        df <- df %>% mutate(across(everything(), as.character))
        return(df)
      }))
      
      # Initialize an empty data.table to accumulate all results
      table_pep_combined_all <- data.table()
      
      # Process in chunks of 5,000,000
      seq_ab_combined <- unique(df_combined$sequence_alignment_aa)
      total_chunks <- ceiling(length(seq_ab_combined) / 5e6) # 5,000,000 per chunk
      rm(df_combined)
      
      for (chunk in seq_len(total_chunks)) {
        cat(sprintf("     Processing chunk %d of %d\n", chunk, total_chunks))
        
        # Define chunk start and end
        start_index <- (chunk - 1) * 5e6 + 1
        end_index <- min(chunk * 5e6, length(seq_ab_combined))
        
        # Extract chunk of sequences
        seq_ab_chunk <- seq_ab_combined[start_index:end_index]
        names_ab_chunk <- paste0("ab_", seq(start_index, end_index))
        
        tryptic_combined <- AAStringSet(setNames(seq_ab_chunk, names_ab_chunk))
        tryptic_combined <- cleave(tryptic_combined, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
        seq_pep_combined <- as.character(unlist(tryptic_combined))
        seq_pep_name <- rep(names(tryptic_combined), lengths(tryptic_combined))
        
        table_pep_combined <- data.table("Sequence" = seq_pep_combined,
                                         "Antibody" = seq_pep_name)
        
        table_pep_combined <- unique(table_pep_combined[nchar(table_pep_combined$Sequence) > 3, ])
        table_pep_combined <- table_pep_combined[!(table_pep_combined$Sequence %in% UniProtNCBI_Tryptic), ]
        
        # Accumulate results in the combined data table
        table_pep_combined_all <- rbind(table_pep_combined_all, table_pep_combined)
        
        # Clean memory
        rm(list = c("table_pep_combined", "seq_pep_combined", "names_ab_chunk", "tryptic_combined"))
        gc()
      }
      
      # After all chunks are processed, table_pep_combined_all contains all data
      # Now write the accumulated data to the output file
      
      # Get Disease_no from disease_summary
      disease_no <- disease_summary[Disease == disease, Disease_no]
      
      # Replace special characters in the disease name with "_"
      disease_cleaned <- gsub("[^A-Za-z0-9]", "_", disease)
      
      # Create output file name
      output_file <- file.path("OAS_tryptic/Disease_index_ab", paste0(disease_no, "_", disease_cleaned, ".csv.gz"))
      fwrite(table_pep_combined_all, file = output_file, compress = "gzip")
      
      # Clean up the combined dataframe
      rm(table_pep_combined_all)
      gc()
      
      # Mark this disease as processed by adding it to the list and saving to file
      processed_diseases <- c(processed_diseases, disease)
      writeLines(processed_diseases, processed_diseases_file)
    }
  }
}






################################################################################
# Section 2: Digest antibodies in Disease-specific files
# Each disease is one data with 2 columns: Sequence, Filename
processed_diseases_file <- "Log_files/Part2_Digestion_processed_diseases2.txt"

# Load the list of already processed diseases, if the file exists
if (file.exists(processed_diseases_file)) {
  processed_diseases <- readLines(processed_diseases_file)
} else {
  processed_diseases <- character(0)
}

# Process files with specific diseases in groups corresponding to each disease
if (nrow(disease_summary) > 0) {
  for (disease in rev(disease_summary$Disease[disease_summary$Disease != "None"])) {
    
    # Check if this disease has already been processed
    if (disease %in% processed_diseases) {
      cat(sprintf("Skipping already processed disease: %s\n", disease))
      next
    }
    
    cat(sprintf("Processing files for disease: %s\n", disease))
    
    # Filter metadata for current disease
    disease_metadata <- metadata[Disease == disease]
    
    # Get file paths for the current disease
    disease_file_paths <- file.path("OAS_full", disease_metadata$Filename)
    
    if (length(disease_file_paths) > 0) {
      # Use rbindlist instead of do.call(rbind, ...)
      df_combined <- rbindlist(lapply(disease_file_paths, function(fp) {
        # Read the CSV file
        df <- read_csv_arrow(fp, skip = 0, col_select = "sequence_alignment_aa")
        df$Filename <- basename(fp)
        # Convert all columns to character
        df <- df %>% mutate(across(everything(), as.character))
        return(df)
      }))
      
      # Initialize an empty data.table to accumulate all results
      table_pep_combined_all <- data.table()
      
      # Process in chunks of 5,000,000
      total_chunks <- ceiling(nrow(df_combined) / 5e6) # 5,000,000 per chunk
      
      for (chunk in seq_len(total_chunks)) {
        cat(sprintf("     Processing chunk %d of %d\n", chunk, total_chunks))
        
        # Define chunk start and end
        start_index <- (chunk - 1) * 5e6 + 1
        end_index <- min(chunk * 5e6, nrow(df_combined))
        
        # Extract chunk of sequences
        seq_ab_chunk <- df_combined$sequence_alignment_aa[start_index:end_index]
        names_ab_chunk <- df_combined$Filename[start_index:end_index]
        
        # Digestion
        tryptic_combined <- AAStringSet(setNames(seq_ab_chunk, names_ab_chunk))
        tryptic_combined <- cleave(tryptic_combined, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
        seq_pep_combined <- as.character(unlist(tryptic_combined))
        seq_pep_name <- rep(names(tryptic_combined), lengths(tryptic_combined))
        
        # Create the data table with tryptic sequences
        table_pep_combined <- data.table("Sequence" = seq_pep_combined,
                                         "Filename" = seq_pep_name)
        
        table_pep_combined <- unique(table_pep_combined[nchar(table_pep_combined$Sequence) > 3, ])
        
        # Remove peptides overlapping with peptides in UniProt & NCBI
        table_pep_combined <- table_pep_combined[!(table_pep_combined$Sequence %in% UniProtNCBI_Tryptic), ]
        
        # Accumulate results in the combined data table
        table_pep_combined_all <- rbind(table_pep_combined_all, table_pep_combined)
        
        # Clean memory
        rm(list = c( "seq_ab_chunk", "names_ab_chunk", "tryptic_combined", 
                    "seq_pep_combined", "seq_pep_name", "table_pep_combined" ))
        gc()
      }
      
      # Clean memory
      rm(df_combined)
      
      # After all chunks are processed, table_pep_combined_all contains all data
      # Now write the accumulated data to the output file
      
      # Get Disease_no from disease_summary
      disease_no <- disease_summary[Disease == disease, Disease_no]
      
      # Replace special characters in the disease name with "_"
      disease_cleaned <- gsub("[^A-Za-z0-9]", "_", disease)
      
      # Create output file name
      output_file <- file.path("OAS_tryptic/Disease_index_filename", paste0(disease_no, "_", disease_cleaned, ".csv.gz"))
      fwrite(table_pep_combined_all, file = output_file, compress = "gzip")
      
      # Clean up the combined dataframe
      rm(table_pep_combined_all)
      gc()
      
      # Mark this disease as processed by adding it to the list and saving to file
      processed_diseases <- c(processed_diseases, disease)
      writeLines(processed_diseases, processed_diseases_file)
    }
  }
}









################################################################################
# Section 3: Digest antibodies in "None" disease files
# save them to files with same names as digested files

processed_none_file <- "Log_files/Part2_Digestion_processed_none_files.txt"

# Load the list of already processed "None" files, if the file exists
if (file.exists(processed_none_file)) {
  processed_none_files <- readLines(processed_none_file)
} else {
  processed_none_files <- character(0)
}

# Function to process each file
process_file <- function(file_path, index, total, output_dir) {
  cat(sprintf("Processing file %d of %d: %s\n", index, total, basename(file_path)))
  
  tryCatch({
    df <- read_csv_arrow(file_path, skip = 0)
    seq_ab <- df$sequence_alignment_aa
    
    # Update the names_ab to include "_ab" and the position number
    names_ab <- paste0(basename(file_path))
    
    tryptic <- AAStringSet(setNames(seq_ab, names_ab))
    tryptic <- cleave(tryptic, enzym="trypsin", missedCleavages=0:1, unique=TRUE)
    seq_pep <- as.character(unlist(tryptic))
    
    # Create the data table with tryptic sequences
    table_pep <- data.table("Sequence" = seq_pep)
    table_pep <- unique(table_pep[nchar(table_pep$Sequence) > 3, ])
    table_pep <- table_pep[!(table_pep$Sequence %in% UniProtNCBI_Tryptic), ]
    
    # Use the same file name for output as the input file
    output_file <- file.path(output_dir, basename(file_path))
    fwrite(table_pep, file = output_file, compress = "gzip")
    
  }, error = function(e) {
    cat(sprintf("Error processing file %d of %d: %s\n", index, total, basename(file_path)))
    cat("Error details: ", conditionMessage(e), "\n")
    failed_files <<- c(failed_files, file_path)
  })
}

# Process files with Disease = "None" individually
if (length(file_paths_none) > 0) {
  total_files_none <- length(file_paths_none)
  for (i in seq_along(file_paths_none)) {
    
    # Check if this file has already been processed
    if (basename(file_paths_none[i]) %in% processed_none_files) {
      cat(sprintf("Skipping already processed file: %s\n", basename(file_paths_none[i])))
      next
    }
    
    process_file(file_paths_none[i], i, total_files_none, "OAS_tryptic/None")
    
    # Mark this file as processed by adding it to the list and saving to file
    processed_none_files <- c(processed_none_files, basename(file_paths_none[i]))
    writeLines(processed_none_files, processed_none_file)
    
    # Clean memory
    rm(list = c("table_pep", "seq_pep", "names_ab", "tryptic", "df"))
    gc()
  }
}
















################################################################################
# Section 4: Save disease peptides to RData files

# Define the directory containing the .csv.gz files
directory_path <- "OAS_tryptic/Disease_index_ab"

# Define the directory to save the RData files
dir.create("OAS_tryptic/RData_peptides", showWarnings = FALSE)
output_directory <- "OAS_tryptic/RData_peptides"

# List all .csv.gz files in the directory
csv_files <- list.files(directory_path, pattern = "*.csv.gz", full.names = TRUE)

# Iterate through each file, read it, extract unique sequences, and save them if not already processed
for (file in rev(csv_files)) {
  # Extract the base name from the filename (e.g., D1_SARS_COV_2.csv.gz becomes D1_SARS_COV_2)
  file_base_name <- sub("\\.csv\\.gz$", "", basename(file))
  
  # Use the base file name directly as the variable name (e.g., D1_SARS_COV_2)
  var_name <- file_base_name
  
  # Create a filename for the output .RData file by replacing the .csv.gz extension with .RData
  output_filename <- file.path(output_directory, paste0(file_base_name, ".RData"))
  
  # Check if the RData file already exists
  if (file.exists(output_filename)) {
    cat(paste("File", output_filename, "already exists. Skipping analysis.\n"))
    next  # Skip to the next file
  }
  
  # message
  cat(paste("Processing file:", file, "\n"))
  
  # Read the compressed CSV file into a data.table
  data <- as.data.table(read_csv_arrow(file))
  
  # Extract the unique sequences
  unique_sequences <- unique(data$Sequence)
  
  # Assign the unique sequences to the dynamically created variable name
  assign(var_name, unique_sequences)
  
  # Save the dynamically created variable to an .RData file
  save(list = var_name, file = output_filename)
  
  # message
  cat(paste("   Finished!", "\n"))
}
