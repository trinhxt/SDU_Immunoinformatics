#-------------------------------------------------------------------------------
# Set a custom library path (optional)
custom_lib <- "/work/Immunoinformatics/R-packages-RStudio/" # Change this to your directory
#custom_lib <- "/work/Immunoinformatics/R-packages-Ubuntu/"
.libPaths(custom_lib)

# Load all packages
cran_packages <- c("arrow", "dplyr", "jsonlite", "seqinr", "data.table", "fs", "stringr", "stringi", "parallel", "DBI", "RSQLite", "tools")
bioc_packages <- c("cleaver", "Biostrings", "GenomeInfoDb")
lapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

################################################################################
# Read URLs from the shell script file
shell_script_file <- "bulk_download.sh"
lines <- readLines(shell_script_file)
urls <- sub("wget ", "", lines)

# Total number of files
total_files <- length(urls)

# Initialize an empty list to store each arrow table and metadata
metadata_list <- list()

# Create "OAS_full" folder if it does not exist
output_folder <- "OAS_full"
if (!dir.exists(output_folder)) {
  dir.create(output_folder)
}

# Function to read and parse the first line of a .csv.gz file
read_metadata <- function(file) {
  first_line <- readLines(gzfile(file), n = 1)
  # Preprocess to replace NaN with null
  first_line <- gsub('\\bNaN\\b', 'null', first_line)
  # Remove leading and trailing quotes, and replace double double-quotes with single double-quotes
  clean_json <- gsub('^"|"$', '', first_line)
  clean_json <- gsub('""', '"', clean_json)
  fromJSON(clean_json)
}

# Function to process each file
process_file <- function(url, index, total) {
  cat(sprintf("Processing file %d of %d\n", index, total))
  
  # Download the file to the current directory
  file_name <- basename(url)
  dest_file <- file.path(file_name)
  
  # Check if the file already exists
  if (!file.exists(dest_file)) {
    tryCatch({
      download.file(url, destfile = dest_file)
    }, error = function(e) {
      cat(sprintf("Error downloading file %d: %s\n", index, e$message))
      return(NULL)  # Return NULL if download fails
    })
  }
  
  # If the file exists, process it
  if (file.exists(dest_file)) {
    tryCatch({
      # Read and parse the first line as JSON metadata
      metadata <- read_metadata(dest_file)
      
      # Convert metadata to a data table
      metadata_df <- as.data.table(t(metadata), stringsAsFactors = FALSE)
      # Add file name
      metadata_df$Filename <- file_name
      # Add metadata_df to metadata_list
      metadata_list <<- append(metadata_list, list(metadata_df))
      
      # Extract data and select specific columns
      df <- read_csv_arrow(dest_file, skip = 1)
      selected_cols <- df %>% select(sequence_alignment_aa, v_call, d_call, j_call, cdr1_aa, cdr2_aa, cdr3_aa)
      
      # Save the selected columns to OAS_full folder
      output_file <- file.path(output_folder, file_name)
      fwrite(selected_cols, file = output_file, compress = "gzip")
      
      # Delete the original file after processing
      file.remove(dest_file)
    }, error = function(e) {
      cat(sprintf("Error processing file %d: %s\n", index, e$message))
      # Optionally, you can choose to keep the file for debugging
    })
  }
}

# Process each file in the list with an indication of progress
lapply(seq_along(urls), function(i) process_file(urls[i], i, total_files))

# Combine all metadata data frames into a single data frame
metadata_df <- rbindlist(metadata_list, fill = TRUE)

# Identify unprocessed files
processed_files <- metadata_df$Filename
unprocessed_urls <- urls[!basename(urls) %in% processed_files]

cat(sprintf("Found %d unprocessed files. Re-processing...\n", length(unprocessed_urls)))

# Re-process unprocessed files
lapply(seq_along(unprocessed_urls), function(i) process_file(unprocessed_urls[i], i, length(unprocessed_urls)))

# Combine all metadata data frames into a single data frame again after reprocessing
metadata_df <- rbindlist(metadata_list, fill = TRUE)

# Combine Author and Subject into a single identifier
metadata_df[, Patient := paste(Author, Subject, sep = "_")]

# Create a unique numeric identifier for each unique combination
metadata_df[, Patient := paste0("P", as.numeric(factor(Patient)))]

# Save metadata
# If there are any list columns, convert them to character or another supported type.
metadata_df[] <- lapply(metadata_df, function(col) {
  if (is.list(col)) {
    return(sapply(col, toString))
  } else {
    return(col)
  }
})

write.csv(metadata_df, "OAS_metadata.csv", row.names = FALSE)
