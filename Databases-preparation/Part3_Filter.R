#-------------------------------------------------------------------------------
# Set a custom library path
custom_lib <- "/work/Immunoinformatics/R-packages-RStudio/"
#custom_lib <- "/work/Immunoinformatics/R-packages-Ubuntu/"
.libPaths(custom_lib)

# Load all packages
cran_packages <- c("arrow", "dplyr", "jsonlite", "seqinr", "data.table", "fs", "stringr", "stringi", "parallel", "DBI", "RSQLite", "tools")
bioc_packages <- c("cleaver", "Biostrings", "GenomeInfoDb")
lapply(c(cran_packages, bioc_packages), require, character.only = TRUE)

################################################################################
# Define the working folder and the log folder
working_folder <- "OAS_tryptic"
log_folder <- "Log_files"
progress_file_loop1 <- file.path(log_folder, "Part3_Filter_Loop1.json")
progress_file_loop2 <- file.path(log_folder, "Part3_Filter_Loop2.json")

# Load progress from previous run if exists for Loop 1
if (file.exists(progress_file_loop1)) {
  completed_diseases_loop1 <- fromJSON(progress_file_loop1)
} else {
  completed_diseases_loop1 <- c()
}

# Load progress from previous run if exists for Loop 2
if (file.exists(progress_file_loop2)) {
  completed_diseases_loop2 <- fromJSON(progress_file_loop2)
} else {
  completed_diseases_loop2 <- c()
}

################################################################################
# Get the list of .csv.gz files in the "Disease" folder within the working folder
disease_files <- dir_ls(file.path(working_folder, "Disease_index_ab"), glob = "*.csv.gz")
# Extract the numeric part of the filenames
disease_files_order <- as.numeric(str_extract(basename(disease_files), "(?<=D)\\d+"))
# Order the files based on this numeric part in reverse order
ordered_indices <- order(disease_files_order, decreasing = TRUE)
disease_files <- disease_files[ordered_indices]

################################################################################
# Get the list of .csv.gz files in the "None" folder within the working folder
disease_none_files <- dir_ls(file.path(working_folder, "None"), glob = "*.csv.gz")

# Calculate the sizes of each "None" file
disease_none_file_sizes <- file_info(disease_none_files)$size

# Combine filenames with their sizes into a data frame
file_size_info <- data.frame(file = disease_none_files, size = disease_none_file_sizes)

# Sort the files by size (descending)
file_size_info <- file_size_info[order(file_size_info$size, decreasing = TRUE), ]

# Initialize chunk list and cumulative size tracker
disease_none_chunks <- list()
current_chunk <- list()
current_chunk_size <- 0
target_chunk_size <- sum(file_size_info$size) / ceiling(length(file_size_info$file) / 200)  # Adjust 200 based on your system constraints

# Distribute files into chunks with approximately equal total sizes
for (i in 1:nrow(file_size_info)) {
  file <- file_size_info$file[i]
  file_size <- file_size_info$size[i]
  
  if (current_chunk_size + file_size > target_chunk_size && length(current_chunk) > 0) {
    # Save the current chunk
    disease_none_chunks <- append(disease_none_chunks, list(current_chunk))
    # Reset for the new chunk
    current_chunk <- list(file)
    current_chunk_size <- file_size
  } else {
    # Add file to the current chunk
    current_chunk <- append(current_chunk, file)
    current_chunk_size <- current_chunk_size + file_size
  }
}

# Append the last chunk if it has any files
if (length(current_chunk) > 0) {
  disease_none_chunks <- append(disease_none_chunks, list(current_chunk))
}

################################################################################
# Process each file in disease_files_to_process
# Limit the disease_files to the specified range for this machine
disease_files_to_process <- disease_files








################################################################################
# Parent Loop 1: Remove overlapping sequences with selected .RData files

# Load all .RData files in the working folder and store loaded object names
oas_files <- dir_ls("OAS_tryptic/RData_peptides", glob = "*.RData")
loaded_objects <- list()  # To keep track of objects loaded from each .RData file

# Load each .RData file and store the names of the loaded objects
for (oas_file in oas_files) {
  tryCatch({
    load(oas_file)
    oas_variable_name <- tools::file_path_sans_ext(basename(oas_file)) # Get variable name from file name
    loaded_objects[[oas_variable_name]] <- oas_variable_name  # Store the loaded object name
    cat(paste("Loaded OAS file:", oas_file, "\n"))
  }, error = function(e) {
    message(paste("Error loading OAS file", basename(oas_file), "-", e$message))
  })
}



################################################################################
for (file in disease_files_to_process) {
  # Check if this file has already been processed in Loop 1
  if (basename(file) %in% completed_diseases_loop1) {
    message(paste("Skipping already processed file in Loop 1:", basename(file)))
    next
  }
  
  message(paste("Processing file for Loop 1:", basename(file)))
  
  # Read and filter the sequences within the current disease file with error handling
  disease_data <- tryCatch({
    read_csv_arrow(file)
  }, error = function(e) {
    message(paste("Ignoring file", basename(file), "due to error reading -", e$message))
    next
  })
  
  # Extract the base name of the current disease file without extension
  disease_file_basename <- sub("\\.csv\\.gz$", "", basename(file))
  
  # Filter oas_files to exclude files with the same base name as the current disease file
  oas_files_to_use <- oas_files[tools::file_path_sans_ext(basename(oas_files)) != disease_file_basename]
  
  # Filter sequences using only the relevant loaded objects
  for (oas_file in oas_files_to_use) {
    oas_var_name <- tools::file_path_sans_ext(basename(oas_file))
    
    if (oas_var_name %in% names(loaded_objects)) {
      message(paste("     Filtering", basename(file), "by comparing with", oas_var_name))
      
      tryCatch({
        # Filter the disease data to remove overlapping sequences
        disease_data <- disease_data %>%
          filter(!(Sequence %in% get(oas_var_name)))
        
      }, error = function(e) {
        message(paste("Error during filtering with", oas_var_name, "-", e$message))
      })
      
      # Clean up
      gc()
    }
  }
  
  # Save the filtered data back to the original location
  fwrite(disease_data, file, compress = "gzip")
  
  # Add the current file to the list of completed diseases for Loop 1
  completed_diseases_loop1 <- c(completed_diseases_loop1, basename(file))
  
  # Save the progress to the Loop 1 progress file
  write_json(completed_diseases_loop1, progress_file_loop1)
  
  # Remove disease_data to save memory
  rm(disease_data); gc()
  
  message(paste("Finished Loop 1 for:", basename(file)))

}

################################################################################
# Unload (remove) all objects loaded from .RData files
for (oas_var_name in names(loaded_objects)) {
  rm(list = loaded_objects[[oas_var_name]], envir = .GlobalEnv)
  cat(paste("Unloaded object:", oas_var_name, "\n"))
}
message("Finished loop 1!")  # Additional message indicating the completion of filtering
# Run garbage collection to free up memory
gc()











################################################################################
# Parent Loop 2: Remove overlapping sequences with disease_none files in balanced chunks
for (file in disease_files_to_process) {
  # Check if this file has already been processed in Loop 2
  if (basename(file) %in% completed_diseases_loop2) {
    message(paste("Skipping already processed file in Loop 2:", basename(file)))
    next
  }
  
  message(paste("Processing file for Loop 2:", basename(file)))
  
  # Read the disease file data again since it might have been modified in Loop 1
  disease_data <- tryCatch({
    read_csv_arrow(file)
  }, error = function(e) {
    message(paste("Ignoring file", basename(file), "due to error reading -", e$message))
    next
  })
  
  for (i in seq_along(disease_none_chunks)) {
    message(paste("Filtering", basename(file), "by comparing with chunk", i, "of", length(disease_none_chunks), "from none-disease data"))
    
    # Read all files in the current chunk at once with error handling
    none_data_list <- lapply(disease_none_chunks[[i]], function(none_file) {
      tryCatch({
        as.data.table(read_csv_arrow(none_file, col_select = "Sequence"))
      }, error = function(e) {
        message(paste("Ignoring file", basename(none_file), "due to error reading -", e$message))
        return(NULL)
      })
    })
    
    # Remove any NULLs from the list (files that failed to read)
    none_data_list <- none_data_list[!sapply(none_data_list, is.null)]
    
    # Combine all none_data into a single data frame using rbindlist from data.table
    combined_none_data <- rbindlist(none_data_list)
    unique_none_sequences <- unique(combined_none_data$Sequence)
    
    # Remove combined_none_data and none_data_list to save memory
    rm(combined_none_data, none_data_list); gc()
    
    # Filter disease_data using the unique sequences from all files in the current chunk
    disease_data <- disease_data %>%
      filter(!(Sequence %in% unique_none_sequences))
    
    # Remove unique_none_sequences to save memory
    rm(unique_none_sequences); gc()
  }
  
  # Save the filtered data back to the original location
  fwrite(disease_data, file, compress = "gzip")
  
  # Add the current file to the list of completed diseases for Loop 2
  completed_diseases_loop2 <- c(completed_diseases_loop2, basename(file))
  
  # Save the progress to the Loop 2 progress file
  write_json(completed_diseases_loop2, progress_file_loop2)
  
  # Remove disease_data to save memory
  rm(disease_data); gc()
  
  message(paste("Finished Loop 2 for:", basename(file)))
  
}

message("Finished loop2!")  # Additional message indicating the completion of filtering
message("All files finished!")
