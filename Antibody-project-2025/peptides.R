# Load necessary packages
library(duckdb)
library(data.table)

# --- Function to process a single DuckDB file ---
process_duckdb_file <- function(duckdb_file) {
  # Connect to the DuckDB database
  con <- dbConnect(duckdb::duckdb(), dbdir = duckdb_file, read_only = TRUE)
  
  # Read data from the "DATDB" table
  dt <- as.data.table(dbGetQuery(con, "SELECT * FROM DATDB"))
  
  # Disconnect from the database
  dbDisconnect(con, shutdown = TRUE)
  
  # Extract disease name from the file path
  disease_name <- basename(duckdb_file)
  
  # Calculate unique sequences per Isotype
  summary_dt <- dt[, .(count = uniqueN(Sequence)), by = Isotype]
  
  # Reshape the data to be wide
  summary_dt_wide <- dcast(summary_dt, . ~ Isotype, value.var = "count", fill = 0)
  
  # Add the Disease column
  summary_dt_wide[, File_db := disease_name]
  
  # Remove the temporary '.' column created by dcast
  summary_dt_wide[, . := NULL]
  
  return(summary_dt_wide)
}

# --- Main script ---

# Define the path to the folder containing the DuckDB files
pepdb_folder <- "pepDB"

# Get a list of all .duckdb files in the folder
duckdb_files <- list.files(pepdb_folder, pattern = "\\.duckdb$", full.names = TRUE)

# Process all files and combine the results
all_summaries <- rbindlist(lapply(duckdb_files, process_duckdb_file), fill = TRUE)

# Set column order to have Disease first
setcolorder(all_summaries, c("File_db", setdiff(names(all_summaries), "File_db")))

# Print the final combined data table
all_summaries[is.na(all_summaries)] <- 0
print(all_summaries)
# Save all_summaries to csv file
write.csv(all_summaries, file = "Isotype_table.csv", row.names = F)

# --- Merge with metadata1.csv ---

# Read metadata1.csv into a data.table
#metadata_dt <- fread("metadata1.csv")

# Merge all_summaries with metadata_dt by File_db
#merged_dt <- merge(all_summaries, metadata_dt, by = "File_db", all.x = TRUE)

# Print the merged data table
#print(merged_dt)
# Save all_summaries to csv file
#write.csv(merged_dt, file = "Isotype_table.csv", row.names = F)
