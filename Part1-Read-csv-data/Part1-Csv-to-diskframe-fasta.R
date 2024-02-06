library(svDialogs)
library(dplyr)
library(disk.frame)
library(data.table)
library(feather)
library(fst)
library(readxl)
setup_disk.frame()
options(future.globals.maxSize = Inf)

# Choose data folder
DataFolder <- dlg_dir(title = "Select csv data folder:", filters = dlg_filters[c("All"), ])$res

# Get list of data files
DataFiles <- list.files(path = DataFolder, pattern = ".csv.gz", all.files = FALSE, 
                        full.names = TRUE, recursive = TRUE, ignore.case = FALSE, 
                        include.dirs = FALSE, no.. = FALSE)

# Read metadata file (OAS-SARS-COV-2-summary.csv) and add columns Filename and PatientID
DataSummary <- read.csv(file = "OAS-SARS-COV-2-summary.csv", header = T)
DataSummary$Filename <- basename(DataSummary$DownloadLink)
DataSummary$PatientID <- paste(DataSummary$DS.Name, DataSummary$Individual, sep = "_")

# Set folder for saving data
setwd(dlg_dir(title = "Select where to save data:", filters = dlg_filters[c("All"), ])$res)


# Assign a temporary folder to save disk frame data
tempFolder <- dlg_dir(title = "Select temporory folder to save disk frame:", filters = dlg_filters[c("All"), ])$res
setwd(tempFolder)

# Read each csv file, save them as a list of diskframes to tempFolder
 Listdata <- lapply(c(1:length(DataFiles)), function(i){
 dat <- csv_to_disk.frame(
    DataFilesFilter[i], 
    outdir = file.path(paste(tempFolder, "/tmp_",i, sep = "")),
    nchunks = 60,
    overwrite = T)
})

# Row-bind list of disk.frames together into one big diskframe
 DatatableDisk <- rbindlist.disk.frame(Listdata,
                                  outdir = (paste(tempFolder, "/OAS-SARS-COV2_Antibodies", sep = "")),
                                  overwrite = T,
                                  by_chunk_id = TRUE,
                                  parallel = TRUE)

# Load OAS data that was saved in hard disk:
OAS_Abs <- disk.frame(dlg_dir(title = "Select disk frame data folder:", filters = dlg_filters[c("All"), ])$res)
# Check column names and row number of OAS_Abs
OAS_colnames <- colnames(OAS_Abs); OAS_colnames
nrow(OAS_Abs)

# Extract some important columns in OAS_Abs (because the whole OAS_Abs data is very large to read in RAM, so we just take some columns)
Datatable  <- OAS_Abs %>% select(sequence_alignment_aa,
                                 v_call,
                                 d_call,
                                 j_call,
                                 fwr1_aa,
                                 cdr1_aa,
                                 fwr2_aa,
                                 cdr2_aa,
                                 fwr3_aa,
                                 cdr3_aa,
                                 fwr4_aa) %>%
                                collect
# Add AntibodyID column
Datatable$AntibodyID <- c(1:nrow(OAS_Abs))

# Make a vector of Data file names that match the 
Datafilename <- lapply(c(1:length(DataFiles)), function(i){
  rep(basename(DataFiles[i]), DataSummary$Unique.Sequences[DataSummary$Filename==basename(DataFiles[i])])
})
# Add data file name to Datatable of antibodies
Datatable$Datafilename <- unlist(Datafilename)

# Correct v_call column: remove any character after *
Datatable$v_call <- gsub("\\*.*","", Datatable$v_call)

# Save data to feather file
write_feather(Datatable, paste(dirname(DataFolder),"/OAS-SARS-COV2_Antibodies_", Sys.Date(), ".feather", sep = ""))

# Load feather datafile
# Datatable <- read_feather(dlg_open(title = "Select feather file:", filters = dlg_filters[c("All"), ])$res, columns = NULL)



#-------------------------------------------------------------------------------
# Save antibody sequences to fasta file

# Extract antibody sequences
Sequences <- as.list(c(Datatable$sequence_alignment_aa))

# Get data file names (without extension)
NameList <- strtrim(Datatable$Datafilename, nchar(as.character(Datatable$Datafilename))-7)

# Make sequence names
Nameseq <- c(paste(NameList, Datatable$AntibodyID, sep = "__"))
                
# Save data to fasta file
library(seqinr)
write.fasta(Sequences,   # sequences
            Nameseq, # name of sequences
            paste(dirname(DataFolder),"/OAS-SARS-COV2_Antibodies_", Sys.Date(), ".fasta", sep = ""), # location to save fasta file  
            open = "w", nbchar = 60, as.string = TRUE)

