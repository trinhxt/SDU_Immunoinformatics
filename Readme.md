Data mining antibody sequences for database searching in bottom-up
proteomics
================

## 0. Introduction

Bottom-up proteomics approaches rely on database searches that compare
experimental values of peptides to theoretical values derived from
protein sequences in a database. While the human body can produce
millions of distinct antibodies, current databases for human antibodies
such as UniProtKB are limited to only 1095 sequences (as of 2024
January). This limitation may hinder the identification of new
antibodies using bottom-up proteomics. Therefore, extending the
databases is an important task for discovering new antibodies.

Herein, we adopted extensive collection of antibody sequences from
[Observed Antibody Space](https://opig.stats.ox.ac.uk/webapps/oas/) for
conducting efficient database searches in publicly available proteomics
data with a focus on the SARS-CoV-2 disease. Thirty million heavy
antibody sequences from 146 SARS-CoV-2 patients in the [Observed
Antibody Space](https://opig.stats.ox.ac.uk/webapps/oas/) were *in
silico* digested to obtain 18 million unique peptides. These peptides
were then used to create six databases (DB1-DB6) for bottom-up
proteomics. We used those databases for searching antibody peptides in
publicly available SARS-CoV-2 human plasma samples in the Proteomics
Identification Database (PRIDE), and we consistently found new antibody
peptides in those samples. The database searching task was done by using
[Fragpipe](https://fragpipe.nesvilab.org/) softwares.

Remaining content of this document is for creating those databases.

## 1. Downloading antibody sequence data

The antibody sequences (heavy chains) were downloaded from The Observed
Antibody Space database
([OAS](https://opig.stats.ox.ac.uk/webapps/oas/oas_unpaired/)). In this
study, we focus on antibodies of SARS-CoV-2, so we choose keywords as
below:

| Keyword      | Note                                                                             | Value      |
|--------------|----------------------------------------------------------------------------------|------------|
| Species      | Species being sequenced                                                          | human      |
| BSource      | Organ or tissue from where the B-cells were collected                            | PBMC       |
| BType        | Type of B-cells                                                                  | \*         |
| Longitudinal | Indicating whether the study was carried out at several time points              | undefined  |
| Age          | Indicating whether the age of the individual was recorded                        | \*         |
| Disease      | Indicating whether the individual was diseased at the time of sequencing         | SARS-CoV-2 |
| Subject      | Indicating whether the sequences can be traced from particular individual        | defined    |
| Vaccine      | Indicating whether the individual was vaccinated                                 | None       |
| Chain        | Specifying heavy or light chain                                                  | Heavy      |
| Isotype      | Isotype of heavy chain if available (IgG, IgM, etc.), Bulk stands for unassigned | \*         |

Table 1. Keywords search in OAS database for SARS-CoV-2.

After choosing keywords and clicking Search button, the website will
give us results: *Your search yielded 30,966,193 unique sequences from 9
studies… A shell-script with the commands to download all the data-units
in this subset of OAS can be downloaded
[here](blob:https://opig.stats.ox.ac.uk/e3480b42-4861-49e6-9b5d-e1852b32baa1).*
We downloaded the shell-script file (`bulk_download.sh`) and saved it to
the directory `Part0-Download-data`.

In order to download all data files, we use terminal (in Mac and Linux,
in Windows we can use terminal of [WSL - Windows Subsystem for
Linux](https://learn.microsoft.com/en-us/windows/wsl/)), navigate the
terminal to where the `bulk_download.sh` file was saved and run the
command:

``` bash
bash bulk_download.sh
```

There are total 990 data file and the size of all data is 13 GB.

## 2. Read antibody sequences data

In this part, we will read the downloaded data (`csv` format) and save
it as `fasta` format (for later *in silico* digestion) and `feather`
format (for later load and use). Full R code and metadata file for this
part are saved in the folder `Part1-Read-csv-data` of this Github
repository.

We will read 990 data files (`csv` format )with the size of 13GB. This
task requires at least 32GB of RAM by using `fread` function of the
`data.table` package. For computers with less than 32GB of RAM, we
recommend to use `diskframe` package as below (the data will be read
from and saved to the hard disk, bypassing storage in RAM to prevent
overloading the computer’s memory). After reading the whole csv data, we
save it as feather format (feather format is designed to make reading
and writing big data frames efficient).

``` r
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
```

Save antibody sequences to `fasta` format:

``` r
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
```

## 3. *In silico* digestion of antibody sequences
