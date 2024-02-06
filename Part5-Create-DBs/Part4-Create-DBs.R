library(svDialogs)
library(openxlsx)
library(readxl)
library(data.table)
library(phylotools)
library(stringr)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(feather)
library(seqinr)
library(stringi)



#-------------------------------------------------------------------------------
SaveFolder <- dlg_dir(title = "Select Save folder:", filters = dlg_filters[c("All"), ])$res
setwd(SaveFolder)

# Read summary file and columns Filename and PatientID
DataSummary <- read.csv(file = "OAS-SARS-COV-2-summary.csv", header = T)
DataSummary$Filename <- basename(DataSummary$DownloadLink)
DataSummary$PatientID <- paste(DataSummary$DS.Name, DataSummary$Individual, sep = "_")

# Read OAS antibody data 
#OAS_antibodies <- as.data.table(read_feather("OAS-SARS-COV2_Antibodies_2023-05-15.feather")); OAS_antibodies

# Read filtered OAS peptide data 
OAS_pep <- as.data.table(read_feather("OAS-SARS-COV2_peptides_non-UniProt_2023-05-16.feather")); OAS_pep

# Total number of antibodies
N_Abs <- uniqueN(OAS_pep$AntibodyID); N_Abs

# Get prequency of peptides in OAS_pep
Table1 <- OAS_pep[, .N , by = Sequence]
colnames(Table1) <- c("Sequence", "No.Abs")
Table1 <- Table1[order(Table1$No.Abs, decreasing=T), ]
Table1$Percent <- round(Table1$No.Abs/N_Abs*100, 6); Table1


Top_pep <- c(100, 10^3, 10^4, 10^5, 10^6, 10^7)


# Check antibody coverage for top peptides
No_Abs_DBs <- sapply(c(1:length(Top_pep)), function(i){
  uniqueN(OAS_pep$AntibodyID[ OAS_pep$Sequence %in% Table1$Sequence[c(1:Top_pep[i])] ] )
})

Table2 <- data.table("DBs"        = c("DB1", "DB2", "DB3", "DB4", "DB5", "DB6"),
                     "OAS_pep_no" = Top_pep, 
                     "Abs_no"     = No_Abs_DBs,
                     "Abs_pct"    = round(No_Abs_DBs/N_Abs*100,3))

write.csv(Table2, file = "P5-DBs-antibody-coverage.csv", row.names = F)





# Extract top peptides and Write fasta files
for (i in c(1:length(Top_pep))) {
  DB <- Table1[c(1:Top_pep[i]),]
  
  # Save peptide sequences to fasta file
  Sequences <- DB$Sequence
  # Get list of peptides and replace I to L (Isoleucine and Leucine have same mass)
  Sequences <- as.list(gsub("I", "L", c(Sequences)))
  
  NameList <- paste("Pep_", c(1:length(Sequences)), sep = "")
  
  Nameseq <- c(paste("sp|", "OAS_antibody_pep_", c(1:length(Sequences)), "|", sep = ""))
  
  # Save data to fasta file
  write.fasta(Sequences,   # sequences
              Nameseq, # name of sequences
              paste("OAS_",i,"_", Sys.Date(), ".fasta", sep = ""), # location to save fasta file  
              open = "w", nbchar = 60, as.string = TRUE)
  
}


#-------------------------------------------------------------------------------
# Combine fasta files (our peptide list and contaminats) together
# Go to Terminal, go to folder containing fasta files, run below command
cat *.fasta > DB6.fasta

