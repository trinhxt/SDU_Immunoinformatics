library(svDialogs)
library(feather)
library(data.table)
library(dplyr)
library(disk.frame)
library(stringr)
library(tidyverse)
library(ggplot2)
library(gtools)
library(scales)
library(gridExtra)
library(openxlsx)
setup_disk.frame()
options(future.globals.maxSize = Inf)

SaveFolder <- dlg_dir(title = "Select Save folder:", filters = dlg_filters[c("All"), ])$res
setwd(SaveFolder)
# Load file OAS-SARS-COV-2-summary.csv
DataInfo <- fread("OAS-SARS-COV-2-summary.csv", header = T)
DataInfo$Filename <- strtrim(basename(DataInfo$DownloadLink), nchar(basename(DataInfo$DownloadLink))-7)

# Read digested OAS data. This is a large data (217 millions peptide sequences), RAM memory may not be sufficient so we use disk.frame package to load data into hard disk instead of RAM.
#OAS_digested <- csv_to_disk.frame(
#  "OAS-SARS-COV2_digested_2023-05-15.txt", 
#  outdir = file.path("D:/DwnlData/Diskframe/OAS-SARS-COV2_Antibodies_digested"),
#  nchunks = 60,
#  overwrite = T)

# If we loaded/saved data before, then we read the saved data:
OAS_digested <- disk.frame("D:/DwnlData/Diskframe/OAS-SARS-COV2_Antibodies_digested")

# Read UNIPROT human data
UNIPROTdata <- fread("UNIPROT_isoform_human_digested.txt", header = T)
print(paste("Number of peptides in UNIPROT = ", nrow(UNIPROTdata)))
print(paste("Number of DIFFERENT peptides in UNIPROT = ", length(unique(UNIPROTdata$Sequence))))

# Collect peptide sequences in OASdata_digested:
OAS1  <- OAS_digested %>% 
  select(Protein_Name, Sequence) %>% 
  collect
print(paste("Number of peptides in OAS-SARS-COV-2 = ", nrow(OAS1)))
print(paste("Number of DIFFERENT peptides in OAS-SARS-COV-2 = ", uniqueN(OAS1$Sequence)))


# Find common sequences in OAS and UNIPROT, and extract OAS data which are not in UNIPROT:
a <- Reduce(intersect, list(OAS1$Sequence, UNIPROTdata$Sequence))
OAS1 <- OAS1[!(OAS1$Sequence %in% a), ]
print(paste("Number of peptides in OAS-SARS-COV-2 and not in UNIPROT = ", nrow(OAS1)))
print(paste("Number of DIFFERENT peptides in OAS-SARS-COV-2 and not in UNIPROT = ", length(unique(OAS1$Sequence))))

# Add DataID column (indicate samples SRA - Data file name, extract all characters before "__"):
OAS1$DataID <- gsub("\\__.*","", OAS1$Protein_Name)

# Add AntibodyID column (Get characters after "__" and convert them to numbers)
OAS1$AntibodyID <- as.numeric(gsub(".*__","", OAS1$Protein_Name))

# Remove the Protein_Name column
OAS1 <- OAS1[,-c("Protein_Name")]; head(OAS1)

# Add BioProject and PatientID columns:
BioProject <- DataInfo %>% select("DS.Name", "Individual", "Filename")
colnames(BioProject) <- c("BioProject", "Individual", "DataID")
BioProject$PatientID <- paste(BioProject$BioProject, BioProject$Individual, sep = "_")
BioProject <- BioProject[,-c("Individual")]
# join 2 data to get BioProject and Patient info
OAS1 <- left_join(x = OAS1, y = BioProject, by = c("DataID"="DataID")); head(OAS1)

# Remove PatientID  having less than 100 peptides
PatientIDfreq <- as.data.table(table(OAS1$PatientID)); colnames(PatientIDfreq) <- c("PatientID", "Freq")
PatientList <- PatientIDfreq$PatientID[PatientIDfreq$Freq > 100]
OAS1  <- OAS1[OAS1$PatientID %in% PatientList,]

# Add v_call (variable gene) column:
OAS_Antibodies <- as.data.table(read_feather(dlg_open(title = "Select Save folder:", filters = dlg_filters[c("All"), ])$res))
#OAS_Antibodies <- as.data.table(read_feather("OAS-SARS-COV2_Antibodies_2023-05-15.feather"))
v_call <- as.data.table(OAS_Antibodies[,c("v_call", "AntibodyID")])
OAS1 <- as.data.table(OAS1)
OAS1 <- merge.data.table(OAS1, v_call, by = "AntibodyID"); OAS1
write_feather(OAS1, paste("OAS-SARS-COV2_peptides_non-UniProt_", Sys.Date(), ".feather", sep = ""))



#-------------------------------------------------------------------------------
# Get gene distribution (v_call)
Table1 <- as.data.table(table(OAS1$v_call)); colnames(Table1) <- c("Gene", "No.peptide.seqs")
Table1 <- Table1[order(Table1$No.peptide.seqs,decreasing=T),]
Table1$Percent <- round(Table1$No.peptide.seqs/sum(Table1$No.peptide.seqs)*100, 5)
Table1$Gene <- factor(Table1$Gene, levels = Table1$Gene)

# Save to csv
write.csv(Table1, file = "P4-OAS-SARS-COV2-Gene_dist.csv", row.names = F)

# Draw plot gene distribution
library(ggplot2)
Dataplot <- Table1
HistPlot <- ggplot(data = Dataplot, aes(x=Gene, y=No.peptide.seqs)) +
  geom_bar(position="stack", stat="identity", alpha=1.0, linewidth=0.02, width=0.75, color="orange", fill="orange") +
  xlab("Gene") + ylab("Number of peptide sequences") +
  guides(fill=guide_legend(ncol = 1, title=" ")) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5),
        legend.spacing.y = unit(20, "pt"),
        legend.position = "right") +
  guides(fill=guide_legend(title="Bio-project"))
HistPlot
# Save plot to pdf
pdf(file = paste(getwd(), "/P4-Gene_dist_all_non_UniProt_peptides_", Sys.Date(), ".pdf", sep = ""), width = 15, height = 6)
HistPlot
dev.off()

# Draw plot as percentage
HistPlot <- ggplot(data = Dataplot, aes(x=Gene, y=Percent)) +
  geom_bar(position="stack", stat="identity", alpha=1.0, linewidth=0.02, width=0.75, color="orange", fill="orange") +
  xlab("Gene") + ylab("Fraction (%)") +
  guides(fill=guide_legend(ncol = 1, title=" ")) +
  theme(axis.text.x=element_text(angle=90,vjust=0.5),
        legend.spacing.y = unit(20, "pt"),
        legend.position = "right") +
  guides(fill=guide_legend(title="Bio-project"))
HistPlot
# Save plot to pdf
pdf(file = paste(getwd(), "/P4-Gene_dist_all_non_UniProt_peptides_percent_", Sys.Date(), ".pdf", sep = ""), width = 15, height = 6)
HistPlot
dev.off()
