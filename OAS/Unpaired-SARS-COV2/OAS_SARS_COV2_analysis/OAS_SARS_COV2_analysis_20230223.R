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
setup_disk.frame()
options(future.globals.maxSize = Inf)
# Set working directory, which is folder containing digested data
setwd(dlg_dir(title = "Select working folder", filters = dlg_filters[c("All"), ])$res)

# Data infor
DataInfo <- fread("OAS-SARS-COV-2-summary.csv", header = T)
DataInfo$Filename <- strtrim(basename(DataInfo$DownloadLink), nchar(basename(DataInfo$DownloadLink))-7)

# Set folder of saving disk frames
tempFolder <- "D:/DwnlData/Diskframe"

# Read digested OAS data
#OASdata_digested <- csv_to_disk.frame(
#  "OAS_SARS_COV2_digested.txt", 
#  outdir = file.path(paste(tempFolder, "/OASdata_digested", sep = "")),
#  nchunks = 60,
#  overwrite = T)
OASdata_digested <- disk.frame("D:/DwnlData/Diskframe/OASdata_digested")
nrow(OASdata_digested)

# Read UNIPROT human data
UNIPROTdata <- fread("UNIPROT_human_digested.txt", header = T)

# Read un-digested OAS data
#OASdata <- csv_to_disk.frame(
#  "OAS_SARS_COV2.txt", 
#  outdir = file.path(paste(tempFolder, "/OASdata", sep = "")),
#  nchunks = 60,
#  overwrite = T)
OASdata <- disk.frame("D:/DwnlData/Diskframe/OASdata")
nrow(OASdata)

# Collect peptide sequences in OASdata_digested
OAS1  <- OASdata_digested %>% 
  select(Protein_Name, Sequence) %>% 
  collect
# number of peptides
nrow(OAS1)
# number of different peptides
length(unique(OAS1$Sequence))


#-------------------------------------------------------------------------------
# Find common sequences in OAS and UNIPROT
a <- Reduce(intersect, list(OAS1$Sequence, UNIPROTdata$Sequence))

# Extract OAS data which are not in UNIPROT
OAS2 <- OAS1[!(OAS1$Sequence %in% a), ]
# number of peptides not in UNIPROT
nrow(OAS2)
# number of different peptides not in UNIPROT
length(unique(OAS2$Sequence))
# number of antibodies that above peptide identify
length(unique(OAS2$AntibodyID))
# Add DataID column (indicate bioprojects)
OAS2$DataID <- gsub("\\__.*","", OAS2$Protein_Name)

# Add AntibodyID column
OAS2$AntibodyID <- as.numeric(gsub(".*__","", OAS2$Protein_Name))

OAS2 <- select(OAS2, -c("Protein_Name"))

# Add name of bio-project
BioProject <- DataInfo %>%
                select("DS.Name", "Individual", "Filename")

colnames(BioProject) <- c("BioProject", "Individual", "DataID")

# join 2 data to get BioProject and Patient info
OAS2 <- left_join(x = OAS2, y = BioProject, by = c("DataID"="DataID"))
# Add PatientID
OAS2$PatientID <- paste(OAS2$BioProject, OAS2$Individual, sep = "_")
# Remove PatientID data having frequency less than 1000
PatientIDfreq <- as.data.table(table(OAS2$PatientID)); colnames(PatientIDfreq) <- c("PatientID", "Freq")
PatientList <- PatientIDfreq$PatientID[PatientIDfreq$Freq > 1000]
OAS2  <- OAS2[OAS2$PatientID %in% PatientList,]

head(OAS2)



#-------------------------------------------------------------------------------
# Extract peptides that appear in all patientID
b <- Reduce(intersect, split(OAS2$Sequence, OAS2$PatientID))

OAS3 <- OAS2[OAS2$Sequence %in% b,]

# number of peptides shared between patientID
nrow(OAS3)
# number of antibodies that above peptides identify
length(unique(OAS3$AntibodyID))
# number of different peptides shared between patientID
length(b)

# Find frequency of each peptide in every PatientID
Freq_PatientID <- dcast(OAS3, PatientID~Sequence, length)
rownames(Freq_PatientID) <- Freq_PatientID$PatientID
Sequence <- colnames(Freq_PatientID)[2:ncol(Freq_PatientID)]

Freq_PatientID <- data.table::transpose(Freq_PatientID)

colnames(Freq_PatientID) <- as.character(Freq_PatientID[1,])
Freq_PatientID <- Freq_PatientID[-c(1),]
Freq_PatientID[, names(Freq_PatientID) := lapply(.SD, as.numeric)]

Freq_PatientID <- cbind(Sequence, Freq_PatientID)


OAS_shared_PatientID <- as.data.table(table(OAS3$Sequence)); colnames(OAS_shared_PatientID) <- c("Sequence", "Freq")
OAS_shared_PatientID <- OAS_shared_PatientID[order(-OAS_shared_PatientID$Freq),]
OAS_shared_PatientID <- left_join(x = OAS_shared_PatientID, y = Freq_PatientID, by = "Sequence")


# Convert frequency to percentage table
# Get total number of antibodies
Total.Abs <- length(unique(OAS3$AntibodyID))
# Get number of antibodies in each PatientID
PatientID <- colnames(OAS_shared_PatientID)[-c(1:2)]
No.Abs.each.PatientID <- sapply(c(1:length(PatientID)), function(i){
  length(unique(OAS3$AntibodyID[OAS3$PatientID == PatientID[i]]))
})

No.Abs <- c(Total.Abs, No.Abs.each.PatientID)

Freq <- as.matrix(OAS_shared_PatientID[,-c(1)])
Percentage <- round(100*sweep(Freq, 2, No.Abs, `/`), digits = 2)

OAS_shared_PatientID_Percent <- cbind(OAS_shared_PatientID[,c(1)], Percentage)

# save Frequency data by PatientID to excel file
library(openxlsx)
ExcelFile <- createWorkbook("Freq_Percent")
addWorksheet(ExcelFile, "Frequency")
writeData(ExcelFile, sheet = 1, OAS_shared_PatientID)
addWorksheet(ExcelFile, "Percent")
writeData(ExcelFile, sheet = 2, OAS_shared_PatientID_Percent)
saveWorkbook(ExcelFile, paste(getwd(), "/Freq_Table_by_PatientID_", Sys.Date(), ".xlsx", sep = ""), overwrite = TRUE)


#-------------------------------------------------------------------------------
# Draw histogram showing percentage appearance of peptides in antibodies

Hist_data <- data.table(SeqNo =  c(1:nrow(OAS_shared_PatientID_Percent)),
                        Percent = OAS_shared_PatientID_Percent$Freq,
                        SeqLength = nchar(OAS_shared_PatientID_Percent$Sequence))

HistPlot <- ggplot(data = Hist_data, aes(x=SeqNo, y=Percent, fill = factor(SeqLength), width=.75)) +
  geom_bar(stat="identity") +
  xlab("Peptide") + ylab("Appear in how many percent of Abs") +
  guides(fill=guide_legend(ncol = 2, title="Peptide \n length")) +
  theme(legend.position = "right") 
HistPlot


# Draw heat map of the OAS_shared_PatientID for seeing peptide sequences percentage in each PatientID
PopTable_plot <- OAS_shared_PatientID_Percent[,-c(1:2)]
colnames(PopTable_plot) <- as.character(c(1:ncol(PopTable_plot)))

PopTable_plot <- PopTable_plot %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
PopTable_plot$rowname <- as.numeric(PopTable_plot$rowname)

PopPlot <- ggplot(PopTable_plot, aes(x = rowname, y = as.numeric(colname), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="grey", high="red") +
  #scale_fill_gradientn(colours=c("grey", "red", "red", "red", "red"),
  #                                values=rescale(c(0, 0.1, 10, 20, 100)),
  #                     guide="colorbar") +
  #ggtitle("Location of peptide sequences in antibody") +
  xlab("Peptide") + ylab("Patient-ID") +
  guides(fill=guide_colourbar(barwidth = 0.5, barheight = 10, nrow = 1, title="% presence", title.position = "top")) +
  theme(legend.title.align=0.5)  +
  theme(legend.position = "right") 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+

PopPlot 

# Combine LocPolt and PopPlot
# build the plots 
PopPlot <- ggplot_gtable(ggplot_build(PopPlot))
HistPlot <- ggplot_gtable(ggplot_build(HistPlot))
# copy the plot width from plot1 to plot2
PopPlot$widths <- HistPlot$widths

# Draw and Save plot to pdf
pdf(file = paste(getwd(), "/CombinePlot_Peptides_by_PatientID_", Sys.Date(), ".pdf", sep = ""), width = 16, height = 10)
CombinePlot <- grid.arrange(HistPlot, PopPlot, heights=c(0.5, 0.5), ncol = 1)
dev.off()









#-------------------------------------------------------------------------------
# Extract peptides that appear in all BioProject
c <- Reduce(intersect, split(OAS2$Sequence, OAS2$BioProject))

OAS4 <- OAS2[OAS2$Sequence %in% c,]

# number of peptides shared between BioProject
nrow(OAS4)
# number of different peptides shared between BioProject
length(c)
# number of antibodies that above peptides identify
length(unique(OAS4$AntibodyID))

# Find frequency of each peptide in every BioProject
Freq_BioProject <- dcast(OAS4, BioProject~Sequence, length)
rownames(Freq_BioProject) <- Freq_BioProject$BioProject
Sequence <- colnames(Freq_BioProject)[2:ncol(Freq_BioProject)]

Freq_BioProject <- data.table::transpose(Freq_BioProject)

colnames(Freq_BioProject) <- as.character(Freq_BioProject[1,])
Freq_BioProject <- Freq_BioProject[-c(1),]
Freq_BioProject[, names(Freq_BioProject) := lapply(.SD, as.numeric)]

Freq_BioProject <- cbind(Sequence, Freq_BioProject)


OAS_shared_BioProject <- as.data.table(table(OAS4$Sequence)); colnames(OAS_shared_BioProject) <- c("Sequence", "Freq")
OAS_shared_BioProject <- OAS_shared_BioProject[order(-OAS_shared_BioProject$Freq),]
OAS_shared_BioProject <- left_join(x = OAS_shared_BioProject, y = Freq_BioProject, by = "Sequence")

# Convert frequency to percentage table
# Get total number of antibodies
Total.Abs <- length(unique(OAS3$AntibodyID))
# Get number of antibodies in each BioProject
BioProject <- colnames(OAS_shared_BioProject)[-c(1:2)]
No.Abs.each.BioProject <- sapply(c(1:length(BioProject)), function(i){
  length(unique(OAS3$AntibodyID[OAS3$BioProject == BioProject[i]]))
})

No.Abs <- c(Total.Abs, No.Abs.each.BioProject)

Freq <- as.matrix(OAS_shared_BioProject[,-c(1)])
Percentage <- round(100*sweep(Freq, 2, No.Abs, `/`), digits = 2)

OAS_shared_BioProject_Percent <- cbind(OAS_shared_BioProject[,c(1)], Percentage)


# save Frequency data by BioProject to excel file
library(openxlsx)
ExcelFile <- createWorkbook("Freq_Percent")
addWorksheet(ExcelFile, "Frequency")
writeData(ExcelFile, sheet = 1, OAS_shared_BioProject)
addWorksheet(ExcelFile, "Percent")
writeData(ExcelFile, sheet = 2, OAS_shared_BioProject_Percent)
saveWorkbook(ExcelFile, paste(getwd(), "/Freq_Table_by_BioProject_", Sys.Date(), ".xlsx", sep = ""), overwrite = TRUE)


#-------------------------------------------------------------------------------
# Draw histogram showing percentage appearance of peptides in antibodies

Hist_data <- data.table(SeqNo =  c(1:100),
                        Percent = OAS_shared_BioProject_Percent$Freq[1:100],
                        SeqLength = nchar(OAS_shared_BioProject_Percent$Sequence[1:100]))

HistPlot <- ggplot(data = Hist_data, aes(x=SeqNo, y=Percent, fill = factor(SeqLength), width=.75)) +
  geom_bar(stat="identity") +
  xlab("Peptide") + ylab("Appear in how many percent of Abs") +
  guides(fill=guide_legend(ncol = 2, title="Peptide \n length")) +
  theme(legend.position = "right") 
HistPlot


# Draw heat map of the OAS_shared_BioProject for seeing peptide sequences percentage in each BioProject
PopTable_plot <- OAS_shared_BioProject_Percent[c(1:100),-c(1:2)]
colnames(PopTable_plot) <- as.character(c(1:ncol(PopTable_plot)))

PopTable_plot <- PopTable_plot %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
PopTable_plot$rowname <- as.numeric(PopTable_plot$rowname)

PopPlot <- ggplot(PopTable_plot, aes(x = rowname, y = colname, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="grey", high="red") +
  #scale_fill_gradientn(colours=c("grey", "red", "red", "red", "red"),
  #                                values=rescale(c(0, 0.1, 10, 20, 100)),
  #                     guide="colorbar") +
  #ggtitle("Location of peptide sequences in antibody") +
  xlab("Peptide") + ylab("BioProject") +
  guides(fill=guide_colourbar(barwidth = 0.5, barheight = 10, nrow = 1, title="% presence", title.position = "top")) +
  theme(legend.title.align=0.5)  +
  theme(legend.position = "right") 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+

PopPlot 

# Combine LocPolt and PopPlot
# build the plots 
PopPlot <- ggplot_gtable(ggplot_build(PopPlot))
HistPlot <- ggplot_gtable(ggplot_build(HistPlot))
# copy the plot width from plot1 to plot2
PopPlot$widths <- HistPlot$widths

# Draw and Save plot to pdf
pdf(file = paste(getwd(), "/CombinePlot_Peptides_by_BioProject_", Sys.Date(), ".pdf", sep = ""), width = 16, height = 10)
CombinePlot <- grid.arrange(HistPlot, PopPlot, heights=c(0.5, 0.5), ncol = 1)
dev.off()



















#-------------------------------------------------------------------------------
# Find location of each peptide in antibody sequence (high computational cost)
# Get OAS data of full sequences
OAS_full <- OASdata %>% select(    sequence_alignment_aa,
                                   fwr1_aa,
                                   cdr1_aa,
                                   fwr2_aa,
                                   cdr2_aa,
                                   fwr3_aa,
                                   cdr3_aa,
                                   fwr4_aa) %>% 
  collect


# Find location of each peptide sequence in every antibody sequence


AbLoc <- lapply(c(1:nrow(OAS_shared_seq)), function(i){
  # Get peptide
  PepSeq <- as.character(OAS_shared_seq$Sequence[i])
  # count how many times this peptide appear in fwr1_aa, cdr1_aa, fwr2_aa, cdr2_aa, fwr3_aa, cdr3_aa, fwr4_aa
  data.table("Freq_Seq" = sum(str_count(OAS_full$sequence_alignment_aa, PepSeq), na.rm = TRUE),
             "fwr1_aa" = sum(str_count(OAS_full$fwr1_aa, PepSeq), na.rm = TRUE),
             "cdr1_aa" = sum(str_count(OAS_full$cdr1_aa, PepSeq), na.rm = TRUE),
             "fwr2_aa" = sum(str_count(OAS_full$fwr2_aa, PepSeq), na.rm = TRUE),
             "cdr2_aa" = sum(str_count(OAS_full$cdr2_aa, PepSeq), na.rm = TRUE),
             "fwr3_aa" = sum(str_count(OAS_full$fwr3_aa, PepSeq), na.rm = TRUE),
             "cdr3_aa" = sum(str_count(OAS_full$cdr3_aa, PepSeq), na.rm = TRUE),
             "fwr4_aa" = sum(str_count(OAS_full$fwr4_aa, PepSeq), na.rm = TRUE))
})
LocTable <- rbindlist(AbLoc)

OAS_shared_seq <- cbind(OAS_shared_seq, LocTable) 



#-------------------------------------------------------------------------------\
# Find frequency of each peptide in every DataID
Freq_DataID <- dcast(OAS_digested, DataID~Sequence, length)
rownames(Freq_DataID) <- Freq_DataID$DataID
Sequence <- colnames(Freq_DataID)[2:ncol(Freq_DataID)]

Freq_DataID <- data.table::transpose(Freq_DataID)

colnames(Freq_DataID) <- as.character(Freq_DataID[1,])
Freq_DataID <- Freq_DataID[-c(1),]
Freq_DataID[, names(Freq_DataID) := lapply(.SD, as.numeric)]

Freq_DataID <- cbind(Sequence, Freq_DataID)


OAS_shared_seq <- left_join(x = OAS_shared_seq, y = Freq_DataID, by = "Sequence")

# save Frequency data to csv file
write.csv(OAS_shared_seq, file = paste(getwd(), "/Frequency_Table_", Sys.Date(), ".csv", sep = ""))
#OAS_shared_seq <- fread("Frequency_Table_2023-02-14.csv")
#OAS_shared_seq <- OAS_shared_seq[,-c(1)]

#-------------------------------------------------------------------------------
# Convert Frequency to Percentage table

Percent_DataID <-  as.data.table(100*mapply('/', OAS_shared_seq[,-c(1:10)], colSums(OAS_shared_seq[,-c(1:10)])))

Percent_Loc <- OAS_shared_seq[,c(3:10)] / OAS_shared_seq$Freq_Seq *100

Percent_Ab <- OAS_shared_seq[,c(2)] / nrow(OASdata) *100

OAS_shared_seq_percent <- cbind(OAS_shared_seq[,c(1)], Percent_Ab, Percent_Loc, Percent_DataID)

# save Percentage data to csv file
write.csv(OAS_shared_seq_percent, file = paste(getwd(), "/Percent_Table_", Sys.Date(), ".csv", sep = ""))

# Save bio-project names
Bio_projects <- data.table("DataID"      = c(1:ncol(Percent_DataID)),
                           "Bio_project" = colnames(Percent_DataID))
write.csv(Bio_projects, file = paste(getwd(), "/DataID_Table_", Sys.Date(), ".csv", sep = ""), row.names = F)


#-------------------------------------------------------------------------------\
# Draw histogram showing percentage appearance of peptides in antibodies

Hist_data <- data.table(SeqNo =  c(1:nrow(OAS_shared_seq_percent)),
                        Percent = OAS_shared_seq_percent$Freq_digest,
                        SeqLength = nchar(OAS_shared_seq_percent$Sequence))

HistPlot <- ggplot(data = Hist_data, aes(x=SeqNo, y=Percent, fill = factor(SeqLength), width=.75)) +
  geom_bar(stat="identity") +
  xlab("Peptide") + ylab("Appear in how many percent of Abs") +
  guides(fill=guide_legend(ncol = 2, title="Peptide \n length")) +
  theme(legend.position = "right") 
HistPlot

#-------------------------------------------------------------------------------
# Draw heat map of the OAS_shared_seq for seeing peptide sequences percentage in each sampleID
PopTable_plot <- Percent_DataID
colnames(PopTable_plot) <- as.character(c(1:ncol(Percent_DataID)))

PopTable_plot <- PopTable_plot %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
PopTable_plot$rowname <- as.numeric(PopTable_plot$rowname)

PopPlot <- ggplot(PopTable_plot, aes(x = rowname, y = as.numeric(colname), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="grey", high="red") +
  #scale_fill_gradientn(colours=c("grey", "red", "red", "red", "red"),
  #                                values=rescale(c(0, 0.1, 10, 20, 100)),
  #                     guide="colorbar") +
  #ggtitle("Location of peptide sequences in antibody") +
  xlab("Peptide") + ylab("Bio-project") +
  guides(fill=guide_colourbar(barwidth = 0.5, barheight = 10, nrow = 1, title="% presence", title.position = "top")) +
  theme(legend.title.align=0.5)  +
  theme(legend.position = "right") 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+

PopPlot  

#-------------------------------------------------------------------------------
# Draw heat map of the OAS_shared_seq for seeing peptide sequences location
LocTable_plot <- Percent_Loc[,-c(1:2)] %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
LocTable_plot$rowname <- as.numeric(LocTable_plot$rowname)

LocationIndicator <- c("fwr1_aa",
                       "cdr1_aa",
                       "fwr2_aa",
                       "cdr2_aa",
                       "fwr3_aa",
                       "cdr3_aa",
                       "fwr4_aa")

LocPlot <- ggplot(LocTable_plot, aes(x = rowname, y = factor(colname, levels = LocationIndicator), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="grey", high="blue") +
  #scale_fill_gradientn(colours=c("grey", "red", "red", "red", "red"),
  #                                values=rescale(c(0, 0.1, 10, 20, 100)),
  #                     guide="colorbar") +
  #ggtitle("Location of peptide sequences in antibody") +
  xlab("Peptide") + ylab("Location") +
  guides(fill=guide_colourbar(barwidth = 0.5, barheight = 5, nrow = 1, title="% presence", title.position = "top")) +
  theme(legend.title.align=0.5)  +
  theme(legend.position = "right") 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+

LocPlot  

#-------------------------------------------------------------------------------
# Combine LocPolt and PopPlot

# build the plots 
PopPlot <- ggplot_gtable(ggplot_build(PopPlot))
LocPlot <- ggplot_gtable(ggplot_build(LocPlot))
HistPlot <- ggplot_gtable(ggplot_build(HistPlot))
# copy the plot width from plot1 to plot2
PopPlot$widths <- LocPlot$widths
HistPlot$widths <- LocPlot$widths
# Draw the combine plot
CombinePlot <- grid.arrange(HistPlot, PopPlot, LocPlot, heights=c(0.4, 0.4, 0.2), ncol = 1)

# Save plot to pdf
pdf(file = paste(getwd(), "/CombinePlot_Common_Peptides_", Sys.Date(), ".pdf", sep = ""), width = 16, height = 10)
# Draw the combine plot
CombinePlot <- grid.arrange(HistPlot, PopPlot, LocPlot, heights=c(0.4, 0.4, 0.2), ncol = 1)
dev.off()
