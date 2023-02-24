library(svDialogs)
# Set working directory, which is folder containing digested data
setwd(dlg_dir(title = "Select working folder", filters = dlg_filters[c("All"), ])$res)

# Read data OAS and UNIPROT
OASdata <- read.table("OAS_SARS_COV2_digested.txt", header = T)
UNIPROTdata <- read.table("UNIPROT_human_digested.txt", header = T)
#load("OAS-SARS-COV-2.RData")
Datatable <- read.csv(file = "OAS_SARS_COV2.csv", header = T)


#-------------------------------------------------------------------------------
# Find common sequences in OAS and UNIPROT
a <- Reduce(intersect, list(OASdata$Sequence, UNIPROTdata$Sequence))

# Find location of those common sequences in OAS data
b <- sapply(c(1:length(a)), function(i){
  which(OASdata$Sequence == a[i])
})

# Mark those common sequences as number 1 in column UNIPROTcheck of OAS data
OASdata$UNIPROTcheck <- 0
OASdata$UNIPROTcheck[unlist(b)] <- 1
OASdata$DataID <- substr(OASdata$Protein_Name, 1, 11)

# Draw histogram showing number of digested sequences in each sample
library(ggplot2)
PlotData <- as.data.frame(table(OASdata$DataID)); colnames(PlotData) <- c("DataID", "Freq")
PlotData$DataID <- factor(PlotData$DataID, levels = stringr::str_sort(unique(PlotData$DataID), numeric = TRUE))

# Create data table for plot
DataIDList <- unique(OASdata$DataID)
FreqTableDataIDList <- lapply(c(1:length(DataIDList)), function(i){
  a <- subset(OASdata, OASdata$DataID == DataIDList[i])
  b <- as.data.frame(table(a$UNIPROTcheck))
  b$DataID <- DataIDList[i]
  b
})
# Combine the above list of data tables into one data table
library(dplyr)
FreqTableDataID <- bind_rows(FreqTableDataIDList, .id = "column_label")
colnames(FreqTableDataID) <- c("No.", "UNIPROTcheck", "Frequency", "DataID")
FreqTableDataID$Percentage <- FreqTableDataID$Frequency*100/sum(FreqTableDataID$Frequency)

# Draw histogram of sequence frequency
library(ggplot2)
Hist_data <- FreqTableDataID
HistPlot <- ggplot(data = Hist_data, aes(x=DataID, y=Frequency, fill=factor(UNIPROTcheck, levels = c("1", "0")), alpha=as.factor(UNIPROTcheck))) +
  geom_bar(position="stack", stat="identity") +
  geom_text(aes(label = Frequency), position = position_stack(vjust = 0.5)) +
  xlab("Sample ID") + ylab("Number of digested peptide sequences") +
  guides(fill=guide_legend(nrow = 1, title="Appear in UNIPROT?")) +
  theme(axis.text.x=element_text(angle=90,hjust=1, vjust = 0.7)) +
  theme(legend.position = "top") +
  scale_alpha_manual(values = c(1, 0.5))
HistPlot

# Save plot to pdf
pdf(file = paste(getwd(), "/Histogram_Peptides_UNIPROTcheck_", Sys.Date(), ".pdf", sep = ""), width = 10, height = 5)
HistPlot
dev.off()

# Subset OAS data that is not in UNIPROT
OASchecked <- subset(OASdata, OASdata$UNIPROTcheck == 0)
library(stringr)
# Add AntibodyID to OAS data
AntibodyID <- str_match(OASchecked$Protein_Name, "_\\s*(.*?)\\s*_")
OASchecked$AntibodyID <- as.numeric(AntibodyID[,2])

# Check how many unique peptide sequences
length(unique(OASchecked$Protein_Name))

#-------------------------------------------------------------------------------
# Find starting point in antibody-sequence that a PepSeq belong to
library(stringr)
StartPosition <- sapply(c(1:nrow(OASchecked)), function(i){
  # Get peptide
  PepSeq <- as.character(OASchecked$Sequence[i])
  # Get antibody IDs of the above peptide
  AbID <- OASchecked$AntibodyID[i]
  # Get a small data table containing above antibody ID and peptide
  DatatablePep <- Datatable[c(AbID), c("sequence_alignment_aa_heavy", "sequence_alignment_aa_light")]
  # Find the starting point of the peptide sequence
  str_locate(DatatablePep, PepSeq)[,1]
})
StartPosition <- t(StartPosition)

# Convert the above list to a data table
StartTable <- as.data.frame((StartPosition))
library(tidyr)
StartTable$V1 <- StartTable$V1 %>% replace_na(0)
StartTable$V2 <- StartTable$V2 %>% replace_na(0)

# Calculate starting position in heavy/light chain and add to OASchecked
OASchecked$StartPosition <- colSums(t(StartTable))

#-------------------------------------------------------------------------------
# Find location of each peptide sequence in the antibody sequence
LocationIndicator <- c("fwr1_aa_heavy", "cdr1_aa_heavy", "fwr2_aa_heavy", "cdr2_aa_heavy", "fwr3_aa_heavy", "cdr3_aa_heavy", "fwr4_aa_heavy",
                        "fwr1_aa_light", "cdr1_aa_light", "fwr2_aa_light", "cdr2_aa_light", "fwr3_aa_light", "cdr3_aa_light", "fwr4_aa_light")
AbLoc <- sapply(c(1:nrow(OASchecked)), function(i){
  # Get peptide
  PepSeq <- as.character(OASchecked$Sequence[i])
  # Get antibody ID of the above peptide
  AbID <- OASchecked$AntibodyID[i]
  # Get a small data table containing above antibody ID and peptide
  DatatablePep <- Datatable[c(AbID), c(LocationIndicator)]
  # Find which column (location in the whole peptide sequence) that a PepSeq belong to
  str_detect(DatatablePep, PepSeq)
})
LocTable <- as.data.frame(t(AbLoc)); colnames(LocTable) <- LocationIndicator
LocTable[] <- lapply(LocTable, function(x) ifelse(x == TRUE, 1, 0))

# Add a column in OASchecked table to indicate position of the pep sequence
OASchecked$SequenceLocation <- as.character(lapply(apply(LocTable == 1, 1, which), names))
# Change character "character(0)" of column SequenceLocation to "Other"
OASchecked$SequenceLocation[OASchecked$SequenceLocation == "character(0)"] <- "Other"
# Get row index of items that have long name such as c("cdr1_aa_heavy", "fwr3_aa_heavy")
seq2posloc <- which(nchar(OASchecked$SequenceLocation) > 13)
seq2pos <- sapply(seq2posloc, function(i){
              string <- OASchecked$SequenceLocation[i]
              b <- str_locate(string, "cdr")
              str_sub(string, start = b[1,1], end = b[1,1]+12)
})
OASchecked$SequenceLocation[seq2posloc] <- seq2pos
# Change character NA of column SequenceLocation to "Other"
library(tidyr)
OASchecked$SequenceLocation <- OASchecked$SequenceLocation %>% replace_na("Other")

# Save OAS data to csv file 
write.csv(OASchecked[,-c(1)], file = paste(getwd(), "/OAS_SARS_COV2_checked_", Sys.Date(), ".csv", sep = ""), row.names = F)


#-------------------------------------------------------------------------------
# Draw stacked histogram: x-asis is starting position and y-axis is frequency of peptide sequence, 
# Create data table for plot
StartPositionList <- unique(OASchecked$StartPosition); StartPositionList <- sort(StartPositionList, decreasing = FALSE)
FreqTableList <- lapply(c(1:length(StartPositionList)), function(i){
  a <- subset(OASchecked, OASchecked$StartPosition == StartPositionList[i])
  b <- as.data.frame(table(a$SequenceLocation))
  b$Start_Position <- StartPositionList[i]
  b
})
# Combine the above list of data tables into one data table
library(dplyr)
FreqTable <- bind_rows(FreqTableList, .id = "column_label")
colnames(FreqTable) <- c("No.", "Location", "Frequency", "Start_Position")
FreqTable$Percentage <- FreqTable$Frequency*100/sum(FreqTable$Frequency)

# Draw histogram of sequence frequency
library(ggplot2)
Hist_data <- FreqTable
#Hist_data$Location <- factor(Hist_data$Location, levels = c(LocationIndicator, "Other"))
HistPlot <- ggplot(data = Hist_data, aes(x=Start_Position, y=Percentage, fill = Location)) +
  geom_bar(position="stack", stat="identity", colour="black", alpha=1.0, linewidth=0.02) +
  xlab("Peptide starting position") + ylab("Percentage (%)") +
  guides(fill=guide_legend(ncol = 1, title="Peptide location")) +
  theme(axis.text.x=element_text(angle=0,vjust=0.5)) +
  theme(legend.position = "right")
HistPlot

# Save plot to pdf
pdf(file = paste(getwd(), "/Histogram_Peptides_Start_Position_", Sys.Date(), ".pdf", sep = ""), width = 12, height = 6)
HistPlot
dev.off()



#-------------------------------------------------------------------------------
# Find shared peptide sequences in all samples
SampleID <- unique(OASdata$DataID)
SampleIDList <- lapply(c(1:length(SampleID)), function(i){
  subset(OASchecked$Sequence, OASchecked$DataID == SampleID[i])
})
a <- Reduce(intersect, SampleIDList)

# Find location of those common peptides in OASchecked
b <- sapply(c(1:length(a)), function(i){
  which(OASchecked$Sequence == a[i])
})

# Subset data table of those common peptide
OAScheckedCommon <- OASchecked[unlist(b),]
PepSeqCommon <- as.data.frame(table(OAScheckedCommon$Sequence)); colnames(PepSeqCommon) <- c("Sequence", "Freq")
PepSeqCommon$SeqLength <- nchar(as.character(PepSeqCommon$Sequence))
PepSeqCommon <- PepSeqCommon[order(-PepSeqCommon$Freq),]
PepSeqCommon$SeqNo <- c(1:nrow(PepSeqCommon))

# Find location of those peptides in antibody sequences
AbLoc <- sapply(c(1:nrow(PepSeqCommon)), function(i){
  # Get peptide
  PepSeq <- as.character(PepSeqCommon$Sequence[i])
  # Find row index the above peptide in OAScheckedCommon data table
  loc <- which(OAScheckedCommon$Sequence == PepSeq)
  # Get antibody ID of the above peptide
  AbID <- OAScheckedCommon$AntibodyID[loc]
  # Get a small data table containing above antibody ID and peptide
  DatatablePep <- Datatable[c(AbID), c(LocationIndicator)]
  # Find which column (location in the whole peptide sequence) that a PepSeq belong to
  str_detect(DatatablePep, PepSeq)
})
LocTable <- as.data.frame(t(AbLoc)); colnames(LocTable) <- LocationIndicator
LocTable[] <- lapply(LocTable, function(x) ifelse(x == TRUE, 1, 0))

a <- unique(OAScheckedCommon[,c(2,3,4,11)])
# Save list of common peptides and their location
write.csv(cbind(PepSeqCommon, LocTable), file = paste(getwd(), "/OAS_SARS_COV2_Common_Peptides_", Sys.Date(), ".csv", sep = ""), row.names = F)

#-------------------------------------------------------------------------------
# Draw heat map of the OAScheckedCommon for seeing peptide sequences frequency and location
library(tidyverse)
library(ggplot2)
library(gtools)

LocTable_plot <- LocTable %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
LocTable_plot$rowname <- as.numeric(LocTable_plot$rowname)


LocPlot <- ggplot(LocTable_plot, aes(x = rowname, y = factor(colname, levels = LocationIndicator), fill = value)) +
  geom_tile() +
  scale_fill_gradient(low="grey", high="red") +
  #ggtitle("Location of peptide sequences in antibody") +
  xlab("Peptide sequence ID") + ylab("Location in antibody sequence") +
  guides(fill=guide_legend(nrow = 1, title="Location\n (1=Yes, 0=No)")) +
  theme(legend.title.align=0.5)  +
  theme(legend.position = "bottom") 
#theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) #+
#scale_x_continuous(breaks = seq(0, nrow(AbSeqTable2), by = 10))
LocPlot  

# Draw histogram of sequence frequency
Hist_data <- PepSeqCommon
HistPlot <- ggplot(data = Hist_data, aes(x=SeqNo, y=Freq, fill = factor(SeqLength))) +
  geom_bar(stat="identity") +
  xlab("Peptide sequence ID") + ylab("Frequency") +
  guides(fill=guide_legend(nrow = 2, title="Sequence length")) +
  theme(legend.position = "top") 
HistPlot

library(gridExtra)
# build the plots 
HistPlot <- ggplot_gtable(ggplot_build(HistPlot))
LocPlot <- ggplot_gtable(ggplot_build(LocPlot))

# copy the plot width from plot1 to plot2
HistPlot$widths <- LocPlot$widths

# Draw the combine plot
CombinePlot <- grid.arrange(HistPlot, LocPlot, ncol = 1, nrow = 2)

# Save plot to pdf
pdf(file = paste(getwd(), "/CombinePlot_Common_Peptides_", Sys.Date(), ".pdf", sep = ""), width = 12, height = 6)
# Draw the combine plot
CombinePlot <- grid.arrange(HistPlot, LocPlot, ncol = 1, nrow = 2)
dev.off()
