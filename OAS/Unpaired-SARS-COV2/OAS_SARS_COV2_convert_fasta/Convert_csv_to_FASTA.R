library(svDialogs)
DataFolder <- dlg_dir(title = "Select csv data folder:", filters = dlg_filters[c("All"), ])$res
setwd(DataFolder)
DataFiles <- list.files(path = DataFolder, pattern = ".csv.gz", all.files = FALSE, 
                        full.names = TRUE, recursive = TRUE, ignore.case = FALSE, 
                        include.dirs = FALSE, no.. = FALSE)
DataSummary <- read.csv(file = "OAS-SARS-COV-2-summary.csv", header = T)
DataSummary$Filename <- basename(DataSummary$DownloadLink)
DataSummary$PatientID <- paste(DataSummary$DS.Name, DataSummary$Individual, sep = "_")
#DataSummary <- subset(DataSummary, (DataSummary$Unique.Sequences > 1000))
nrow(unique(DataSummary[,c(2,9)]))

# Choose only data file that have number of unique sequences > 1000 and < 300000
a <- sapply(c(1:nrow(DataSummary)), function(i){
  which(basename(DataFiles) == DataSummary$Filename[i])
})






# Summary of data
# Draw histogram of antibody sequences in each SRA data

DataPlot <- data.frame("SampleID" = c(1:nrow(DataSummary)),
                       "No.AbSeq" = DataSummary$Unique.Sequences)
library(ggplot2)
HistPlot <- ggplot(data = DataPlot, aes(x=SampleID, y=No.AbSeq, fill = "orange")) +
  geom_bar(position="stack", stat="identity", colour="orange", alpha=1.0, linewidth=0.02) +
  xlab("SRA file") + ylab("Number of antibody sequences") +
  guides(fill=guide_legend(ncol = 1, title=" ")) +
  theme(axis.text.x=element_text(angle=0,vjust=0.5)) +
  theme(legend.position = "right")
HistPlot
# Save plot to pdf
pdf(file = paste(getwd(), "/Histogram_SRA_", Sys.Date(), ".pdf", sep = ""), width = 12, height = 6)
HistPlot
dev.off()



# Draw histogram of antibody sequences in each Bio-Project
library(dplyr)
DataPlot <- DataSummary %>% 
                group_by(DS.Name) %>% 
                summarise(No.AbSeq = sum(Unique.Sequences))
colnames(DataPlot) <- c("SampleID", "No.AbSeq")

library(ggplot2)
HistPlot <- ggplot(data = DataPlot, aes(x=SampleID, y=No.AbSeq, fill = "grey")) +
  geom_bar(position="stack", stat="identity", colour="grey", alpha=1.0, linewidth=0.02) +
  xlab("Bio-Project") + ylab("Number of antibody sequences") +
  guides(fill=guide_legend(ncol = 1, title=" ")) +
  theme(axis.text.x=element_text(angle=0,vjust=0.5)) +
  theme(legend.position = "right")
HistPlot
# Save plot to pdf
pdf(file = paste(getwd(), "/Histogram_BioProject_", Sys.Date(), ".pdf", sep = ""), width = 12, height = 6)
HistPlot
dev.off()



# Draw histogram of antibody sequences in each Patient
library(dplyr)
DataPlot <- DataSummary %>% 
  group_by(PatientID) %>% 
  summarise(No.AbSeq = sum(Unique.Sequences))
DataPlot$SampleID <- c(1:nrow(DataPlot))

library(ggplot2)
HistPlot <- ggplot(data = DataPlot, aes(x=SampleID, y=No.AbSeq, fill = "cyan")) +
  geom_bar(position="stack", stat="identity", colour="cyan", alpha=1.0, linewidth=0.02) +
  xlab("Patient") + ylab("Number of antibody sequences") +
  guides(fill=guide_legend(ncol = 1, title=" ")) +
  theme(axis.text.x=element_text(angle=0,vjust=0.5)) +
  theme(legend.position = "right")
HistPlot
# Save plot to pdf
pdf(file = paste(getwd(), "/Histogram_Patient_", Sys.Date(), ".pdf", sep = ""), width = 12, height = 6)
HistPlot
dev.off()




#-------------------------------------------------------------------------------
library(dplyr)
library(disk.frame)
library(data.table)
setup_disk.frame()

tempFolder <- "D:/DwnlData/Diskframe"
DataFilesFilter <- DataFiles[a]
nfiles <- length(DataFilesFilter)

Listdata <- lapply(c(1:nfiles), function(i){
  dat <- csv_to_disk.frame(
    DataFilesFilter[i], 
    outdir = file.path(paste(tempFolder, "/tmp_",i, sep = "")),
    nchunks = 60,
    overwrite = T)
})

# Load disk frames which were saved on hard disk
#Listdata <- lapply(c(1:nfiles), function(i){
  
#  disk.frame(paste("D:/DwnlData/Diskframe/tmp_", i, sep = ""))
  
#})


# row-bind two disk.frames
library(fst)
DatatableDisk <- rbindlist.disk.frame(Listdata,
                                  outdir = (paste(tempFolder, "/CombinedData", sep = "")),
                                  overwrite = T,
                                  #by_chunk_id = TRUE,
                                  parallel = TRUE)

Datatable  <- DatatableDisk %>% select(sequence_alignment_aa,
                                       fwr1_aa,
                                       cdr1_aa,
                                       fwr2_aa,
                                       cdr2_aa,
                                       fwr3_aa,
                                       cdr3_aa,
                                       fwr4_aa) %>%
                                collect

Datatable$AntibodyID <- c(1:nrow(DatatableDisk))

Datafilename <- lapply(c(1:nfiles), function(i){
  rep(basename(DataFilesFilter[i]), DataSummary$Unique.Sequences[i])
})

Datatable$Datafilename <- unlist(Datafilename)



#-------------------------------------------------------------------------------


Sequences <- as.list(c(Datatable$sequence_alignment_aa))

NameList <- strtrim(Datatable$Datafilename, nchar(as.character(Datatable$Datafilename))-7)

Nameseq <- c(paste(NameList, Datatable$AntibodyID, sep = "__"))
                
# Save data to fasta file
library(seqinr)
write.fasta(Sequences,   # sequences
            Nameseq, # name of sequences
            paste(dirname(DataFolder),"/OAS_SARS_COV2_", Sys.Date(), ".fasta", sep = ""), # location to save fasta file  
            open = "w", nbchar = 60, as.string = TRUE)

# Save data to feather file
library(feather)
write_feather(Datatable, paste(dirname(DataFolder),"/OAS_SARS_COV2_", Sys.Date(), ".feather", sep = ""))

# Load feather datafile
#Datatable <- read_feather(paste(dirname(DataFolder),"/OAS_SARS_COV2_", Sys.Date(), ".feather", sep = ""), columns = NULL)
