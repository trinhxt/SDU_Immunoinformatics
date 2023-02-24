library(svDialogs)
DataFolder <- dlg_dir(title = "Select csv data folder:", filters = dlg_filters[c("All"), ])$res
setwd(DataFolder)
DataFiles <- list.files(path = DataFolder, pattern = ".csv.gz", all.files = FALSE, 
                        full.names = TRUE, recursive = TRUE, ignore.case = FALSE, 
                        include.dirs = FALSE, no.. = FALSE)
DataSummary <- read.csv(file = "OAS-SARS-COV-2-summary.csv", header = T)
library(data.table)
Listdata <- lapply(c(1:length(DataFiles)), function(i){
  dat <- fread((DataFiles[i]), header = TRUE, skip = 1)
  Datafilename <- substr(basename(DataFiles[i]), 1, nchar(basename(DataFiles[i]))-7)
  cbind(Datafilename, dat)
})
Datatable <- rbindlist( Listdata )
Datatable$AntibodyID <- c(1:nrow(Datatable))
# Draw histogram of antibody sequences in each sample
library(ggplot2)
PlotData <- as.data.frame(table(Datatable$Datafilename))
PlotData$DataID <- strtrim(PlotData$Var1, nchar(as.character(PlotData$Var1))-13)
PlotData$DataID <- factor(PlotData$DataID, levels = stringr::str_sort(unique(PlotData$DataID), numeric = TRUE))
HistPlot <- ggplot(data = PlotData, aes(x=DataID, y=Freq)) +
  geom_bar(position="stack", stat="identity", fill="orange", alpha=0.8) +
  geom_text(aes(label = Freq), vjust = -0.4) +
  xlab("Sample ID") + ylab("Number of antibody sequences") +
  theme(axis.text.x=element_text(angle=90, hjust=0, vjust = 0.5))
HistPlot
# Save plot to pdf
pdf(file = paste(DataFolder, "/Histogram_Samples_Antibodies_", Sys.Date(), ".pdf", sep = ""), width = 10, height = 5)
HistPlot
dev.off()

Sequences <- as.list(c(Datatable$sequence_alignment_aa_heavy, Datatable$sequence_alignment_aa_light))

NameList <- strtrim(Datatable$Datafilename, nchar(as.character(Datatable$Datafilename))-13)
Nameseq <- c(paste(NameList,Datatable$AntibodyID, "heavy ", sep = "_"), paste(NameList,Datatable$AntibodyID,"light ", sep = "_"))
                
# Save data to fasta file
library(seqinr)
write.fasta(Sequences,   # sequences
            Nameseq, # name of sequences
            paste(DataFolder,"/OAS_SARS_COV2_", Sys.Date(), ".fasta", sep = ""), # location to save fasta file  
            open = "w", nbchar = 60, as.string = TRUE)

# Save data to csv file
#write.csv(Datatable, file = paste(DataFolder, "/OAS_SARS_COV2_", Sys.Date(), ".csv", sep = ""), row.names = F)
write.table(Datatable, file = paste(DataFolder, "/OAS_SARS_COV2_", Sys.Date(), ".txt", sep = ""), row.names = F)
# Save data to .RData file
# save(Datatable, file = "OAS-SARS-COV-2.RData")