library(svDialogs)
DataFolder <- dlg_dir(title = "Select fasta data folder:", filters = dlg_filters[c("All"), ])$res
DataFiles <- list.files(path = DataFolder, pattern = ".csv.gz", all.files = FALSE, 
                        full.names = TRUE, recursive = TRUE, ignore.case = FALSE, 
                        include.dirs = FALSE, no.. = FALSE)

library(data.table)
Listdata <- lapply(c(1:length(DataFiles)), function(i){
  dat <- fread((DataFiles[i]), header = TRUE, skip = 1)
  Datafilename <- substr(basename(DataFiles[i]), 1, nchar(basename(DataFiles[i]))-7)
  cbind(Datafilename, dat)
})
Datatable <- rbindlist( Listdata )


Sequences <- as.list(c(Datatable$sequence_alignment_aa_heavy, Datatable$sequence_alignment_aa_light))

Nameseq <- c(paste(Datatable$Datafilename,"_heavy ", sep = ""), paste(Datatable$Datafilename,"_light ", sep = ""))
                

library(seqinr)
write.fasta(Sequences,   # sequences
            Nameseq, # name of sequences
            paste(DataFolder,"/OAS_SARS_COV2.fasta", sep = ""), # location to save fasta file  
            open = "w", nbchar = 60, as.string = TRUE)
