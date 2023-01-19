library(svDialogs)
DataFolder <- dlg_dir(title = "Select fasta data folder:", filters = dlg_filters[c("All"), ])$res
DataFiles <- list.files(path = DataFolder, pattern = ".csv.gz", all.files = FALSE, 
                        full.names = TRUE, recursive = TRUE, ignore.case = FALSE, 
                        include.dirs = FALSE, no.. = FALSE)

library(data.table)
Listdata <- lapply(DataFiles, fread, sep=",")

Listdata <- lapply(c(1:length(DataFiles)), function(i){
  dat <- fread((DataFiles[i]), header = TRUE, skip = 1)
  Datafilename <- substr(basename(DataFiles[i]), 1, nchar(basename(DataFiles[i]))-7)
  cbind(Datafilename, dat)
})

Datatable <- rbindlist( Listdata )

Sequences <- as.list(Datatable$sequence_alignment_aa_heavy)

Nameseq <- paste(Datatable$Datafilename, Datatable$sequence_id_heavy, Datatable$sequence_id_light, sep = " ")

#OAS_dat <- fread(DataFiles[1], header = TRUE, skip = 1)

library(seqinr)
write.fasta(Sequences, 
            Nameseq, # name of sequences
            paste(DataFolder,"/OAS_SARS_COV2.fasta", sep = ""), # location to save fasta file  
            open = "w", nbchar = 60, as.string = TRUE)
