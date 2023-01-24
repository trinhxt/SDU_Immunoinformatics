library(svDialogs)
# Set working directory, which is folder containing digested data
setwd(dlg_dir(title = "Select working folder", filters = dlg_filters[c("All"), ])$res)

# Read data OAS and UNIPROT
OASdata <- read.table("OAS_SARS_COV2_digested_Mass400to6000.txt", header = T)
UNIPROTdata <- read.table("Uniprot-human-reviewed-20230103_digested_Mass400to6000.txt", header = T)

# Find common sequences in OAS and UNIPROT
a <- Reduce(intersect, list(OASdata$Sequence, UNIPROTdata$Sequence))

# Find location of those common sequences in OAS data
b <- sapply(c(1:length(a)), function(i){
  which(OASdata$Sequence == a[i])
})

# Mark those common sequences as number 1 in column UNIPROTcheck of OAS data
OASdata$UNIPROTcheck <- 0
OASdata$UNIPROTcheck[unlist(b)] <- 1

# Find location of those common sequences in UNIPROT data
c <- sapply(c(1:length(a)), function(i){
  which(UNIPROTdata$Sequence == a[i])
})
d <- lapply(c, `[[`, 1)

# Assign Protein name from UNIPROT data to OAS data
OASdata$UNIPROTname <- NA; 
for (i in c(1:length(d))) {
  e1 <- which(OASdata$Sequence == UNIPROTdata$Sequence[d[[i]]])
  OASdata$UNIPROTname[e1] <- UNIPROTdata$Protein_Name[d[[i]]]
}

# Make result table and save to csv file:
ResultTable <- data.frame("Sequence"      = OASdata$Sequence,
                          "OAS_name"      = OASdata$Protein_Name,
                          "UNIPROT_name"  = OASdata$UNIPROTname,
                          "OAS_check"     = 1,
                          "UNIPROT_check" = OASdata$UNIPROTcheck)

write.csv(ResultTable, file = paste(getwd(), "/OAS_UNIPROT_checked.csv", sep = ""))
