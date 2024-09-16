################################################################################
## Section 0.0: Set a custom library path (optional)
custom_lib <- "/work/Immunoinformatics/R-packages-RStudio/" # Change this to your directory
#custom_lib <- "/work/Immunoinformatics/R-packages-Ubuntu/"
if (!dir.exists(custom_lib)) {
  dir.create(custom_lib, recursive = TRUE)
}
.libPaths(custom_lib)


################################################################################
## Section 0.1: List of CRAN and Bioconductor packages
cran_packages <- c("arrow", "dplyr", "jsonlite", "seqinr", "data.table", "fs", "stringr", "stringi", "duckdb", "DBI", "RSQLite", "plotly")
bioc_packages <- c("cleaver", "Biostrings", "GenomeInfoDb")

# Function to install missing CRAN packages
install_missing_cran <- function(packages, lib) {
  installed <- rownames(installed.packages(lib.loc = lib))
  for (pkg in packages) {
    if (!pkg %in% installed) {
      install.packages(pkg, lib = lib, dependencies = TRUE)
    }
  }
}

################################################################################
# Section 0.2: Install missing CRAN packages to the custom library
install_missing_cran(cran_packages, custom_lib)
# Function to install missing Bioconductor packages
install_missing_bioc <- function(packages, lib) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = lib)
  }
  installed <- rownames(installed.packages(lib.loc = lib))
  for (pkg in packages) {
    if (!pkg %in% installed) {
      BiocManager::install(pkg, lib = lib, update = TRUE, ask = FALSE)
    }
  }
}


################################################################################
# Section 0.3: Install missing Bioconductor packages to the custom library
install_missing_bioc(bioc_packages, custom_lib)



################################################################################
# Section 0.4: Load installed packages from the custom library
# Ensure the custom library path is set again if starting a new session
.libPaths(custom_lib)
# Check the library paths
print(.libPaths())
# Load all packages
lapply(c(cran_packages, bioc_packages), require, character.only = TRUE)
