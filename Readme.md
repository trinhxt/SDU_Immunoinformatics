# OASpepDB: Human Disease CDR3 Antibody Peptides Database

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R: >= 4.2](https://img.shields.io/badge/R-%3E%3D%204.2-276DC3.svg)](https://www.r-project.org/)
[![Python: >= 3.9](https://img.shields.io/badge/Python-%3E%3D%203.9-3776AB.svg)](https://www.python.org/)
[![Database: DuckDB & Parquet](https://img.shields.io/badge/Database-DuckDB%20%7C%20Parquet-FFF000.svg)](https://duckdb.org/)

A curated repository of **65,510,795** non-redundant, **disease-exclusive antibody CDR3 peptides** across **25 human disease cohorts**, optimized for bottom-up immunoproteomics and liquid chromatography–tandem mass spectrometry (LC-MS/MS) database searching.

---

## 1. Scientific Overview & Rationale

Direct identification of circulating, disease-specific antibodies from patient biofluids (e.g., plasma, serum, cerebrospinal fluid) via bottom-up mass spectrometry requires sequence databases that capture hypervariable complementarity-determining regions (specifically CDR-H3). 

Standard reference proteomes such as UniProtKB/Swiss-Prot contain fewer than 40,000 immunoglobulin entries, omitting somatic hypermutations and V(D)J combinatorial diversity. Conversely, searching raw next-generation sequencing repertoires (>2 billion reads) causes severe database expansion, prolonged search times, and uncalibrated false discovery rates.

**OASpepDB** addresses this trade-off through a multi-tier negative subtraction and clonotype curation pipeline:

1. **Repertoire Aggregation**: Ingestion of 14,433 paired and unpaired human B-cell repertoire datasets from the [Observed Antibody Space (OAS)](https://opig.stats.ox.ac.uk/webapps/oas/).
2. **Tier-1 Negative Subtraction (Healthy Repertoire Background)**: Identification and subtraction of all CDR3 sequences observed across **8,089 healthy control repertoires** to eliminate common germline and non-disease background antibodies.
3. **Tier-2 & Tier-3 Negative Subtraction (Human Reference Proteomes)**: *In silico* tryptic digestion of UniProtKB/Swiss-Prot (canonical and isoforms) and NCBI RefSeq (GRCh38.p14) to purge any peptide fragments matching the human background proteome.
4. **Disease-Exclusive Curation**: Retention of sequences unique to each of the 25 clinical conditions, with clonal sharing quantified across distinct patients ($N \ge 1$) and read depths.
5. **Empirical FDR Control**: Integration of a calibrated ~1% non-human Camelid VHH entrapment library (NCBI GenBank) formatted with isobaric $I \rightarrow L$ substitution for target-decoy mass spectrometry validation.

---

## 2. Repository Structure

```text
SDU_Immunoinformatics/
├── OASpepDB/                               # Core application and data assets
│   ├── DBquery.R                           # Interactive Shiny web application (DuckDB OLAP)
│   ├── DBquery-Windows.bat                 # Launch script for Windows
│   ├── DBquery-Mac.command                 # Launch script for macOS
│   ├── DBquery-Linux.sh                    # Launch script for Linux
│   ├── Data/
│   │   ├── OAS_metadata.csv                # Curated metadata catalog (14,433 studies)
│   │   ├── cRAP.fasta                      # Common Repository of Adventitious Proteins
│   │   ├── entrapment_cassettes.fasta      # Camelid VHH entrapment sequences for empirical FDR
│   │   └── bulk_download_human_unpaired.sh # Bulk download script for OAS repertoires
│   └── Scripts/
│       ├── 01_download_OAS.py              # Step 01: Parallel download of OAS repertoires
│       ├── 02_build_healthy_cdr3_cache.py  # Step 02: Multi-threaded healthy CDR3 indexer
│       ├── 03_insilico_digest.R            # Step 03: Human proteome in silico digestion
│       ├── 04_build_disease_db.py          # Step 04: Disease-exclusive database builder
│       ├── 05_fetch_ncbi_entrapment.py     # Step 05: Camelid VHH entrapment library generator
│       └── generate_unified_alluvial_svg.R # Dynamic alluvial diagram generator
├── Archived/                               # Prior versions, deprecated files, and notes
├── .gitignore                              # Git exclusion rules for large data artifacts
├── LICENSE                                 # MIT License
└── Readme.md                               # Project documentation
```

---

## 3. Database Generation Pipeline

To reproduce the database construction from raw repertoire data, execute the pipeline scripts sequentially:

### Prerequisites

```bash
# Python dependencies
pip install pyarrow duckdb pandas biopython psutil requests

# R dependencies
R -e "install.packages(c('shiny', 'bslib', 'DT', 'duckdb', 'DBI', 'arrow', 'dplyr', 'htmltools', 'zip'))"
R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install('Biostrings')"
```

### Step 1: Download Raw Repertoire Data from OAS
Run the automated downloader using the metadata catalog:
```bash
python OASpepDB/Scripts/01_download_OAS.py --metadata OASpepDB/Data/OAS_metadata.csv --out-dir /path/to/OAS_raw
```

### Step 2: Index Healthy Control CDR3 Cache (Tier-1 Subtraction)
Extracts and indexes CDR3 sequences from 8,089 healthy repertoires (`Disease == 'None'`) into a compressed Parquet cache:
```bash
python OASpepDB/Scripts/02_build_healthy_cdr3_cache.py \
    --metadata OASpepDB/Data/OAS_metadata.csv \
    --data-dir /path/to/OAS_raw \
    --out-parquet /path/to/CDR3_db/healthy_cdr3_cache.parquet \
    --workers 16
```

### Step 3: Digest Human Reference Proteomes (Tier-2 & Tier-3 Subtraction)
Generates an *in silico* tryptic peptidome (Trypsin/P cleavage, 0–2 missed cleavages, peptide length 6–45 aa) from UniProtKB/Swiss-Prot and NCBI RefSeq:
```bash
Rscript OASpepDB/Scripts/03_insilico_digest.R
```
*Output: `OASpepDB/Data/negative_human_reference_peptides.parquet`.*

### Step 4: Construct Disease-Exclusive Hive Database
Filters disease-exclusive sequences, quantifies patient sharing ($N$ patients) and read coverage, and partitions data by clinical cohort:
```bash
python OASpepDB/Scripts/04_build_disease_db.py \
    --metadata OASpepDB/Data/OAS_metadata.csv \
    --data-dir /path/to/OAS_raw \
    --healthy-cache /path/to/CDR3_db/healthy_cdr3_cache.parquet \
    --negative-ref OASpepDB/Data/negative_human_reference_peptides.parquet \
    --out-dir /path/to/CDR3_db \
    --workers 8
```
*Output directory structure (Hive-partitioned):*
```text
CDR3_db/
├── Disease=COVID-19/
│   └── BSource=PBMC/
│       └── BType=Memory-B-Cells/
│           └── Isotype=IGHG1/*.parquet
└── Disease=HIV/...
```

### Step 5: Generate Camelid VHH Entrapment Library (Empirical FDR)
Queries NCBI Entrez E-utilities for published Camelid VHH sequences, producing isobaric ($I \rightarrow L$) micro-cassettes for target-decoy calibration:
```bash
python OASpepDB/Scripts/05_fetch_ncbi_entrapment.py
```
*Output: `OASpepDB/Data/entrapment_cassettes.fasta`.*

---

## 4. Running the Interactive Explorer

The repository includes a standalone R Shiny application utilizing an in-memory DuckDB engine for rapid querying and multi-format data export.

### 1-Click Launchers:
* **Windows**: Double-click `OASpepDB/DBquery-Windows.bat`
* **macOS**: Double-click `OASpepDB/DBquery-Mac.command`
* **Linux**: Run `bash OASpepDB/DBquery-Linux.sh`

### Command Line / R Console:
```r
# From the repository root
shiny::runApp("OASpepDB", launch.browser = TRUE, port = 8080)

# Or within the OASpepDB directory
setwd("OASpepDB")
shiny::runApp("DBquery.R", launch.browser = TRUE, port = 8080)
```

---

## 5. Web Application Features

1. **Query & Export (Cohort Stratification Dashboard)**:
   - **Clinical Hierarchy Visualization**: Interactive alluvial diagram illustrating sequence distribution across 4 clinical domains, 25 disease cohorts, tissue sources, B-cell subsets, and heavy-chain isotypes.
   - **Clonal Stringency Thresholds**:
     - *Min Patients ($N \ge 1, 2, 3, 5$)*: Restricts results to convergent, multi-patient clonotypes.
     - *Min Reads ($\ge 1, 5, 10, 50$)*: Filters low-abundance sequences to ensure robust sequencing evidence.
   - **Export Formats**:
     - `OAS_<Disease>_<Date>.fasta`: Search-ready database containing disease-exclusive CDR3 micro-cassettes, common laboratory contaminants (cRAP), and ~1% Camelid VHH entrapment controls ($I \rightarrow L$ converted).
     - `OAS_<Disease>_<Date>.parquet`: Tabular dataset containing clonotype identifiers, patient counts, read metrics, full amino acid sequences, and source accessions.
     - `OAS_<Disease>_<Date>.txt`: Audit manifest documenting applied filters and query metadata.

2. **Reverse Lookup (Peptide-to-Clonotype Mapping)**:
   - Queries identified experimental tryptic peptides against the complete 65.5-million sequence database via vectorized DuckDB SQL.
   - Returns matching clonotypes, disease specificity, isotype classifications, complete CDR3 cassettes, and source study records.

---

## 6. Cohort Distribution Summary

| Clinical Category | Diseases / Cohorts Included | Non-Redundant CDR3 Peptides |
| :--- | :--- | :--- |
| **Infectious Diseases** | COVID-19, HIV, CMV-EBV, Influenza, Sepsis, RSV, West Nile, Dengue, Ebola, Hepatitis B, Cholera, Rotavirus, Enterovirus, Lyme, S. pneumoniae, C. difficile | **39,400,000+** |
| **Allergy & Airways** | Allergy-SIT, Allergy-NoSIT, Asthma | **8,900,000+** |
| **Autoimmune & Neuro** | SLE (Lupus), Multiple Sclerosis, Rheumatoid Arthritis, Myasthenia Gravis | **13,300,000+** |
| **Hematology & Tumors** | CLL (Chronic Lymphocytic Leukemia), Melanoma | **3,900,000+** |
| **Total Database** | **25 Disease Cohorts** | **65,510,795** |

---

## 7. Data Availability

* **Catalog and Reference Annotations**: Included directly in this repository under [`OASpepDB/Data/`](OASpepDB/Data/).
* **Raw Repertoire Sequences**: Publicly hosted by the [Observed Antibody Space](https://opig.stats.ox.ac.uk/webapps/oas/).
* **Pre-compiled Parquet Database**: The complete partitioned database (~18 GB compressed Parquet across 25 cohorts) is archived on **Zenodo** (*DOI: [10.5281/zenodo.10561456](https://doi.org/10.5281/zenodo.10561456)*). Extract the archive to `D:/OAS/unpaired/CDR3_db` (or specify a custom path via `DB_ROOT` in `DBquery.R`).

---

## 8. Citation & Contact

If you use OASpepDB in your research, please cite:

> **SDU Immunoinformatics Group** (2026). *OASpepDB: A High-Throughput Disease-Exclusive Antibody CDR3 Peptidome Database for Bottom-Up Immunoproteomics*. University of Southern Denmark (SDU).

For inquiries, issue reports, or data submissions, please submit an issue via GitHub or contact the maintainers.
