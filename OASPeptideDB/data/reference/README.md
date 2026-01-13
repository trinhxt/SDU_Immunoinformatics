# Reference tryptic peptide set (UniProt + NCBI)

This directory contains instructions for generating the reference tryptic
peptide set used in **Part2 (digestion)** of the PeptideDB-OAS pipeline.

⚠️ Large reference data files are **not included** in this repository.

---

## Purpose

The reference peptide set is used to **filter out peptides** that are present
in common human proteomes, reducing background and non–antibody-specific
signals.

Specifically, peptides derived from **human UniProt and NCBI RefSeq proteins**
are removed from the antibody-derived peptide pool.

---

## Source data (human proteome)

The reference was built from the following FASTA files:

1. **UniProt (reviewed + unreviewed)**
   - File: `UniProt_TR_SP_Human_2024_07_25.fasta`
   - Source: https://www.uniprot.org
   - Includes: Swiss-Prot + TrEMBL human entries

2. **NCBI RefSeq (human)**
   - File: `NCBI_RefSeq_Human_2024_07_25.fasta`
   - Source: https://www.ncbi.nlm.nih.gov/refseq/

---

## Processing steps

1. Combine UniProt and NCBI human protein FASTA files
2. Perform **in-silico tryptic digestion**
   - Enzyme: Trypsin
   - Missed cleavages: 0–1
3. Extract **distinct peptide sequences**
4. Save as an R object:

```r
UniProtNCBI_Tryptic  # character vector of unique tryptic peptides
