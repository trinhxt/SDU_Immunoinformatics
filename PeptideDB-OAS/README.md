# PeptideDB-OAS

A reproducible R pipeline (Part1–Part6) for building a **peptide-level antibody database**
from the Observed Antibody Space (OAS), enriched with:

- Disease presence per peptide
- Antibody uniqueness per peptide
- Binary annotation of peptide occurrence within CDR3 regions

The pipeline is designed for **large-scale datasets** (hundreds of millions of rows),
using **Arrow + DuckDB**, and produces a **Zenodo-ready Parquet database**.

---

## Overview of the pipeline

| Part | Description |
|-----:|-------------|
| Part 1 | Download OAS files, extract minimal antibody columns, build metadata |
| Part 2 | In-silico tryptic digestion → peptide–antibody Parquet database |
| Part 3 | Disease presence (distinct diseases per antibody / peptide) |
| Part 4 | Peptide uniqueness (distinct antibodies per peptide) |
| Part 5 | Peptide ∈ CDR3 annotation (row-consistent logic) |
| Part 6 | Build a single Zenodo-shareable PeptideDB folder |

---

## Repository structure

```text
PeptideDB-OAS/
├─ scripts/        # Part1–Part6 pipeline steps
├─ R/              # Shared helper functions
├─ config/         # Configuration templates
├─ docs/           # Pipeline and data model documentation
├─ renv.lock       # Reproducible R environment
├─ run_all.R       # Single entry point
└─ README.md
