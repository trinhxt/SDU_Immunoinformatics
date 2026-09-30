#!/usr/bin/env python3
# ==============================================================================
# OASpepDB: Step 04 - Pure Disease-Exclusive Database Construction Engine
# ==============================================================================
# Target Journal: Nature Communications / Nature Biotechnology
# Author / Lab: SDU_Immunoinformatics
# Date: 2026-09-28
#
# Scientific Purpose:
# -------------------
# Constructs the pure, disease-exclusive antibody database (CDR3_db) from 6,344
# human unpaired OAS disease files (spanning 25 disease cohorts). Every candidate
# read is subjected to a stringent Triple-Tier Negative Subtraction Filter:
#
#   1. Tier 1 (Healthy Repertoire Filter):
#      Rejects any read whose CDR3 amino acid sequence (cdr3_aa) exists in the
#      in-memory index of 294,280,746 unique healthy CDR3s (healthy_cdr3_cache.parquet).
#
#   2. Minimal-Flank Micro-Cassette Construction & Proteolytic Boundary Calibration:
#      Synthesizes the proteolytic minimal-flank micro-cassette:
#        Micro-Cassette = [Minimal_FR3_Flank] + [CDR3] + [FR4] + [Isotype_Tail]
#      The flank eliminates redundant upstream FR3 tryptic cleavage sites (e.g., ...LQMNSLR,
#      ...YYCAR) to prevent false-positive framework PSMs in bottom-up MS searches, while
#      preserving terminal junction cleavage. The isotype tail matches physiological constant
#      domain cleavage (+ASTK for IgG/IgM/Bulk, +ASPTSPK for IgA, +GQPK for Lambda, native J for Kappa).
#
#   3. Tier 2 & 3 (Human Proteome Reference Filter):
#      In-silico tryptic digestion is performed on each micro-cassette. Any candidate
#      whose CDR3-spanning tryptic peptide matches a peptide in the human reference
#      proteome (UniProt Swiss-Prot + NCBI RefSeq, 3,063,177 tryptic peptides) is
#      immediately eliminated to prevent false cross-annotation of host proteins.
#
#   4. Global Patient Traceability & Clonal Aggregation:
#      Subject IDs across studies are disambiguated using study-prefixed Patient_ID
#      ({Author}_{Subject}) to prevent cross-study patient collisions. For each
#      unique cassette, distinct patient recurrence (N_Patients) and full provenance
#      (Patient_IDs) are computed alongside total read redundancy (Redundancy).
#
#   5. Standardized 4-Tier Hive Partitioning:
#      Disease names are standardized with Disease_short (e.g. COVID-19, CLL, HIV,
#      Celiac-Disease, AL-Amyloidosis) to guarantee slash-free, URI-safe directory
#      structures while embedding the complete clinical designation (Disease_long).
#      Output directory hierarchy:
#        Disease=<Disease_short>/BSource=<BSource>/BType=<BType>/Isotype=<Isotype>/data_0.parquet
#      compressed with Zstandard (ZSTD) Level 9 for sub-0.05s DuckDB partition pruning.
# ==============================================================================

import argparse
import gzip
import hashlib
import os
import sys
import time
import urllib.parse
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd
import psutil
import pyarrow as pa
import pyarrow.csv as pacsv
import pyarrow.parquet as pq


# ------------------------------------------------------------------------------
# Helper Functions: Biological Cleavage & Hive Path Formulation
# ------------------------------------------------------------------------------
def quote_hive(val: str) -> str:
    """
    Escapes reserved path characters for standard Apache Hive partitioning compatibility.
    
    Example: 'Celiac-Disease' -> 'Celiac-Disease'
    """
    return urllib.parse.quote(str(val).strip(), safe="-_.")


def get_isotype_tail(isotype: str, chain: str, v_call: str = "", j_call: str = "") -> str:
    """
    Determines the physiological C-terminal tryptic extension tail for the micro-cassette.
    
    Rules:
      - Light Kappa (IGK): Native unextended J-terminus (empty string "")
      - Light Lambda (IGL): +GQPK (first 4 aa of IGLC constant region)
      - Heavy IgA (IGHA): +ASPTSPK (first 7 aa of IGHA constant region)
      - Heavy IgG, IgM, IgD, IgE, Bulk, or unspecified Heavy: +ASTK
    """
    iso_up = str(isotype).upper()
    chain_up = str(chain).upper()
    v_up = str(v_call).upper()
    j_up = str(j_call).upper()

    # Distinguish light chain classes (Kappa vs Lambda)
    if chain_up == "LIGHT" or "IGK" in v_up or "IGL" in v_up:
        if "IGL" in v_up or "IGL" in j_up:
            return "GQPK"
        elif "IGK" in v_up or "IGK" in j_up:
            return ""
        else:
            # Default light chain fallback to Kappa (predominant in human repertoires)
            return ""

    # Heavy chain isotype assignment
    if iso_up == "IGHA":
        return "ASPTSPK"
    else:
        # Standard IgG, IgM, IgD, IgE, Bulk, or default Heavy
        return "ASTK"


FRAMEWORK_BLACKLIST = {
    "WGQGTLVTVSSASTK",
    "GTLVTVSSASTK",
    "GTLVTVSSASPTSPK",
    "GTTVTVSSASTK",
    "GSLVTVSSASTK",
    "LTFGGGTK",
    "SLTFGGGTK",
    "AEDTATYYCAR",
    "AEDTALYYCAR",
    "SEDTATYYCAR",
    "VEDTAVYYCAR",
}


def build_minimal_cassette(fr3: str, cdr3: str, fwr4: str, tail: str) -> tuple[str, int]:
    """
    Constructs a minimal-flank micro-cassette engineered for bottom-up mass spectrometry.

    Eliminates redundant upstream Framework 3 tryptic cleavage sites (e.g. ...LQMNSLR,
    ...YYCAR) that generate non-specific shared framework PSMs during database searches:
      - If CDR3 starts with basic residues at the junction (e.g. AR, AK, VR, TR),
        tryptic cleavage occurs precisely after position 106. No upstream FR3 flank
        is required; the cassette retains only the conserved C104 anchor.
      - If CDR3 lacks a basic residue at position 106 (e.g. AS, AT, AG), cleavage
        relies on the C-terminal-most basic residue of FR3 (e.g. ...AYMELSSLR).
        Only the terminal non-cleavable segment of FR3 (e.g. SEDTAVYYC) is preserved.

    Parameters:
      fr3: Full Framework 3 sequence from AIRR annotation (~25 aa)
      cdr3: Core CDR3 loop sequence
      fwr4: Framework 4 sequence
      tail: Constant region extension tail (e.g. ASTK, ASPTSPK, GQPK)

    Returns:
      tuple of (minimal_cassette_sequence, cdr3_start_index)
    """
    prefix = cdr3[:2].upper() if len(cdr3) >= 2 else ""

    # Case 1: Tryptic cleavage site exists at the N-terminal anchor of CDR3 (AR, AK, etc.)
    if prefix in ("AR", "AK", "VR", "TR", "PR", "GR", "IR", "LR"):
        min_cas = "C" + cdr3 + fwr4 + tail
        cdr3_start = 1
        return min_cas, cdr3_start

    # Case 2: No basic residue at junction; find the final K/R cleavage in FR3
    last_cleavage = -1
    for idx in range(len(fr3)):
        if fr3[idx] in ("K", "R"):
            if idx + 1 < len(fr3) and fr3[idx + 1] == "P":
                continue
            last_cleavage = idx

    if last_cleavage != -1 and last_cleavage + 1 < len(fr3):
        minimal_fr3 = fr3[last_cleavage + 1 :]
    else:
        # Fallback: retain at most the last 9 residues of FR3
        minimal_fr3 = fr3[-9:] if len(fr3) >= 9 else fr3

    min_cas = minimal_fr3 + cdr3 + fwr4 + tail
    cdr3_start = len(minimal_fr3)
    return min_cas, cdr3_start


def extract_cdr3_tryptic_peptides(
    cassette: str, cdr3_start: int, cdr3_len: int
) -> list[str]:
    """
    Performs in-silico tryptic cleavage (cleaves at carboxyl side of Lys/Arg,
    unless followed by Proline) and extracts ONLY bona fide CDR3-bearing tryptic
    peptides conforming to 4-Pillar Nature Portfolio standards:
      1. Proteotypic length: 7 to 40 amino acids.
      2. Framework blacklist exclusion: Rejects pure framework/tail fragments.
      3. Hypervariable overlap: Covers at least 5 consecutive amino acids in
         the CDR3 core, or >= 50% of the loop length for short CDR3s.

    Parameters:
      cassette: Micro-cassette amino acid sequence
      cdr3_start: 0-based start index of the CDR3 loop within the cassette
      cdr3_len: Length of the core CDR3 loop

    Returns:
      List of unique qualified peptide sequences spanning the CDR3 core.
    """
    n = len(cassette)
    cdr3_end = cdr3_start + cdr3_len

    peptides = []
    start = 0
    for i in range(n):
        if cassette[i] in ("K", "R"):
            if i + 1 < n and cassette[i + 1] == "P":
                continue  # Proline rule: cleavage inhibited
            peptides.append((cassette[start : i + 1], start, i + 1))
            start = i + 1
    if start < n:
        peptides.append((cassette[start:], start, n))

    qualified = []
    for p, p_start, p_end in peptides:
        # Criterion 1: Mass spectrometry proteotypic length window (7 - 40 aa)
        if not (7 <= len(p) <= 40):
            continue

        # Criterion 2: Framework blacklist exclusion
        if p in FRAMEWORK_BLACKLIST:
            continue

        # Criterion 3: CDR3 hypervariable overlap gate
        overlap_start = max(p_start, cdr3_start)
        overlap_end = min(p_end, cdr3_end)
        overlap_len = max(0, overlap_end - overlap_start)

        if overlap_len >= 5 or (cdr3_len > 0 and overlap_len / cdr3_len >= 0.5):
            qualified.append(p)

    return qualified


# ------------------------------------------------------------------------------
# Loading Reference Indexes
# ------------------------------------------------------------------------------
def load_healthy_cdr3_cache(cache_path: Path) -> set[str]:
    """
    Loads the Tier 1 negative background cache (294,280,746 unique healthy CDR3s)
    from Apache Parquet into an in-memory Python hash set for instant O(1) lookups.
    """
    print(f"\n[STEP 1/5] Loading Tier 1 Healthy CDR3 Cache from Parquet...")
    print(f"  Source file: {cache_path}")
    t0 = time.time()

    if not cache_path.is_file():
        raise FileNotFoundError(f"Healthy CDR3 cache not found at: {cache_path}")

    healthy_set = set()
    pf = pq.ParquetFile(cache_path)
    total_rows = pf.metadata.num_rows
    loaded_rows = 0

    # Stream in 10M record batches to monitor memory and progress
    for batch in pf.iter_batches(batch_size=10_000_000, columns=["cdr3_aa"]):
        col_list = batch["cdr3_aa"].to_pylist()
        healthy_set.update(col_list)
        loaded_rows += len(col_list)
        pct = (loaded_rows / total_rows) * 100
        ram_gb = psutil.virtual_memory().used / (1024**3)
        print(
            f"  -> Ingested {loaded_rows:,}/{total_rows:,} records ({pct:5.1f}%) | "
            f"RAM: {ram_gb:.1f} GB | Elapsed: {time.time()-t0:.1f}s",
            end="\r",
            flush=True,
        )

    elapsed = time.time() - t0
    final_ram_gb = psutil.virtual_memory().used / (1024**3)
    print(f"\n  Tier 1 Cache Loaded: {len(healthy_set):,} unique healthy CDR3s in {elapsed:.1f}s.")
    print(f"  System RAM Usage: {final_ram_gb:.2f} GB (Ready for thread-safe O(1) lookups).")
    return healthy_set


def load_negative_human_reference(ref_path: Path) -> set[str]:
    """
    Loads the Tier 2 & 3 negative human reference peptides (3,063,177 tryptic
    peptides from UniProt Swiss-Prot + NCBI RefSeq) into an in-memory hash set.
    """
    print(f"\n[STEP 2/5] Loading Tier 2 & 3 Human Proteome Reference Peptides...")
    print(f"  Source file: {ref_path}")
    t0 = time.time()

    if not ref_path.is_file():
        raise FileNotFoundError(f"Human reference peptides file not found at: {ref_path}")

    table = pq.read_table(ref_path, columns=["peptide"])
    peptides = table["peptide"].to_pylist()
    ref_set = set(peptides)

    elapsed = time.time() - t0
    print(f"  Tier 2 & 3 Reference Loaded: {len(ref_set):,} human tryptic peptides in {elapsed:.2f}s.")
    return ref_set


# ------------------------------------------------------------------------------
# Partition Worker Function
# ------------------------------------------------------------------------------
def process_single_partition(
    partition_key: tuple[str, str, str, str],
    file_records: list[dict],
    healthy_set: set[str],
    ref_set: set[str],
    data_dir: Path,
    out_root: Path,
) -> dict:
    """
    Processes all disease files associated with a single (Disease_short, BSource, BType, Isotype)
    partition. Applies Tier 1, 2, and 3 negative subtraction filters, builds micro-cassettes,
    aggregates clonal statistics with cross-study Patient_ID deduplication, and writes a
    compressed Parquet partition file.
    
    Returns a dictionary of execution metrics for the partition.
    """
    disease_short, bsource, btype, isotype = partition_key
    total_reads_scanned = 0
    tier1_dropped_count = 0
    tier23_dropped_count = 0

    # Retrieve full clinical disease name from metadata records
    disease_long = file_records[0]["Disease_long"] if file_records else disease_short

    # In-memory clonal aggregation map: cassette -> clonal_record
    # Structure:
    # {
    #   'cdr3_aa': str,
    #   'fr3_tail': str,
    #   'fwr4_aa': str,
    #   'isotype_tail': str,
    #   'tryptic_peptides': str,
    #   'v_call': str,
    #   'd_call': str,
    #   'j_call': str,
    #   'v_ident_sum': float,
    #   'v_ident_count': int,
    #   'sequence_alignment_aa': str,
    #   'redundancy': int,
    #   'patient_ids': set[str]  # Disambiguated cross-study Patient_IDs
    # }
    clones = {}

    for record in file_records:
        fname = record["Filename"]
        patient_id = record["Patient_ID"]  # Study-disambiguated unique patient identifier
        chain = record["Chain"]
        filepath = data_dir / fname

        if not filepath.is_file():
            continue

        try:
            # Fast C++ SIMD CSV reader via PyArrow
            table = pacsv.read_csv(filepath)
            n_rows = len(table)
            total_reads_scanned += n_rows

            if n_rows == 0:
                continue

            # Extract columnar arrays as Python lists for low-overhead row iteration
            cdr3_col = table["cdr3_aa"].to_pylist()
            fr3_col = table["fr3_tail"].to_pylist()
            fwr4_col = table["fwr4_aa"].to_pylist()
            v_call_col = table["v_call"].to_pylist()
            d_call_col = table["d_call"].to_pylist()
            j_call_col = table["j_call"].to_pylist()
            v_ident_col = table["v_identity"].to_pylist()
            seq_align_col = table["sequence_alignment_aa"].to_pylist()
            redundancy_col = table["Redundancy"].to_pylist()

            for i in range(n_rows):
                cdr3 = cdr3_col[i]
                if not cdr3 or not (5 <= len(cdr3) <= 40) or not cdr3.isalpha():
                    continue

                # --------------------------------------------------------------
                # FILTER 1 (Tier 1 Healthy Subtraction)
                # --------------------------------------------------------------
                if cdr3 in healthy_set:
                    tier1_dropped_count += 1
                    continue

                fr3 = fr3_col[i] or ""
                fwr4 = fwr4_col[i] or ""
                v_call = v_call_col[i] or ""
                j_call = j_call_col[i] or ""
                d_call = d_call_col[i] or ""
                v_ident = float(v_ident_col[i]) if v_ident_col[i] is not None else 100.0
                seq_align = seq_align_col[i] or ""
                red = int(redundancy_col[i]) if redundancy_col[i] is not None else 1

                # Framework sanity checks
                if len(fr3) < 20 or len(fwr4) < 5:
                    continue

                # Determine constant extension tail
                tail = get_isotype_tail(isotype, chain, v_call, j_call)

                # Assemble bottom-up MS engineered Minimal-Flank Cassette
                cassette, cdr3_start_idx = build_minimal_cassette(fr3, cdr3, fwr4, tail)

                # --------------------------------------------------------------
                # FILTER 2 & 3 (Tier 2 & 3 Proteome Reference Subtraction)
                # --------------------------------------------------------------
                tryptic_cdr3_peps = extract_cdr3_tryptic_peptides(
                    cassette, cdr3_start_idx, len(cdr3)
                )

                # Disqualify cassettes yielding no detectable proteotypic CDR3 tryptic peptide
                if not tryptic_cdr3_peps:
                    continue

                # Check if any bona fide CDR3 peptide matches the human reference proteome
                if any(p in ref_set for p in tryptic_cdr3_peps):
                    tier23_dropped_count += 1
                    continue

                # --------------------------------------------------------------
                # Clonal Aggregation within Partition
                # --------------------------------------------------------------
                tryptic_str = ";".join(tryptic_cdr3_peps)
                if cassette in clones:
                    entry = clones[cassette]
                    entry["redundancy"] += red
                    entry["patient_ids"].add(patient_id)
                    entry["data_files"].add(fname)
                    entry["v_ident_sum"] += v_ident
                    entry["v_ident_count"] += 1
                    # Retain the longest complete VH alignment as representative
                    if len(seq_align) > len(entry["sequence_alignment_aa"]):
                        entry["sequence_alignment_aa"] = seq_align
                else:
                    clones[cassette] = {
                        "cdr3_aa": cdr3,
                        "fr3_tail": fr3,
                        "fwr4_aa": fwr4,
                        "isotype_tail": tail,
                        "tryptic_peptides": tryptic_str,
                        "v_call": v_call,
                        "d_call": d_call,
                        "j_call": j_call,
                        "v_ident_sum": v_ident,
                        "v_ident_count": 1,
                        "sequence_alignment_aa": seq_align,
                        "redundancy": red,
                        "patient_ids": {patient_id},
                        "data_files": {fname},
                    }

        except Exception as e:
            print(f"[WARN] Error reading {fname}: {e}", file=sys.stderr)
            continue

    n_unique_cassettes = len(clones)
    out_file_size = 0

    # --------------------------------------------------------------------------
    # Write Partition to Apache Parquet (ZSTD Level 9)
    # --------------------------------------------------------------------------
    if n_unique_cassettes > 0:
        clonotype_ids_list = []
        cassettes_list = []
        cdr3_list = []
        cdr3_len_list = []
        fr3_list = []
        fwr4_list = []
        tail_list = []
        tryptic_list = []
        v_call_list = []
        d_call_list = []
        j_call_list = []
        v_ident_list = []
        seq_align_list = []
        redundancy_list = []
        n_patients_list = []
        data_files_list = []
        disease_long_list = []

        for cas, d in clones.items():
            # Deterministic, unique MS Accession ID (e.g. OAS_COVID-19_8C967109)
            cid = f"OAS_{quote_hive(disease_short)}_{hashlib.md5(cas.encode()).hexdigest()[:8].upper()}"
            clonotype_ids_list.append(cid)
            cassettes_list.append(cas)
            cdr3_list.append(d["cdr3_aa"])
            cdr3_len_list.append(len(d["cdr3_aa"]))
            fr3_list.append(d["fr3_tail"])
            fwr4_list.append(d["fwr4_aa"])
            tail_list.append(d["isotype_tail"])
            tryptic_list.append(d["tryptic_peptides"])
            v_call_list.append(d["v_call"])
            d_call_list.append(d["d_call"])
            j_call_list.append(d["j_call"])
            avg_v_ident = round(d["v_ident_sum"] / d["v_ident_count"], 3) if d["v_ident_count"] > 0 else 100.0
            v_ident_list.append(avg_v_ident)
            seq_align_list.append(d["sequence_alignment_aa"])
            redundancy_list.append(d["redundancy"])
            n_patients_list.append(len(d["patient_ids"]))
            data_files_list.append(";".join(sorted(d["data_files"])))
            disease_long_list.append(disease_long)

        arrow_table = pa.Table.from_arrays(
            [
                pa.array(clonotype_ids_list, type=pa.string()),
                pa.array(cassettes_list, type=pa.string()),
                pa.array(cdr3_list, type=pa.string()),
                pa.array(cdr3_len_list, type=pa.int16()),
                pa.array(fr3_list, type=pa.string()),
                pa.array(fwr4_list, type=pa.string()),
                pa.array(tail_list, type=pa.string()),
                pa.array(tryptic_list, type=pa.string()),
                pa.array(v_call_list, type=pa.string()),
                pa.array(d_call_list, type=pa.string()),
                pa.array(j_call_list, type=pa.string()),
                pa.array(v_ident_list, type=pa.float32()),
                pa.array(seq_align_list, type=pa.string()),
                pa.array(redundancy_list, type=pa.int64()),
                pa.array(n_patients_list, type=pa.int32()),
                pa.array(data_files_list, type=pa.string()),
                pa.array(disease_long_list, type=pa.string()),
            ],
            names=[
                "clonotype_id",
                "cdr3_cassette",
                "cdr3_aa",
                "cdr3_len",
                "fr3_tail",
                "fwr4_aa",
                "isotype_tail",
                "tryptic_peptides",
                "v_call",
                "d_call",
                "j_call",
                "v_identity",
                "sequence_alignment_aa",
                "Redundancy",
                "N_Patients",
                "Data_file",
                "Disease_long",
            ],
        )

        # Build 4-tier Hive directory structure using clean Disease_short
        partition_dir = (
            out_root
            / f"Disease={quote_hive(disease_short)}"
            / f"BSource={quote_hive(bsource)}"
            / f"BType={quote_hive(btype)}"
            / f"Isotype={quote_hive(isotype)}"
        )
        partition_dir.mkdir(parents=True, exist_ok=True)
        out_parquet = partition_dir / "data_0.parquet"

        # Write Parquet with ZSTD level 9
        pq.write_table(
            arrow_table,
            out_parquet,
            compression="zstd",
            compression_level=9,
            use_dictionary=True,
        )
        out_file_size = out_parquet.stat().st_size

    return {
        "partition": partition_key,
        "n_files": len(file_records),
        "reads_scanned": total_reads_scanned,
        "tier1_dropped": tier1_dropped_count,
        "tier23_dropped": tier23_dropped_count,
        "unique_cassettes": n_unique_cassettes,
        "file_size": out_file_size,
    }


# ------------------------------------------------------------------------------
# Main Execution Pipeline
# ------------------------------------------------------------------------------
def main():
    script_dir = Path(__file__).resolve().parent
    default_meta = str((script_dir.parent / "Data" / "OAS_metadata.csv").resolve())
    default_neg_ref = str((script_dir.parent / "Data" / "negative_human_reference_peptides.parquet").resolve())

    # Dynamically select default paths
    candidate_data = Path(r"D:\OAS\unpaired\download")
    default_data_dir = str(candidate_data) if candidate_data.is_dir() else str((script_dir.parent / "Data" / "raw").resolve())

    candidate_cache = Path(r"D:\OAS\unpaired\cache\healthy_cdr3_cache.parquet")
    default_healthy = str(candidate_cache) if candidate_cache.is_file() else str((script_dir.parent / "Data" / "healthy_cdr3_cache.parquet").resolve())

    candidate_out = Path(r"D:\OAS\unpaired\CDR3_db")
    default_out_dir = str(candidate_out) if candidate_out.is_dir() else str((script_dir.parent / "CDR3_db").resolve())

    parser = argparse.ArgumentParser(
        description="OASpepDB Step 04: Build Pure Disease-Exclusive Antibody Database (CDR3_db)."
    )
    parser.add_argument(
        "--metadata",
        "--meta-csv",
        dest="metadata",
        type=str,
        default=default_meta,
        help="Path to OAS_metadata.csv catalog",
    )
    parser.add_argument(
        "--data-dir",
        "--input-dir",
        dest="data_dir",
        type=str,
        default=default_data_dir,
        help="Directory containing the disease CSV.gz files",
    )
    parser.add_argument(
        "--healthy-cache",
        type=str,
        default=default_healthy,
        help="Path to Tier 1 Healthy CDR3 cache Parquet file",
    )
    parser.add_argument(
        "--negative-ref",
        type=str,
        default=default_neg_ref,
        help="Path to Tier 2 & 3 Human Proteome negative reference Parquet file",
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=default_out_dir,
        help="Root output directory for Hive-partitioned CDR3_db",
    )
    parser.add_argument(
        "--workers",
        "--threads",
        dest="workers",
        type=int,
        default=min(8, os.cpu_count() or 4),
        help="Number of concurrent worker threads (default: min(8, CPU count))",
    )
    parser.add_argument(
        "--disease",
        type=str,
        default=None,
        help="Optional: Filter execution to a single specific disease (e.g. 'COVID-19'). If omitted, processes all diseases.",
    )

    args = parser.parse_args()

    meta_path = Path(args.metadata)
    data_dir = Path(args.data_dir)
    healthy_cache_path = Path(args.healthy_cache)
    neg_ref_path = Path(args.negative_ref)
    out_dir = Path(args.out_dir)

    print("=" * 80)
    print("  OASpepDB: PURE DISEASE-EXCLUSIVE ANTIBODY DATABASE BUILDER (CDR3_db)")
    print("=" * 80)
    print(f"Catalog Metadata:   {meta_path}")
    print(f"Raw Disease Data:   {data_dir}")
    print(f"Tier 1 Healthy:     {healthy_cache_path}")
    print(f"Tier 2 & 3 Proteome:{neg_ref_path}")
    print(f"Output Database:    {out_dir}")
    print(f"Worker Threads:     {args.workers} (CPU cores available: {os.cpu_count()})")
    print(f"System RAM:         {psutil.virtual_memory().total / (1024**3):.1f} GB")
    print("=" * 80)

    # 1. Load Tier 1 Cache
    healthy_set = load_healthy_cdr3_cache(healthy_cache_path)

    # 2. Load Tier 2 & 3 Reference
    ref_set = load_negative_human_reference(neg_ref_path)

    # 3. Parse Metadata and Partition Disease Cohort
    print(f"\n[STEP 3/5] Parsing metadata and organizing 4-tier Hive partitions...")
    df_meta = pd.read_csv(meta_path, keep_default_na=False)
    
    # Isolate non-healthy files (Disease != 'None')
    disease_df = df_meta[df_meta["Disease"] != "None"].copy()
    if args.disease:
        disease_df = disease_df[disease_df["Disease_short"].str.upper() == args.disease.strip().upper()].copy()
        print(f"  Filtering exclusively for Disease cohort: '{args.disease}' ({len(disease_df):,} files).")
    else:
        print(f"  Identified {len(disease_df):,} disease files across {disease_df['Disease_short'].nunique()} disease cohorts.")

    # Group files by 4-tier partition tuple (Disease_short, BSource, BType, Isotype)
    partition_groups = defaultdict(list)
    total_expected_files = 0

    for _, row in disease_df.iterrows():
        pkey = (row["Disease_short"], row["BSource"], row["BType"], row["Isotype"])
        partition_groups[pkey].append({
            "Filename": row["Filename"],
            "Patient_ID": str(row["Patient_ID"]).strip(),
            "Disease_long": str(row["Disease_long"]).strip(),
            "Chain": str(row["Chain"]).strip(),
        })
        total_expected_files += 1

    total_partitions = len(partition_groups)
    print(f"  Organized into {total_partitions:,} distinct Hive partitions (using standardized Disease_short).")

    # 4. Multi-Threaded Partition Processing
    print(f"\n[STEP 4/5] Executing multi-threaded partition ETL ({args.workers} parallel workers)...")
    start_time = time.time()
    
    completed_partitions = 0
    total_reads_scanned = 0
    total_tier1_dropped = 0
    total_tier23_dropped = 0
    total_unique_disease_cassettes = 0
    total_written_bytes = 0

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_partition = {
            executor.submit(
                process_single_partition,
                pkey,
                files,
                healthy_set,
                ref_set,
                data_dir,
                out_dir,
            ): pkey
            for pkey, files in partition_groups.items()
        }

        for future in as_completed(future_to_partition):
            res = future.result()
            completed_partitions += 1
            total_reads_scanned += res["reads_scanned"]
            total_tier1_dropped += res["tier1_dropped"]
            total_tier23_dropped += res["tier23_dropped"]
            total_unique_disease_cassettes += res["unique_cassettes"]
            total_written_bytes += res["file_size"]

            # Real-time console progress
            elapsed = time.time() - start_time
            part_per_sec = completed_partitions / elapsed if elapsed > 0 else 0
            pct = (completed_partitions / total_partitions) * 100
            rem_partitions = total_partitions - completed_partitions
            eta_sec = rem_partitions / part_per_sec if part_per_sec > 0 else 0
            ram_used_gb = psutil.virtual_memory().used / (1024**3)

            print(
                f"  [{pct:5.1f}%] {completed_partitions:,}/{total_partitions:,} partitions | "
                f"Disease Cassettes: {total_unique_disease_cassettes:,} | "
                f"Scanned: {total_reads_scanned:,} | "
                f"Speed: {part_per_sec:4.1f} part/s | "
                f"RAM: {ram_used_gb:.1f} GB | "
                f"ETA: {eta_sec:3.0f}s",
                end="\r",
                flush=True,
            )

    total_elapsed = time.time() - start_time

    # 5. Final Consolidated Summary
    print(f"\n\n" + "=" * 80)
    print("  PURE DISEASE-EXCLUSIVE DATABASE (CDR3_db) BUILT SUCCESSFULLY!")
    print("=" * 80)
    print(f"  Total Partitions Built:     {completed_partitions:,} / {total_partitions:,}")
    print(f"  Disease Files Processed:    {total_expected_files:,} files")
    print(f"  Total NGS Reads Scanned:    {total_reads_scanned:,} reads")
    print(f"  Filtered by Tier 1 Healthy: {total_tier1_dropped:,} reads ({total_tier1_dropped/total_reads_scanned*100:.1f}%)" if total_reads_scanned > 0 else "  Filtered by Tier 1: 0")
    print(f"  Filtered by Tier 2/3 Human: {total_tier23_dropped:,} reads")
    print(f"  Disease Micro-Cassettes:    {total_unique_disease_cassettes:,} non-redundant clonotypes")
    print(f"  Total Database Disk Size:   {total_written_bytes / (1024**2):.2f} MB ({total_written_bytes / (1024**3):.2f} GB)")
    print(f"  Total Execution Time:       {total_elapsed:.1f}s ({total_elapsed/60:.2f} min)")
    print(f"  Storage Format:             4-tier Hive Apache Parquet (Zstandard Level 9)")
    print(f"  Target Directory:           {out_dir}")
    print("=" * 80)


if __name__ == "__main__":
    main()
