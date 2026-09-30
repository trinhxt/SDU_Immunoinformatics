#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
OAS Paired Metadata Extractor
================================================================================
Fast, lightweight tool to scan Line 1 JSON metadata from all paired antibody
files (*.csv.gz) in OAS paired repository and compile OAS_paired_metadata.csv.

Default Paths:
  Input raw files:  D:\\OAS\\paired\\raw
  Output metadata:  D:\\OAS\\paired\\OAS_paired_metadata.csv
================================================================================
"""

import os
import sys
import csv
import json
import gzip
import time
import argparse
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed


# Canonical order of columns for clean readability
CANONICAL_COLUMNS = [
    "Run",
    "Link",
    "Author",
    "Species",
    "BSource",
    "BType",
    "Longitudinal",
    "Disease",
    "Subject",
    "Age",
    "Vaccine",
    "Chain",
    "Unique sequences",
    "Isotype",
    "Filename",
]


def extract_header_from_file(filepath: Path) -> dict:
    """
    Reads only Line 1 from a gzip file, parses CSV-wrapped JSON metadata.
    Returns dict of metadata with 'Filename' key, or None if invalid.
    """
    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="ignore") as f:
            line1 = f.readline().strip()
            if not line1:
                return None
            
            # Line 1 is a CSV row enclosing the JSON object
            reader = csv.reader([line1])
            json_str = next(reader)[0]
            data = json.loads(json_str)
            data["Filename"] = filepath.name
            return data
    except Exception as e:
        print(f"[WARN] Failed to parse {filepath.name}: {e}", file=sys.stderr)
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Extract Line 1 JSON metadata from OAS paired .csv.gz files."
    )
    parser.add_argument(
        "--raw-dir",
        type=str,
        default=r"D:\OAS\paired\raw",
        help="Directory containing raw paired .csv.gz files (default: D:\\OAS\\paired\\raw)",
    )
    parser.add_argument(
        "--out-csv",
        type=str,
        default=r"D:\OAS\paired\OAS_paired_metadata.csv",
        help="Path for output metadata CSV (default: D:\\OAS\\paired\\OAS_paired_metadata.csv)",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=16,
        help="Number of concurrent worker threads (default: 16)",
    )
    args = parser.parse_args()

    raw_dir = Path(args.raw_dir)
    out_csv = Path(args.out_csv)

    print("=" * 80)
    print("  OAS PAIRED METADATA EXTRACTOR")
    print("=" * 80)
    print(f"Source Directory: {raw_dir}")
    print(f"Output Metadata:  {out_csv}")
    print(f"Worker Threads:   {args.workers}")
    print("=" * 80)

    if not raw_dir.is_dir():
        print(f"[ERROR] Directory not found: {raw_dir}", file=sys.stderr)
        sys.exit(1)

    # 1. Discover all .csv.gz files
    gz_files = sorted(list(raw_dir.glob("*.csv.gz")))
    total_files = len(gz_files)
    print(f"[INFO] Found {total_files:,} .csv.gz files in {raw_dir}")

    if total_files == 0:
        print("[WARN] No .csv.gz files found to process.")
        sys.exit(0)

    start_time = time.time()
    records = []
    discovered_keys = set()

    # 2. Extract metadata in parallel
    print(f"[INFO] Extracting Line 1 metadata using {args.workers} threads...")
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        future_to_file = {
            executor.submit(extract_header_from_file, f): f for f in gz_files
        }
        done_count = 0
        for future in as_completed(future_to_file):
            result = future.result()
            done_count += 1
            if result:
                records.append(result)
                discovered_keys.update(result.keys())
            
            if done_count % 100 == 0 or done_count == total_files:
                pct = (done_count / total_files) * 100
                print(f"  Progress: {done_count:,}/{total_files:,} files ({pct:.1f}%)", end="\r")
    
    elapsed = time.time() - start_time
    print(f"\n[INFO] Extraction finished in {elapsed:.2f} seconds ({len(records):,}/{total_files:,} valid).")

    # 3. Determine column ordering: Canonical columns first, then any extra discovered keys
    extra_keys = sorted([k for k in discovered_keys if k not in CANONICAL_COLUMNS])
    final_columns = [k for k in CANONICAL_COLUMNS if k in discovered_keys] + extra_keys

    # 4. Sort records by Run (numeric if possible), then Filename
    def sort_key(rec):
        run_val = rec.get("Run", 0)
        try:
            return (int(run_val), rec.get("Filename", ""))
        except (ValueError, TypeError):
            return (99999999, rec.get("Filename", ""))

    records.sort(key=sort_key)

    # 5. Write to CSV
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=final_columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(records)

    file_size_kb = out_csv.stat().st_size / 1024
    print(f"[SUCCESS] Successfully saved metadata to {out_csv} ({file_size_kb:.1f} KB)")

    # 6. Summary breakdown
    diseases = set(r.get("Disease", "Unknown") for r in records)
    subjects = set(r.get("Subject", "Unknown") for r in records)
    species = set(r.get("Species", "Unknown") for r in records)
    total_unique_seqs = sum(
        int(r.get("Unique sequences", 0))
        for r in records
        if str(r.get("Unique sequences", "")).isdigit()
    )

    print("-" * 80)
    print("METADATA SUMMARY:")
    print(f"  Total Paired Files:     {len(records):,}")
    print(f"  Total Unique Sequences: {total_unique_seqs:,}")
    print(f"  Distinct Species:       {len(species)} ({', '.join(sorted(species))})")
    print(f"  Distinct Diseases:      {len(diseases)} ({', '.join(sorted(diseases)[:10])}{'...' if len(diseases) > 10 else ''})")
    print(f"  Distinct Subjects:      {len(subjects):,}")
    print("=" * 80)


if __name__ == "__main__":
    main()
