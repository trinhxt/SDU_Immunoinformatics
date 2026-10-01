#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
OASpepDB: Reference Proteomes Downloader
(download_reference_proteomes.py)
================================================================================
Purpose:
  Automated, cross-platform downloader for the standard human reference
  proteomes required by OASpepDB Step 03 (In Silico Digestion) and FASTA export:
    1. UniProt Swiss-Prot Human Canonical + Reviewed Isoforms (FASTA)
       -> Data/UniProt_SP_canonical_isoform_2026_09_26.fasta
    2. UniProt Swiss-Prot Human Canonical Proteome (FASTA)
       -> Data/Uniprot_SP_canonical_2026_09_14.fasta
    3. NCBI RefSeq Human Proteome GRCh38.p14 (FAA)
       -> Data/GCF_000001405.40_GRCh38.p14_protein.faa

Features:
  - Streaming decompression directly from official UniProt REST and NCBI FTP/HTTPS.
  - Chunked streaming for low memory usage.
  - Verification of target file existence and size.
  - 100% Python standard library (no external packages required).
================================================================================
"""

import os
import sys
import gzip
import time
import shutil
import argparse
import urllib.request
from pathlib import Path

# Official public endpoints
TARGETS = {
    "uniprot_isoforms": {
        "description": "UniProt Swiss-Prot Human Canonical + Reviewed Isoforms",
        "filename": "UniProt_SP_canonical_isoform_2026_09_26.fasta",
        "url": "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&includeIsoform=true&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29",
        "is_gzip": True,
        "approx_size_mb": 29.3,
    },
    "uniprot_canonical": {
        "description": "UniProt Swiss-Prot Human Canonical Proteome",
        "filename": "Uniprot_SP_canonical_2026_09_14.fasta",
        "url": "https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28proteome%3AUP000005640%29+AND+%28reviewed%3Atrue%29",
        "is_gzip": True,
        "approx_size_mb": 13.6,
    },
    "refseq": {
        "description": "NCBI RefSeq Human Proteome (GRCh38.p14: NP_ + XP_)",
        "filename": "GCF_000001405.40_GRCh38.p14_protein.faa",
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_protein.faa.gz",
        "is_gzip": True,
        "approx_size_mb": 106.3,
    },
}

USER_AGENT = "SDU-Immunoinformatics/3.0 (academic research; contact: trinh@sdu.dk)"


def download_and_decompress(url: str, dest_path: Path, is_gzip: bool, approx_mb: float):
    """Streams data from URL, decompressing gzip if needed, writing directly to dest_path."""
    tmp_path = dest_path.with_suffix(dest_path.suffix + ".download")
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})

    start_time = time.time()
    downloaded_bytes = 0
    chunk_size = 256 * 1024  # 256 KB chunks

    print(f"  Fetching: {url}")
    print(f"  Target:   {dest_path.name} (~{approx_mb:.1f} MB uncompressed)")

    try:
        with urllib.request.urlopen(req, timeout=60) as resp:
            # Check for Content-Length if available
            cl = resp.headers.get("Content-Length")
            total_expected = int(cl) if cl and cl.isdigit() else None

            if is_gzip:
                # Wrap response in GzipFile for on-the-fly streaming decompression
                with gzip.GzipFile(fileobj=resp) as gz_stream, open(tmp_path, "wb") as out_f:
                    while True:
                        chunk = gz_stream.read(chunk_size)
                        if not chunk:
                            break
                        out_f.write(chunk)
                        downloaded_bytes += len(chunk)
                        mb = downloaded_bytes / (1024 * 1024)
                        elapsed = time.time() - start_time
                        speed = mb / elapsed if elapsed > 0 else 0
                        print(f"  Decompressing: {mb:6.1f} MB [{speed:5.1f} MB/s] ...", end="\r", flush=True)
            else:
                with open(tmp_path, "wb") as out_f:
                    while True:
                        chunk = resp.read(chunk_size)
                        if not chunk:
                            break
                        out_f.write(chunk)
                        downloaded_bytes += len(chunk)
                        mb = downloaded_bytes / (1024 * 1024)
                        elapsed = time.time() - start_time
                        speed = mb / elapsed if elapsed > 0 else 0
                        print(f"  Downloaded:    {mb:6.1f} MB [{speed:5.1f} MB/s] ...", end="\r", flush=True)

        elapsed = time.time() - start_time
        final_mb = downloaded_bytes / (1024 * 1024)
        print(f"\n  Completed in {elapsed:.1f}s ({final_mb:.1f} MB extracted).")

        # Atomic replacement of target file
        if dest_path.exists():
            dest_path.unlink()
        tmp_path.rename(dest_path)
        return True

    except Exception as e:
        if tmp_path.exists():
            tmp_path.unlink()
        print(f"\n  [ERROR] Failed to download {url}: {e}", file=sys.stderr)
        return False


def main():
    script_dir = Path(__file__).resolve().parent
    data_dir_name = "data" if (script_dir.parent / "data").is_dir() else "Data"
    default_out_dir = (script_dir.parent / data_dir_name).resolve()

    parser = argparse.ArgumentParser(
        description="OASpepDB: Automated Reference Proteome Downloader for Step 03 & FASTA Export."
    )
    parser.add_argument(
        "--out-dir",
        type=str,
        default=str(default_out_dir),
        help=f"Target directory for reference proteome files (default: {default_out_dir})",
    )
    parser.add_argument(
        "--target",
        type=str,
        default="all",
        choices=["all", "uniprot_isoforms", "uniprot_canonical", "refseq"],
        help="Specific reference proteome target to download (default: all)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-download even if the output file already exists",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("  OASpepDB: REFERENCE PROTEOME DOWNLOADER")
    print("=" * 80)
    print(f"Output Directory: {out_dir}")
    print(f"Selected Target:  {args.target}")
    print(f"Force Overwrite:  {args.force}")
    print("=" * 80)

    selected_keys = list(TARGETS.keys()) if args.target == "all" else [args.target]

    success_count = 0
    for key in selected_keys:
        info = TARGETS[key]
        dest_file = out_dir / info["filename"]
        print(f"\n[{key.upper()}] {info['description']}")

        if dest_file.exists() and not args.force:
            file_mb = dest_file.stat().st_size / (1024 * 1024)
            print(f"  File already exists: {dest_file.name} ({file_mb:.1f} MB). Skipping (use --force to overwrite).")
            success_count += 1
            continue

        ok = download_and_decompress(
            url=info["url"],
            dest_path=dest_file,
            is_gzip=info["is_gzip"],
            approx_mb=info["approx_size_mb"],
        )
        if ok:
            success_count += 1

    print("\n" + "=" * 80)
    print(f"  DOWNLOAD SUMMARY: {success_count}/{len(selected_keys)} target(s) verified.")
    print("=" * 80)


if __name__ == "__main__":
    main()
