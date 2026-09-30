#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
OASpepDB: Step 01 - Automated OAS Raw Repertoire Downloader & Trimming Engine
================================================================================
Parallel downloader for human antibody repertoires from the Observed Antibody Space (OAS).
Applies stream-based Biological Quality Control (productive, no stop codon, in-frame VJ)
and retains minimal essential immunoproteomics fields to reduce storage by ~70%.
"""

import os
import sys
import io
import gzip
import json
import csv
import time
import argparse
import urllib.request
import urllib.error
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

# Ensure UTF-8 console output
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# Essential immunoproteomics fields preserved after Biological QC
TRIMMED_COLUMNS = [
    "cdr3_aa", "fr3_tail", "fwr4_aa", "v_call", "d_call", "j_call",
    "v_identity", "sequence_alignment_aa", "Redundancy"
]

csv_lock = Lock()


def format_time(seconds: float) -> str:
    """Formats duration in seconds into HH:MM:SS or MM:SS."""
    if seconds < 0 or seconds > 86400 * 30:
        return "--:--"
    m, s = divmod(int(seconds), 60)
    h, m = divmod(m, 60)
    return f"{h:02d}h {m:02d}m {s:02d}s" if h > 0 else f"{m:02d}m {s:02d}s"


def is_file_complete(filepath: str) -> bool:
    """Checks whether the trimmed file exists and contains valid filtered header."""
    if not os.path.exists(filepath) or os.path.getsize(filepath) < 50:
        return False
    try:
        with gzip.open(filepath, "rt", encoding="utf-8", errors="ignore") as f:
            return f.readline().strip().startswith("cdr3_aa,")
    except Exception:
        return False


def download_file(url: str, temp_path: str, max_retries: int = 3) -> bool:
    """Downloads a raw repertoire archive from OAS with HTTP Range resume capability."""
    headers = {"User-Agent": "SDU-Immunoinformatics/3.0 (academic research)"}
    for attempt in range(1, max_retries + 1):
        try:
            existing_size = os.path.getsize(temp_path) if os.path.exists(temp_path) else 0
            req_headers = dict(headers)
            mode = "wb"

            if existing_size > 0:
                req_headers["Range"] = f"bytes={existing_size}-"
                mode = "ab"

            req = urllib.request.Request(url, headers=req_headers)
            with urllib.request.urlopen(req, timeout=60) as resp:
                if existing_size > 0 and resp.status == 200:
                    mode = "wb"

                with open(temp_path, mode) as f_out:
                    while chunk := resp.read(256 * 1024):
                        f_out.write(chunk)

            if os.path.exists(temp_path) and os.path.getsize(temp_path) > 50:
                return True

        except urllib.error.HTTPError as e:
            if e.code == 416:  # Requested range already fully downloaded
                return True
            time.sleep(attempt)
        except Exception:
            time.sleep(attempt * 2)

    return False


def trim_file(raw_path: str, final_path: str, filename: str) -> tuple:
    """
    Parses a raw OAS CSV.gz file, applies Biological Quality Control,
    and streams the trimmed records directly to the destination archive.
    """
    raw_size = os.path.getsize(raw_path)
    temp_trimmed = final_path + ".tmp.gz"

    try:
        with gzip.open(raw_path, "rt", encoding="utf-8", errors="ignore") as f_in:
            line1 = f_in.readline().strip()
            if not line1:
                return False, raw_size, 0

            # Line 2: Column header
            line2 = f_in.readline().strip()
            cols = next(csv.reader([line2]))
            col_idx = {c: i for i, c in enumerate(cols)}

            # Mandatory biological QC columns
            req_cols = ["productive", "stop_codon", "vj_in_frame", "fwr3_aa", "cdr3_aa", "fwr4_aa"]
            for c in req_cols:
                if c not in col_idx:
                    return False, raw_size, 0

            idx_p = col_idx["productive"]
            idx_sc = col_idx["stop_codon"]
            idx_vj = col_idx["vj_in_frame"]
            idx_f3 = col_idx["fwr3_aa"]
            idx_c3 = col_idx["cdr3_aa"]
            idx_f4 = col_idx["fwr4_aa"]
            min_cols = max(idx_p, idx_sc, idx_vj, idx_f3, idx_c3, idx_f4) + 1

            idx_vc = col_idx.get("v_call", -1)
            idx_dc = col_idx.get("d_call", -1)
            idx_jc = col_idx.get("j_call", -1)
            idx_vi = col_idx.get("v_identity", -1)
            idx_seq = col_idx.get("sequence_alignment_aa", -1)
            idx_red = col_idx.get("Redundancy", -1)

            clean_name = filename[:-3] if filename.endswith(".gz") else filename
            with open(temp_trimmed, "wb") as raw_out:
                with gzip.GzipFile(filename=clean_name, mode="wb", fileobj=raw_out) as gz_out:
                    tw = io.TextIOWrapper(gz_out, encoding="utf-8", newline="", write_through=False)
                    writer = csv.writer(tw)
                    writer.writerow(TRIMMED_COLUMNS)

                    for row in csv.reader(f_in):
                        if len(row) < min_cols:
                            continue

                        # Biological QC filters
                        if row[idx_p] != "T" or row[idx_sc] != "F" or row[idx_vj] != "T":
                            continue

                        fwr3 = row[idx_f3].strip()
                        cdr3 = row[idx_c3].strip()
                        fwr4 = row[idx_f4].strip()

                        if not fwr3 or not cdr3 or not fwr4:
                            continue

                        # Minimum length checks
                        if len(fwr3) < 2 or len(fwr4) < 2 or len(cdr3) < 3:
                            continue

                        fr3_tail = fwr3[-2:]
                        v_call = row[idx_vc].strip() if idx_vc != -1 and len(row) > idx_vc else ""
                        d_call = row[idx_dc].strip() if idx_dc != -1 and len(row) > idx_dc else ""
                        j_call = row[idx_jc].strip() if idx_jc != -1 and len(row) > idx_jc else ""
                        v_id_val = row[idx_vi].strip() if idx_vi != -1 and len(row) > idx_vi else ""
                        seq_aa = row[idx_seq].strip() if idx_seq != -1 and len(row) > idx_seq else ""
                        redundancy = row[idx_red].strip() if idx_red != -1 and len(row) > idx_red else "1"

                        writer.writerow([
                            cdr3, fr3_tail, fwr4, v_call, d_call, j_call,
                            v_id_val, seq_aa, redundancy
                        ])

                    tw.flush()

        trimmed_size = os.path.getsize(temp_trimmed)
        os.replace(temp_trimmed, final_path)
        return True, raw_size, trimmed_size

    except Exception:
        if os.path.exists(temp_trimmed):
            try:
                os.remove(temp_trimmed)
            except Exception:
                pass
        return False, raw_size, 0


def process_url(entry: tuple, data_dir: str) -> tuple:
    """Full lifecycle for one repertoire URL: Download -> Trim & Filter -> Remove raw temp."""
    filename, url = entry
    final_path = os.path.join(data_dir, filename)

    if is_file_complete(final_path):
        return filename, "SKIP", 0, 0

    temp_raw = final_path + ".download.tmp"

    try:
        ok = download_file(url, temp_raw)
        if not ok:
            return filename, "DOWNLOAD_FAIL", 0, 0

        success, raw_sz, trim_sz = trim_file(temp_raw, final_path, filename)
        if not success:
            return filename, "TRIM_FAIL", raw_sz, 0

        return filename, "OK", raw_sz, trim_sz

    finally:
        if os.path.exists(temp_raw):
            try:
                os.remove(temp_raw)
            except Exception:
                pass


def get_urls_from_script(sh_path: str) -> list:
    """Parses download target URLs from the shell download script."""
    urls = []
    if not os.path.exists(sh_path):
        return urls
    with open(sh_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            parts = line.strip().split()
            for p in parts:
                if p.startswith("http://") or p.startswith("https://"):
                    filename = p.split("/")[-1].strip()
                    urls.append((filename, p))
                    break
    return urls


def load_completed_catalog(meta_csv: str, data_dir: str) -> set:
    """
    Reads the metadata catalog to identify already processed files and avoid duplicate work.
    Never destructively overwrites or alters the catalog file.
    """
    done_set = set()
    if not os.path.exists(meta_csv):
        return done_set

    try:
        with open(meta_csv, "r", encoding="utf-8", errors="ignore") as f:
            reader = csv.reader(f)
            header = next(reader, None)
            fn_idx = header.index("Filename") if header and "Filename" in header else -1
            for row in reader:
                if row and fn_idx != -1 and fn_idx < len(row):
                    fn = row[fn_idx].strip()
                    if is_file_complete(os.path.join(data_dir, fn)):
                        done_set.add(fn)
    except Exception as e:
        print(f"Notice: reading existing metadata {meta_csv} ({e})")

    return done_set


def main():
    script_dir = Path(__file__).resolve().parent
    default_meta = str((script_dir.parent / "Data" / "OAS_metadata.csv").resolve())
    default_sh = str((script_dir.parent / "Data" / "bulk_download_human_unpaired.sh").resolve())
    default_out = r"D:\OAS\human_unpaired"

    parser = argparse.ArgumentParser(
        description="OASpepDB Step 01: Automated OAS Raw Repertoire Downloader & Biological Trimmer"
    )
    # Support both --metadata and --meta-csv for seamless CLI compatibility
    parser.add_argument(
        "--metadata", "--meta-csv",
        dest="metadata",
        default=default_meta,
        help="Path to OAS_metadata.csv reference catalog (default: OASpepDB/Data/OAS_metadata.csv)"
    )
    # Support both --out-dir and --data-dir for seamless CLI compatibility
    parser.add_argument(
        "--out-dir", "--data-dir",
        dest="out_dir",
        default=default_out,
        help="Directory to save downloaded trimmed .csv.gz files (default: D:/OAS/human_unpaired)"
    )
    parser.add_argument(
        "--sh",
        default=default_sh,
        help="Path to shell download script containing repertoire URLs"
    )
    parser.add_argument(
        "--threads", "--workers",
        dest="threads",
        type=int,
        default=8,
        help="Number of concurrent download worker threads (default: 8)"
    )
    parser.add_argument(
        "--fresh",
        action="store_true",
        help="Clear existing trimmed files and restart from scratch"
    )
    args = parser.parse_args()

    data_dir = args.out_dir
    meta_csv = args.metadata
    sh_path = args.sh if os.path.exists(args.sh) else default_sh

    if args.fresh:
        print("[!] Fresh run requested: purging downloaded files in output directory...")
        if os.path.exists(data_dir):
            for f in os.listdir(data_dir):
                try:
                    os.remove(os.path.join(data_dir, f))
                except Exception:
                    pass
        print("    Directory cleaned.\n")

    os.makedirs(data_dir, exist_ok=True)

    # Identify already completed files safely (read-only inspection of metadata)
    completed_files = load_completed_catalog(meta_csv, data_dir)
    all_targets = get_urls_from_script(sh_path)
    total = len(all_targets)

    if total == 0:
        print(f"Error: No URLs found in script {sh_path}")
        return

    # Filter pending targets
    pending = [item for item in all_targets if not is_file_complete(os.path.join(data_dir, item[0]))]
    already_done = total - len(pending)

    print("=" * 80)
    print("  OASpepDB: STEP 01 - AUTOMATED REPERTOIRE DOWNLOAD & BIOLOGICAL TRIMMING")
    print("=" * 80)
    print(f"  Target Download Script: {sh_path}")
    print(f"  Metadata Catalog:       {meta_csv}")
    print(f"  Output Directory:       {data_dir}")
    print(f"  Total Repertoires:      {total:,}")
    print(f"  Already Processed:      {already_done:,}")
    print(f"  Pending to Process:     {len(pending):,}")
    print(f"  Worker Threads:         {args.threads}")
    print("=" * 80)

    if not pending:
        print("\nAll repertoire files are already downloaded, trimmed, and verified complete!")
        return

    t_start = time.time()
    done_count = 0
    total_raw = 0
    total_trimmed = 0

    executor = ThreadPoolExecutor(max_workers=args.threads)
    try:
        futures = {executor.submit(process_url, item, data_dir): item for item in pending}
        for idx, future in enumerate(as_completed(futures), 1):
            fn, status, r_sz, t_sz = future.result()
            elapsed = time.time() - t_start
            speed = idx / elapsed if elapsed > 0 else 0
            eta = format_time((len(pending) - idx) / speed) if speed > 0 else "--:--"

            if status == "OK":
                done_count += 1
                total_raw += r_sz
                total_trimmed += t_sz
                current_total = already_done + done_count
                pct = (current_total / total) * 100
                r_mb = r_sz / (1024 * 1024)
                t_mb = t_sz / (1024 * 1024)
                print(f"[{current_total:5d}/{total:,}] ({pct:5.1f}%) {fn[:32]:32s} | {r_mb:5.1f}MB -> {t_mb:4.1f}MB | ETA: {eta}")
            elif status != "SKIP":
                print(f"Warning: processing failure for {fn}: {status}")

    except KeyboardInterrupt:
        print("\nProcess interrupted by user. Shutting down gracefully...")
        executor.shutdown(wait=False, cancel_futures=True)
        sys.exit(0)
    finally:
        executor.shutdown(wait=True)

    print("-" * 80)
    print(f"Download & Trimming completed in: {format_time(time.time() - t_start)}")
    print(f"New Repertoire Files Processed:  {done_count:,}")
    print(f"Estimated Raw Size:              {total_raw / (1024**3):.2f} GB")
    print(f"Compressed Trimmed Size:         {total_trimmed / (1024**3):.2f} GB (approx. ~70% reduction)")
    print("=" * 80)


if __name__ == "__main__":
    main()
