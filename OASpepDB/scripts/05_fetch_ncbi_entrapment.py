#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
================================================================================
OASpepDB: Step 05 - NCBI Entrez Camelid VHH Entrapment Library Generator
(05_fetch_ncbi_entrapment.py)
================================================================================
Scientific Rationale:
  In bottom-up LC-MS/MS immunoproteomics, evaluating false discovery rates (FDR)
  using standard reversed decoy databases can underestimate error rates when
  searching hypervariable antibody repertoires. OASpepDB supports an empirical
  target-decoy entrapment benchmark strategy using authentic Camelid VHH
  (single-domain heavy-chain antibody) sequences:
    1. Camelid VHH repertoires possess antibody-like amino acid composition and
       hypervariable CDR3 loops, mirroring human antibody mass spectra properties.
    2. Because they derive from non-human Camelidae species (e.g. Vicugna pacos,
       Camelus dromedarius), authentic identifications in human patient biofluids
       represent definitive false discoveries.
    3. Isobaric normalization (I -> L) ensures uniform search space compatibility.
    4. Formatted as minimal-flank micro-cassettes with constant tail mimicry (+ASTK)
       under standard PE=4 (Predicted / Non-human entrapment) SV=1 flags.

Usage:
  python OASpepDB/Scripts/05_fetch_ncbi_entrapment.py --target-count 2000
  python OASpepDB/Scripts/05_fetch_ncbi_entrapment.py --target-count 20000 --out-file custom_entrapment.fasta
================================================================================
"""

import os
import re
import sys
import json
import time
import argparse
import urllib.request
import urllib.parse
from pathlib import Path

USER_AGENT = "SDU-Immunoinformatics/3.0 (academic research; contact: trinh@sdu.dk)"


def main():
    script_dir = Path(__file__).resolve().parent
    data_dir_name = "data" if (script_dir.parent / "data").is_dir() else "Data"
    default_out = (script_dir.parent / data_dir_name / "entrapment_cassettes.fasta").resolve()

    parser = argparse.ArgumentParser(
        description="OASpepDB: Fetch and format authentic Camelid VHH entrapment controls from NCBI Entrez."
    )
    parser.add_argument(
        "--target-count",
        type=int,
        default=2000,
        help="Target number of unique Camelid VHH micro-cassettes to extract (default: 2000)",
    )
    parser.add_argument(
        "--out-file",
        type=str,
        default=str(default_out),
        help=f"Output FASTA file path (default: {default_out})",
    )
    parser.add_argument(
        "--retmax",
        type=int,
        default=25000,
        help="Maximum candidate records to query via NCBI ESearch (default: 25000)",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=200,
        help="Number of records per EFetch request (default: 200)",
    )
    parser.add_argument(
        "--delay",
        type=float,
        default=0.35,
        help="Delay in seconds between NCBI API requests to respect polite rate limits (default: 0.35)",
    )
    args = parser.parse_args()

    target_count = args.target_count
    out_path = Path(args.out_file)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("  OASpepDB: NCBI CAMELID VHH ENTRAPMENT LIBRARY GENERATOR")
    print("=" * 80)
    print(f"Target Cassettes:  {target_count:,}")
    print(f"Output FASTA:      {out_path}")
    print(f"ESearch RetMax:    {args.retmax:,}")
    print(f"Batch Fetch Size:  {args.batch_size}")
    print("=" * 80)

    # 1. Query NCBI Entrez Protein
    print("\n[STEP 1/4] Querying NCBI Entrez E-utilities for published Camelid VHH records...")
    term = 'Camelidae[Organism] AND (immunoglobulin OR antibody OR VHH OR nanobody OR "single domain")'
    esearch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=protein&term={urllib.parse.quote(term)}&retmax={args.retmax}&retmode=json"
    )

    req = urllib.request.Request(esearch_url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as e:
        print(f"[FATAL] NCBI ESearch failed: {e}", file=sys.stderr)
        sys.exit(1)

    id_list = data.get("esearchresult", {}).get("idlist", [])
    print(f"  Found {len(id_list):,} candidate Camelid antibody protein records in NCBI.")

    if not id_list:
        print("[FATAL] No records returned from NCBI query.", file=sys.stderr)
        sys.exit(1)

    # 2. Fetch records in batches and parse micro-cassettes
    print(f"\n[STEP 2/4] Streaming FASTA records and assembling minimal-flank micro-cassettes...")
    valid_entries = []
    seen_cassettes = set()
    batch_size = args.batch_size

    for i in range(0, len(id_list), batch_size):
        if len(valid_entries) >= target_count:
            break

        batch_ids = id_list[i:i + batch_size]
        efetch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=protein&id={','.join(batch_ids)}&rettype=fasta&retmode=text"
        )

        fasta_text = ""
        for attempt in range(3):
            try:
                req = urllib.request.Request(efetch_url, headers={"User-Agent": USER_AGENT})
                with urllib.request.urlopen(req, timeout=30) as resp:
                    fasta_text = resp.read().decode("utf-8", errors="ignore")
                break
            except Exception as e:
                time.sleep(1.0 * (attempt + 1))
                if attempt == 2:
                    print(f"  Warning: Batch {i//batch_size + 1} fetch failed ({e}), skipping batch.")

        if not fasta_text:
            continue

        blocks = fasta_text.strip().split(">")[1:]
        for block in blocks:
            if len(valid_entries) >= target_count:
                break

            lines = block.strip().split("\n")
            header = lines[0]
            seq = "".join(lines[1:]).replace(" ", "").upper()

            parts = header.split(None, 1)
            acc = parts[0]
            desc = parts[1] if len(parts) > 1 else ""

            org_match = re.search(r"\[([^\]]+)\]", desc)
            organism = org_match.group(1).replace(" ", "_") if org_match else "Camelidae"

            # Locate FR3 C-terminus: [YFHW][YFALITMS]C (IMGT Cys104 anchor)
            m_c = list(re.finditer(r"(?:[YFHW][YFALITMS]C)[A-Z]", seq))
            # Locate FR4: W[GQRK]...V[TS][VLA]SS (IMGT Trp118 anchor to J-terminus)
            m_w = list(re.finditer(r"W[GQRK][A-Z]{2,6}V[TS][VLA]S[SA]", seq))

            if m_c and m_w:
                c_pos = m_c[-1].start() + 2  # 0-indexed position of Cys104
                w_pos = m_w[-1].start()      # 0-indexed position of Trp118

                # CDR3 hypervariable loop between Cys104 and Trp118
                if w_pos > c_pos + 4:
                    cdr3 = seq[c_pos + 1:w_pos]
                    fr4 = m_w[-1].group()
                    # 11 aa minimal flank ending at Cys104
                    flank = seq[max(0, c_pos - 10):c_pos + 1]

                    # Assemble micro-cassette with +ASTK isotype extension tail
                    cassette = (flank + cdr3 + fr4 + "ASTK").replace("I", "L")

                    if 25 <= len(cassette) <= 80 and cassette not in seen_cassettes:
                        seen_cassettes.add(cassette)
                        valid_entries.append({
                            "acc": acc,
                            "organism": organism,
                            "cassette": cassette,
                            "cdr3_len": len(cdr3),
                        })

        pct = (len(valid_entries) / target_count) * 100
        print(f"  Progress: {len(valid_entries):,}/{target_count:,} ({pct:5.1f}%) cassettes assembled...", end="\r", flush=True)
        time.sleep(args.delay)

    print(f"\n[STEP 3/4] Successfully assembled {len(valid_entries):,} authentic Camelid VHH cassettes.")

    # 4. Serialize to FASTA
    print(f"[STEP 4/4] Writing FASTA records to {out_path}...")
    with open(out_path, "w", encoding="utf-8") as f:
        for entry in valid_entries:
            f.write(f">ENTRAPMENT|NCBI_{entry['acc']} Species={entry['organism']} Locus=VHH PE=4 SV=1\n")
            f.write(f"{entry['cassette']}\n")

    print("\n" + "=" * 80)
    print("  CAMELID VHH ENTRAPMENT LIBRARY GENERATED SUCCESSFULLY!")
    print("=" * 80)
    print(f"  Output File:      {out_path}")
    print(f"  Total Cassettes:  {len(valid_entries):,}")
    print(f"  Isobaric State:   Isobaric I -> L Converted")
    print(f"  Annotation Flag:  PE=4 SV=1 (Empirical FDR Entrapment Benchmark)")
    print("=" * 80)


if __name__ == "__main__":
    main()
