#!/usr/bin/env python3
"""
================================================================================
OASpepDB: Step 05 - NCBI Entrez Camelid VHH Entrapment Library Generator
================================================================================
Fetch authentic published Camelid VHH (Nanobody) sequences from NCBI GenBank.
Converts them into minimal-flank micro-cassettes with Isobaric I->L conversion
for Empirical FDR benchmarking in OASpepDB (PE=4 SV=1 standard).
"""

from pathlib import Path
import urllib.request
import urllib.parse
import json
import re
import time
import os
import sys

def main():
    target_count = 2000
    script_dir = Path(__file__).resolve().parent
    out_path = str((script_dir.parent / "Data" / "entrapment_cassettes.fasta").resolve())
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    
    print("[1/4] Querying NCBI Entrez E-utilities for published Camelid VHH records...")
    term = 'Camelidae[Organism] AND (immunoglobulin OR antibody OR VHH OR nanobody OR "single domain")'
    esearch_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
        f"db=protein&term={urllib.parse.quote(term)}&retmax=25000&retmode=json"
    )
    
    req = urllib.request.Request(esearch_url, headers={"User-Agent": "SDU-Immunoinformatics/3.0 (academic research)"})
    with urllib.request.urlopen(req) as resp:
        data = json.loads(resp.read().decode())
    
    id_list = data["esearchresult"]["idlist"]
    print(f"      Found {len(id_list)} published candidate records in NCBI Protein.")
    
    print("[2/4] Fetching FASTA records in batches and parsing micro-cassette architecture...")
    valid_entries = []
    seen_cassettes = set()
    batch_size = 200
    
    for i in range(0, len(id_list), batch_size):
        if len(valid_entries) >= target_count:
            break
            
        batch_ids = id_list[i:i + batch_size]
        efetch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=protein&id={','.join(batch_ids)}&rettype=fasta&retmode=text"
        )
        
        try:
            req = urllib.request.Request(efetch_url, headers={"User-Agent": "SDU-Immunoinformatics/3.0 (academic research)"})
            with urllib.request.urlopen(req) as resp:
                fasta_text = resp.read().decode("utf-8", errors="ignore")
        except Exception as e:
            print(f"      Warning: Batch {i//batch_size + 1} fetch failed ({e}), retrying after delay...")
            time.sleep(1)
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
            
            # Locate FR3 C-terminus: [YFHW][YFALITMS]C (IMGT Cys104)
            m_c = list(re.finditer(r"(?:[YFHW][YFALITMS]C)[A-Z]", seq))
            # Locate FR4: W[GQRK]...V[TS][VLA]SS (IMGT Trp118 to J-terminus)
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
                            "cdr3_len": len(cdr3)
                        })
                        
        print(f"      Progress: {len(valid_entries)}/{target_count} authentic Camelid VHH cassettes extracted...")
        time.sleep(0.35)  # Enforce polite NCBI rate limit
        
    print(f"[3/4] Successfully assembled {len(valid_entries)} authentic NCBI Camelid VHH cassettes.")
    
    print(f"[4/4] Writing to {out_path}...")
    with open(out_path, "w", encoding="utf-8") as f:
        for idx, entry in enumerate(valid_entries, 1):
            f.write(f">ENTRAPMENT|NCBI_{entry['acc']} Species={entry['organism']} Locus=VHH PE=4 SV=1\n")
            f.write(f"{entry['cassette']}\n")
            
    print(f"[DONE] Saved {len(valid_entries)} NCBI Entrapment micro-cassettes to {out_path}!")

if __name__ == "__main__":
    main()
