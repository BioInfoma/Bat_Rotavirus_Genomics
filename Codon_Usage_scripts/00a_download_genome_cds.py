#!/usr/bin/env python3
"""
00a_download_genome_cds.py
Download CDS files from NCBI Assembly FTP for 7 bat species
(3 direct genome + 4 proxy genomes), parse CDS, count codons,
compute RSCU, and save per-species codon tables.

Usage:
    python3 00a_download_genome_cds.py

Output:
    c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables/{species}_codon.csv
"""

import urllib.request
import gzip
import os
import sys
import re
from collections import Counter

# ---- Genetic code ----
CT = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','TGA':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R','AGA':'R','AGG':'R',
    'AGT':'S','AGC':'S','GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}
STOPS = {'TAA','TAG','TGA'}
FAM = {
    'Phe':['TTT','TTC'],'Leu':['TTA','TTG','CTT','CTC','CTA','CTG'],
    'Ile':['ATT','ATC','ATA'],'Val':['GTT','GTC','GTA','GTG'],
    'Ser':['TCT','TCC','TCA','TCG','AGT','AGC'],'Pro':['CCT','CCC','CCA','CCG'],
    'Thr':['ACT','ACC','ACA','ACG'],'Ala':['GCT','GCC','GCA','GCG'],
    'Tyr':['TAT','TAC'],'His':['CAT','CAC'],'Gln':['CAA','CAG'],
    'Asn':['AAT','AAC'],'Lys':['AAA','AAG'],'Asp':['GAT','GAC'],'Glu':['GAA','GAG'],
    'Cys':['TGT','TGC'],'Arg':['CGT','CGC','CGA','CGG','AGA','AGG'],'Gly':['GGT','GGC','GGA','GGG']
}

# ---- Species to download ----
SPECIES = [
    # (output_name, assembly_acc, assembly_name)
    ("Rhinolophus_hipposideros", "GCF_964194185.1", "mRhiHip2.hap1.1"),
    ("Rousettus_aegyptiacus",    "GCF_014176215.1", "mRouAeg1.p"),
    ("Molossus_molossus",        "GCF_014108415.1", "mMolMol1.p"),
    # Proxy species
    ("Rhinolophus_ferrumequinum", "GCF_004115265.2", "mRhiFer1_v1.p"),
    ("Hipposideros_armiger",     "GCF_001890085.2", "ASM189008v1"),
    ("Saccopteryx_bilineata",    "GCF_036850765.1", "mSacBil1_pri_phased_curated"),
    ("Pipistrellus_kuhlii",      "GCF_014108245.1", "mPipKuh1.p"),
]

OUT_DIR = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables"
DL_DIR = "/workspace/cds_downloads"

def build_url(acc, name):
    prefix = acc[:3]  # GCF or GCA
    parts = acc[4:].split(".")
    p1, p2, p3 = parts[0][:3], parts[0][3:6], parts[0][6:]
    folder = f"{acc}_{name}"
    return f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{p1}/{p2}/{p3}/{folder}/{folder}_cds_from_genomic.fna.gz"

def download_cds(acc, name, outpath):
    url = build_url(acc, name)
    print(f"  Downloading: {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
        with open(outpath, 'wb') as f:
            f.write(data)
        print(f"  Saved: {outpath} ({len(data)/1e6:.1f} MB)")
        return True
    except Exception as e:
        print(f"  ERROR: {e}")
        return False

def parse_cds_fasta(gz_path):
    """Parse CDS FASTA, return list of sequences (uppercase, no gaps)"""
    seqs = []
    current_id = None
    current_seq = []
    with gzip.open(gz_path, 'rt') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    seqs.append(''.join(current_seq))
                current_id = line
                current_seq = []
            else:
                current_seq.append(line.upper())
    if current_seq:
        seqs.append(''.join(current_seq))
    return seqs

def count_codons(seq):
    """Count codons in a CDS sequence"""
    seq = seq.upper().replace('N', '')
    n = len(seq)
    if n % 3 != 0:
        return Counter()
    codons = [seq[i:i+3] for i in range(0, n-2, 3)]
    # Remove stop codon at end if present
    if codons and codons[-1] in STOPS:
        codons = codons[:-1]
    # Filter out codons with N or invalid chars
    codons = [c for c in codons if all(ch in 'ATGC' for ch in c)]
    return Counter(codons)

def compute_rscu(total_counts):
    """Compute RSCU from aggregated codon counts"""
    rscu = {}
    for aa, codons in FAM.items():
        tot = sum(total_counts.get(c, 0) for c in codons)
        if tot > 0:
            for c in codons:
                rscu[c] = total_counts.get(c, 0) * len(codons) / tot
        else:
            for c in codons:
                rscu[c] = 0.0
    return rscu

def process_species(name, acc, asm_name):
    print(f"\n{'='*60}")
    print(f"Processing: {name} ({acc})")
    print(f"{'='*60}")

    gz_path = os.path.join(DL_DIR, f"{acc}_cds_from_genomic.fna.gz")

    # Download if not already present
    if not os.path.exists(gz_path) or os.path.getsize(gz_path) < 1000:
        if not download_cds(acc, asm_name, gz_path):
            print(f"  FAILED to download {name}")
            return False

    # Parse CDS
    print(f"  Parsing CDS FASTA...")
    seqs = parse_cds_fasta(gz_path)
    print(f"  Found {len(seqs)} CDS sequences")

    # Count codons across all CDS
    total_counts = Counter()
    n_valid = 0
    for seq in seqs:
        if len(seq) >= 100 and len(seq) % 3 == 0:  # Skip very short or non-multiple-of-3
            cnt = count_codons(seq)
            if cnt:
                total_counts.update(cnt)
                n_valid += 1

    print(f"  {n_valid} valid CDS, {sum(total_counts.values())} total codons")

    # Compute RSCU
    rscu = compute_rscu(total_counts)

    # Save codon table
    out_path = os.path.join(OUT_DIR, f"{name}_codon.csv")
    with open(out_path, 'w') as f:
        f.write("codon,aa,count,rscu\n")
        for codon in sorted(CT.keys()):
            aa = CT[codon]
            cnt = total_counts.get(codon, 0)
            r = rscu.get(codon, 0.0)
            f.write(f"{codon},{aa},{cnt},{r:.6f}\n")

    print(f"  Saved codon table: {out_path}")
    return True

def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(DL_DIR, exist_ok=True)

    results = []
    for name, acc, asm_name in SPECIES:
        success = process_species(name, acc, asm_name)
        results.append((name, acc, success))

    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for name, acc, success in results:
        status = "OK" if success else "FAILED"
        print(f"  {name:40s} {acc:20s} {status}")

    failed = [r for r in results if not r[2]]
    if failed:
        print(f"\n{len(failed)} species failed!")
        sys.exit(1)
    else:
        print(f"\nAll {len(results)} species processed successfully!")

if __name__ == "__main__":
    main()
