#!/usr/bin/env python3
"""
00d_build_transcriptome_table.py
================================
Process Trinity transcriptome assembly output through TransDecoder to extract
CDS sequences, then build a host codon-usage table (RSCU + raw counts).

This script is run AFTER the Trinity HPC job completes and results are downloaded.

Usage:
    python 00d_build_transcriptome_table.py <species_name> <trinity_fasta_path>

Example:
    python 00d_build_transcriptome_table.py Rhinolophus_simulator c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/pipeline-xxx/outputs/Trinity.fasta

Requirements:
    - TransDecoder installed at /workspace/TransDecoder/
    - Biopython for FASTA parsing
"""

import sys
import os
import subprocess
from collections import Counter
from pathlib import Path

# ---- Genetic code (standard table 1) ----
GENETIC_CODE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

STOP_CODONS = {'TAA', 'TAG', 'TGA'}
SYNONYMOUS_CODONS = {}  # aa -> list of codons
for codon, aa in GENETIC_CODE.items():
    if aa != '*':
        SYNONYMOUS_CODONS.setdefault(aa, []).append(codon)
# Sort codons within each AA for consistency
for aa in SYNONYMOUS_CODONS:
    SYNONYMOUS_CODONS[aa].sort()

# All 64 codon names (sorted)
ALL_CODONS = sorted(GENETIC_CODE.keys())


def run_transdecoder(trinity_fasta, work_dir):
    """Run TransDecoder LongOrfs + Predict on Trinity assembly."""
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    transdecoder_dir = "/workspace/TransDecoder"
    long_orfs = os.path.join(transdecoder_dir, "util", "TransDecoder.LongOrfs")
    predict = os.path.join(transdecoder_dir, "util", "TransDecoder.Predict")

    # Step 1: LongOrfs
    print(f"  Running TransDecoder.LongOrfs...")
    cmd1 = f"perl {long_orfs} -t {trinity_fasta} --output_dir {work_dir}/long_orfs"
    result = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  LongOrs stderr: {result.stderr[:500]}")
        # LongOrfs may return non-zero but still produce output

    # Step 2: Predict (must use same output_dir as LongOrfs)
    print(f"  Running TransDecoder.Predict...")
    cmd2 = f"perl {predict} -t {trinity_fasta} --output_dir {work_dir}/long_orfs"
    result = subprocess.run(cmd2, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Predict stderr: {result.stderr[:500]}")

    # Find CDS file
    cds_file = work_dir / "long_orfs" / f"{Path(trinity_fasta).name}.transdecoder.cds"
    if not cds_file.exists():
        # Try alternative locations
        candidates = list(work_dir.rglob("*.transdecoder.cds"))
        if candidates:
            cds_file = candidates[0]
        else:
            raise FileNotFoundError(f"TransDecoder CDS file not found in {work_dir}")

    print(f"  CDS file: {cds_file}")
    return str(cds_file)


def parse_fasta(fasta_path):
    """Parse FASTA file, yield (header, sequence) tuples."""
    header = None
    seq_parts = []
    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header is not None:
                    yield header, ''.join(seq_parts)
                header = line[1:]
                seq_parts = []
            else:
                seq_parts.append(line.upper())
    if header is not None:
        yield header, ''.join(seq_parts)


def count_codons_in_cds(cds_fasta):
    """Count all codons across all CDS sequences."""
    total_counts = Counter()
    n_seqs = 0
    n_valid = 0

    for header, seq in parse_fasta(cds_fasta):
        n_seqs += 1
        seq = seq.upper().replace('U', 'T')
        n = len(seq)
        # Trim to multiple of 3
        if n % 3 != 0:
            seq = seq[:n - (n % 3)]
            n = len(seq)
        if n < 3:
            continue

        # Check for internal stop codons — skip sequences with internal stops
        codons = [seq[i:i+3] for i in range(0, n, 3)]
        # Remove terminal stop codon
        if codons and codons[-1] in STOP_CODONS:
            codons = codons[:-1]

        # Skip if internal stop codons
        has_internal_stop = any(c in STOP_CODONS for c in codons)
        if has_internal_stop:
            continue

        # Count codons
        for c in codons:
            if c in GENETIC_CODE:
                total_counts[c] += 1
        n_valid += 1

    print(f"  Total CDS: {n_seqs}, valid (no internal stops): {n_valid}")
    print(f"  Total codons: {sum(total_counts.values()):,}")
    return total_counts


def compute_rscu(counts):
    """Compute RSCU from codon counts."""
    rscu = {}
    for aa, codons in SYNONYMOUS_CODONS.items():
        total = sum(counts.get(c, 0) for c in codons)
        n_syn = len(codons)
        if total == 0 or n_syn == 0:
            for c in codons:
                rscu[c] = 0.0
        else:
            for c in codons:
                rscu[c] = counts.get(c, 0) * n_syn / total
    return rscu


def save_codon_table(species, counts, rscu, output_path):
    """Save codon table as CSV with columns: codon, aa, count, rscu."""
    rows = []
    for codon in ALL_CODONS:
        aa = GENETIC_CODE[codon]
        rows.append({
            'codon': codon,
            'aa': aa,
            'count': counts.get(codon, 0),
            'rscu': rscu.get(codon, 0.0)
        })

    import csv
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['codon', 'aa', 'count', 'rscu'])
        writer.writeheader()
        writer.writerows(rows)
    print(f"  Saved: {output_path}")


def main():
    if len(sys.argv) < 3:
        print("Usage: python 00d_build_transcriptome_table.py <species_name> <trinity_fasta_path>")
        print("Example: python 00d_build_transcriptome_table.py Rhinolophus_simulator /path/to/Trinity.fasta")
        sys.exit(1)

    species = sys.argv[1]
    trinity_fasta = sys.argv[2]

    if not os.path.exists(trinity_fasta):
        print(f"ERROR: Trinity FASTA not found: {trinity_fasta}")
        sys.exit(1)

    output_dir = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables"
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{species}_codon.csv")

    work_dir = f"/workspace/transdecoder_{species}"

    print(f"=== Building codon table for {species} ===")
    print(f"  Trinity input: {trinity_fasta}")

    # Step 1: Run TransDecoder
    print("\n--- TransDecoder ---")
    cds_file = run_transdecoder(trinity_fasta, work_dir)

    # Step 2: Count codons
    print("\n--- Codon counting ---")
    counts = count_codons_in_cds(cds_file)

    # Step 3: Compute RSCU
    print("\n--- RSCU computation ---")
    rscu = compute_rscu(counts)

    # Step 4: Save
    print("\n--- Saving ---")
    save_codon_table(species, counts, rscu, output_path)

    # Also save to /workspace for R file.copy compatibility
    workspace_copy = f"/workspace/{species}_codon.csv"
    import shutil
    shutil.copy2(output_path, workspace_copy)

    print(f"\n=== Done: {species} ===")
    print(f"  Codon table: {output_path}")
    print(f"  Total codons: {sum(counts.values()):,}")


if __name__ == '__main__':
    main()
