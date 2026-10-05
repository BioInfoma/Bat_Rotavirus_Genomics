#!/usr/bin/env python3
"""
00b_build_ehelvum_table.py
Extract CDS sequences from E. helvum figshare genome + GFF annotation,
count codons, compute RSCU, and save codon table.

Streams the GFF from the ZIP file (10.5 GB uncompressed) without extracting it.
Uses Biopython SeqIO.index for memory-efficient genome access.

Usage:
    python3 00b_build_ehelvum_table.py

Output:
    c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables/Eidolon_helvum_codon.csv
"""

import os
import sys
import gzip
import zipfile
import subprocess
from collections import Counter, defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

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

GENOME_PATH = "/workspace/ehelvum/Eidolon_helvum_draft_genome.fa"
ZIP_PATH = "/workspace/ehelvum/Eidolon_helvum_annotation_LMU.zip"
GFF_NAME = "Eidolon_helvum_QMUgenome_MAKER_annot.renamed.noHash_iprscanAnnot_blastpAnnot.gff"
OUT_PATH = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables/Eidolon_helvum_codon.csv"
CDS_CACHE = "/workspace/ehelvum/ehelvum_cds_coords.tsv"

def extract_cds_coords():
    """Stream GFF from ZIP, extract CDS coordinates grouped by transcript"""
    if os.path.exists(CDS_CACHE):
        print(f"  Using cached CDS coordinates: {CDS_CACHE}")
        return
    
    print("  Streaming GFF from ZIP and extracting CDS coordinates...")
    cds_by_transcript = defaultdict(list)
    
    with zipfile.ZipFile(ZIP_PATH, 'r') as zf:
        with zf.open(GFF_NAME) as f:
            for line in f:
                line = line.decode('utf-8')
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 9:
                    continue
                if parts[2] != 'CDS':
                    continue
                scaffold = parts[0]
                start = int(parts[3])
                end = int(parts[4])
                strand = parts[6]
                phase = int(parts[7]) if parts[7] != '.' else 0
                attrs = parts[8]
                
                # Extract Parent ID
                parent = None
                for attr in attrs.split(';'):
                    if attr.strip().startswith('Parent='):
                        parent = attr.strip().split('=')[1]
                        break
                if not parent:
                    continue
                
                cds_by_transcript[parent].append({
                    'scaffold': scaffold,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phase': phase,
                })
    
    print(f"  Found {len(cds_by_transcript)} transcripts with CDS")
    
    # Save to cache file
    with open(CDS_CACHE, 'w') as f:
        for transcript, exons in cds_by_transcript.items():
            for ex in exons:
                f.write(f"{transcript}\t{ex['scaffold']}\t{ex['start']}\t{ex['end']}\t{ex['strand']}\t{ex['phase']}\n")
    
    print(f"  Cached to {CDS_CACHE}")
    return cds_by_transcript

def load_cds_coords():
    """Load cached CDS coordinates"""
    cds_by_transcript = defaultdict(list)
    with open(CDS_CACHE, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            transcript = parts[0]
            cds_by_transcript[transcript].append({
                'scaffold': parts[1],
                'start': int(parts[2]),
                'end': int(parts[3]),
                'strand': parts[4],
                'phase': int(parts[5]),
            })
    return cds_by_transcript

def extract_cds_sequences(cds_by_transcript, genome_index):
    """Extract CDS sequences from genome, concatenate exons per transcript"""
    total_counts = Counter()
    n_transcripts = 0
    n_valid = 0
    
    for transcript_id, exons in cds_by_transcript.items():
        n_transcripts += 1
        if n_transcripts % 5000 == 0:
            print(f"  Processing transcript {n_transcripts}/{len(cds_by_transcript)}...")
        
        # Sort exons by position
        if exons[0]['strand'] == '+':
            exons_sorted = sorted(exons, key=lambda x: x['start'])
        else:
            exons_sorted = sorted(exons, key=lambda x: x['start'], reverse=True)
        
        # Extract and concatenate exon sequences
        cds_seq = ''
        for ex in exons_sorted:
            try:
                scaffold_seq = genome_index[ex['scaffold']]
                exon_seq = str(scaffold_seq[ex['start']-1:ex['end']].seq)
                if ex['strand'] == '-':
                    exon_seq = str(Seq(exon_seq).reverse_complement())
                cds_seq += exon_seq
            except (KeyError, IndexError):
                break
        
        if not cds_seq or len(cds_seq) < 100:
            continue
        if len(cds_seq) % 3 != 0:
            continue
        
        # Count codons
        codons = [cds_seq[i:i+3] for i in range(0, len(cds_seq)-2, 3)]
        if codons and codons[-1] in STOPS:
            codons = codons[:-1]
        codons = [c for c in codons if all(ch in 'ATGC' for ch in c)]
        
        if codons:
            total_counts.update(codons)
            n_valid += 1
    
    print(f"  {n_valid}/{n_transcripts} valid transcripts, {sum(total_counts.values())} total codons")
    return total_counts

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

def main():
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    
    # Step 1: Extract CDS coordinates from GFF
    if not os.path.exists(CDS_CACHE):
        extract_cds_coords()
    
    cds_by_transcript = load_cds_coords()
    print(f"  Loaded {len(cds_by_transcript)} transcripts")
    
    # Step 2: Index genome FASTA
    print("  Indexing genome FASTA (may take a few minutes)...")
    genome_index = SeqIO.index(GENOME_PATH, "fasta")
    print(f"  Genome indexed: {len(genome_index)} scaffolds")
    
    # Step 3: Extract CDS sequences and count codons
    print("  Extracting CDS sequences and counting codons...")
    total_counts = extract_cds_sequences(cds_by_transcript, genome_index)
    
    # Step 4: Compute RSCU
    rscu = compute_rscu(total_counts)
    
    # Step 5: Save codon table
    with open(OUT_PATH, 'w') as f:
        f.write("codon,aa,count,rscu\n")
        for codon in sorted(CT.keys()):
            aa = CT[codon]
            cnt = total_counts.get(codon, 0)
            r = rscu.get(codon, 0.0)
            f.write(f"{codon},{aa},{cnt},{r:.6f}\n")
    
    print(f"  Saved codon table: {OUT_PATH}")
    print(f"  Total codons: {sum(total_counts.values())}")

if __name__ == "__main__":
    main()
