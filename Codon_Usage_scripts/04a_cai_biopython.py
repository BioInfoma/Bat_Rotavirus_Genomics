#!/usr/bin/env python3
# ============================================================
# 04a_cai_biopython.py
# CAI (Sharp & Li 1987) of each viral sequence against its own
# host's codon usage table, computed with Biopython's
# Bio.SeqUtils.CodonAdaptationIndex — the standard, citable
# implementation (Cock et al. 2009, Bioinformatics).
#
# Reference handling: each host codon-usage table (counts) is
# expanded into a pseudo-sequence with exactly those codon
# multiplicities and passed to CodonAdaptationIndex(), which then
# computes the Sharp & Li w values (count/max per synonymous
# family, 0.5 for codons absent from the reference) entirely
# inside Biopython. Viral sequences are preprocessed exactly as
# in the rest of the pipeline (uppercase, gap-stripped, trimmed
# to the last complete codon, terminal stop removed) and scored
# with CodonAdaptationIndex.calculate(), which excludes ATG, TGG
# and stop codons per Sharp & Li.
#
# Inputs:
#   output/tables/{seg}_cds.fasta, {seg}_metadata_parsed.csv
#   host_codon_tables/*.csv + host_table_map.csv
# Output:
#   output/tables/cai_values_{seg}.csv  (accession, host, CAI, GC, full_name)
# ============================================================

import csv
import os
from Bio.SeqUtils import CodonAdaptationIndex

TAB = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/output_tables"
FASTA_DIR = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/data/fastas"
META_DIR = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/data/metadata"
CODON_DIR = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/host_codon_tables"
SEGMENTS = ["VP4", "VP6", "VP7", "NSP4"]
LEN_THRESH = {"VP4": 2300, "VP6": 1150, "VP7": 950, "NSP4": 500}
STOPS = {"TAA", "TAG", "TGA"}


def trim_cds(seq):
    """Uppercase, strip gaps, trim to last complete codon, drop terminal stop,
    and drop any codon containing ambiguous (non-ACGT) bases. Dropping whole
    codons preserves the reading frame and matches the rest of the pipeline,
    where ambiguous codons never enter codon counts."""
    seq = seq.upper().replace("-", "")
    n = len(seq) - (len(seq) % 3)
    seq = seq[:n]
    codons = [seq[i:i+3] for i in range(0, n, 3)]
    codons = [c for c in codons if set(c) <= set("ACGT")]
    if codons and codons[-1] in STOPS:
        codons = codons[:-1]
    return "".join(codons)


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line)
    if name:
        seqs[name] = "".join(buf)
    return seqs


def load_host_indices():
    """Build one Biopython CodonAdaptationIndex per host table file."""
    # host species -> table file
    mapping = {}
    with open(os.path.join(CODON_DIR, "host_table_map.csv")) as fh:
        for row in csv.DictReader(fh):
            mapping[row["host_species"]] = row["codon_table"]

    indices = {}
    for host, fname in mapping.items():
        path = os.path.join(CODON_DIR, fname)
        if not os.path.exists(path):
            print(f"  WARNING: missing table {fname} for {host}")
            continue
        if fname in indices:  # table shared by proxy hosts
            indices[host] = indices[fname]
            continue
        parts = []
        total = 0
        with open(path) as fh:
            for row in csv.DictReader(fh):
                c = int(row["count"])
                if c > 0:
                    parts.append(row["codon"].upper() * c)
                    total += c
        ref_seq = "".join(parts)
        idx = CodonAdaptationIndex([ref_seq])   # Sharp & Li w values
        indices[host] = idx
        indices[fname] = idx                    # cache by filename too
        print(f"  Built CAI reference for {host} from {fname} ({total:,} codons)")
    return indices


def main():
    host_idx = load_host_indices()
    for seg in SEGMENTS:
        seqs = read_fasta(os.path.join(FASTA_DIR, f"{seg}_cds.fasta"))
        meta = {}
        with open(os.path.join(META_DIR, f"{seg}_metadata_parsed.csv")) as fh:
            for row in csv.DictReader(fh):
                meta[row["accession"]] = row

        out_rows, n_skip = [], 0
        for sid, raw in seqs.items():
            seq = trim_cds(raw)
            if len(seq) < LEN_THRESH[seg]:
                continue
            m = meta.get(sid) or meta.get(sid.split("_NODE_")[0])
            if m is None:
                n_skip += 1
                continue
            host = m["host"]
            idx = host_idx.get(host)
            cai = idx.calculate(seq) if idx is not None else ""
            out_rows.append({"accession": sid, "host": host, "CAI": cai,
                             "GC": m["GC"], "full_name": m.get("full_name", "")})

        out_path = os.path.join(TAB, f"cai_values_{seg}.csv")
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["accession", "host", "CAI", "GC", "full_name"])
            w.writeheader()
            w.writerows(out_rows)
        vals = [r["CAI"] for r in out_rows if r["CAI"] != ""]
        print(f"{seg}: {len(out_rows)} sequences scored "
              f"(median CAI {sorted(vals)[len(vals)//2]:.4f}); "
              f"{n_skip} without metadata skipped")


if __name__ == "__main__":
    main()
