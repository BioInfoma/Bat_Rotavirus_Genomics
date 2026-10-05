#!/usr/bin/env python3
"""
01a_parse_genbank.py
Parse GenBank files for VP6, VP7, NSP4 to extract CDS sequences and host info.
Also parse study sequences from four_segment.fasta.txt.
Save CDS as FASTA and metadata as CSV for downstream R analysis.

Usage:
    python3 01a_parse_genbank.py

Output:
    c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/output/tables/{segment}_cds.fasta  — CDS sequences per segment
    c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/output/tables/{segment}_metadata_parsed.csv  — parsed metadata per segment
"""

import os
import re
import csv
from Bio import SeqIO

# ---- Paths ----
INPUT_DIR = "/mnt/user-uploads/CU_bat_only"
STUDY_FASTA = "/mnt/user-uploads/four_segment.fasta.txt"
OUT_DIR = "c:/Users/USER/Documents/persornal_project/CU_bat_only/bat_rva_cub_reproducibility/output/tables"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- GenBank files ----
GB_FILES = {
    "VP6":  os.path.join(INPUT_DIR, "VP6/bat_VP6_records_2.gb"),
    "VP7":  os.path.join(INPUT_DIR, "VP7/VP7_records.gb"),
    "NSP4": os.path.join(INPUT_DIR, "NSP4/NSP4_records.gb"),
}

# ---- Expected CDS lengths ----
EXPECTED_LEN = {
    "VP4":  (2300, 2400),   # ~2349 bp, threshold ≥2300
    "VP6":  (1150, 1250),   # ~1194 bp
    "VP7":  (950, 1000),    # ~981 bp
    "NSP4": (500, 560),     # ~525 bp
}

def extract_host(record):
    """Extract host from GenBank source feature"""
    for feat in record.features:
        if feat.type == "source":
            host = feat.qualifiers.get("host", ["Unknown"])[0]
            # Clean up host name - remove common name in parentheses
            host = re.sub(r'\s*\(.*?\)\s*', '', host).strip()
            return host
    return "Unknown"

def extract_country(record):
    """Extract country from GenBank source feature"""
    for feat in record.features:
        if feat.type == "source":
            country = feat.qualifiers.get("geo_loc_name", ["Unknown"])[0]
            return country.split(":")[0]  # Remove region detail
    return "Unknown"

def extract_year(record):
    """Extract collection year from GenBank source feature"""
    for feat in record.features:
        if feat.type == "source":
            date = feat.qualifiers.get("collection_date", ["Unknown"])[0]
            # Extract year from date string
            year_match = re.search(r'(20\d{2}|19\d{2})', date)
            if year_match:
                return year_match.group(1)
    return "Unknown"

def extract_cds_sequence(record):
    """Extract CDS sequence from GenBank record"""
    for feat in record.features:
        if feat.type == "CDS":
            cds_seq = feat.extract(record.seq)
            return str(cds_seq)
    # If no CDS feature, return the full sequence (for study sequences)
    return str(record.seq)

def parse_genbank(segment, gb_path):
    """Parse a GenBank file and return list of (accession, seq, host, country, year, full_name)"""
    records = []
    for record in SeqIO.parse(gb_path, "genbank"):
        acc = f"{record.id}"
        cds = extract_cds_sequence(record)
        host = extract_host(record)
        country = extract_country(record)
        year = extract_year(record)
        full_name = record.description
        records.append({
            "accession": acc,
            "full_name": full_name,
            "host": host,
            "country": country,
            "year": year,
            "segment": segment,
            "sequence": cds,
            "length": len(cds),
            "source": "ref",
        })
    return records

def parse_study_fasta(fasta_path):
    """Parse the four_segment.fasta.txt file with section headers"""
    segments = {}
    current_segment = None
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # Check for section headers (e.g., "VP7:", "VP4:", "VP6:", "NSP4:", "NSP2:")
            if re.match(r'^(VP[467]|NSP[24]):', line):
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    segments.setdefault(current_segment, []).append({
                        "accession": current_id,
                        "full_name": current_id,
                        "host": "Eidolon helvum",
                        "country": "Nigeria",
                        "year": "Unknown",
                        "segment": current_segment,
                        "sequence": seq,
                        "length": len(seq),
                        "source": "Study",
                    })
                current_segment = line.split(":")[0]
                current_id = None
                current_seq = []
            elif line.startswith(">"):
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    segments.setdefault(current_segment, []).append({
                        "accession": current_id,
                        "full_name": current_id,
                        "host": "Eidolon helvum",
                        "country": "Nigeria",
                        "year": "Unknown",
                        "segment": current_segment,
                        "sequence": seq,
                        "length": len(seq),
                        "source": "Study",
                    })
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    
    # Don't forget the last sequence
    if current_id and current_seq:
        seq = ''.join(current_seq)
        segments.setdefault(current_segment, []).append({
            "accession": current_id,
            "full_name": current_id,
            "host": "Eidolon helvum",
            "country": "Nigeria",
            "year": "Unknown",
            "segment": current_segment,
            "sequence": seq,
            "length": len(seq),
            "source": "Study",
        })
    
    return segments

def parse_vp4_fasta(fasta_path):
    """Parse the VP4 matched FASTA file"""
    records = []
    current_id = None
    current_seq = []
    
    with open(fasta_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id and current_seq:
                    seq = ''.join(current_seq)
                    # Parse header: >Bat|strain|accession
                    parts = current_id.split("|")
                    host = parts[0] if len(parts) > 0 else "Bat"
                    strain = parts[1] if len(parts) > 1 else current_id
                    acc = parts[2] if len(parts) > 2 else current_id
                    records.append({
                        "accession": acc,
                        "full_name": f"RVA/Bat/{strain}",
                        "host": "Bat",  # Will be updated from metadata
                        "country": "Unknown",
                        "year": "Unknown",
                        "segment": "VP4",
                        "sequence": seq,
                        "length": len(seq),
                        "source": "ref",
                    })
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
    
    if current_id and current_seq:
        seq = ''.join(current_seq)
        parts = current_id.split("|")
        host = parts[0] if len(parts) > 0 else "Bat"
        strain = parts[1] if len(parts) > 1 else current_id
        acc = parts[2] if len(parts) > 2 else current_id
        records.append({
            "accession": acc,
            "full_name": f"RVA/Bat/{strain}",
            "host": "Bat",
            "country": "Unknown",
            "year": "Unknown",
            "segment": "VP4",
            "sequence": seq,
            "length": len(seq),
            "source": "ref",
        })
    
    return records

def load_metadata_csv(segment):
    """Load the segment-specific metadata CSV to get constellation group and host"""
    csv_map = {
        "VP4":  "/mnt/user-uploads/CU_bat_only/VP4_metadata.csv",
        "VP6":  "/mnt/user-uploads/CU_bat_only/VP6/VP6_metadata.csv",
        "VP7":  "/mnt/user-uploads/CU_bat_only/VP7/VP7_metadata.csv",
        "NSP4": "/mnt/user-uploads/CU_bat_only/NSP4/NSP4_metadata.csv",
    }
    meta = {}
    if not os.path.exists(csv_map[segment]):
        return meta
    
    # VP4_metadata.csv has no header — use fixed column positions
    if segment == "VP4":
        with open(csv_map[segment], 'r') as f:
            reader = csv.reader(f)
            for row in reader:
                if len(row) < 5:
                    continue
                acc = row[0].strip()
                full_name = row[1].strip() if len(row) > 1 else ""
                host = row[2].strip() if len(row) > 2 else ""
                gc = row[3].strip() if len(row) > 3 else ""
                meta[acc] = {
                    "GC": gc,
                    "host_csv": host,
                    "full_name": full_name,
                }
        return meta
    
    # Other CSVs have headers
    with open(csv_map[segment], 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            # Get accession (handle different column names)
            acc = row.get("accession") or row.get("Ascession") or ""
            acc = acc.strip()
            if not acc:
                continue
            # Get constellation group (handle different column names)
            gc = row.get("GC") or row.get("Genotyoe_constellation") or ""
            gc = gc.strip()
            # Get host (handle different column names)
            host = row.get("host") or row.get("Bat_species") or ""
            host = host.strip()
            # Get full name
            full_name = row.get("full_name") or row.get("Full_name") or ""
            full_name = full_name.strip()
            
            meta[acc] = {
                "GC": gc,
                "host_csv": host,
                "full_name": full_name,
            }
    return meta

def normalize_host(h):
    """Normalize host species names"""
    h = h.strip()
    
    # First, map full common names to scientific names
    common_names = {
        'straw-coloured fruit bat': 'Eidolon helvum',
        'straw-colored fruit bat': 'Eidolon helvum',
        'african fruit bat': 'Eidolon helvum',
        'lesser horseshoe bat': 'Rhinolophus hipposideros',
    }
    h_lower = h.lower()
    for common, sci in common_names.items():
        if common in h_lower:
            return sci
    
    # Fix typos
    h = re.sub(r'(?i)Eidolom', 'Eidolon helvum', h)
    h = re.sub(r'(?i)Macronycteris gigas', 'Hipposideros gigas', h)
    
    # Remove trailing " bat" (from GenBank /host= qualifier)
    h = re.sub(r'\s+bat$', '', h, flags=re.IGNORECASE)
    
    # Map generic "bat" to Unknown
    h = re.sub(r'(?i)^bat$', 'Unknown bat', h)
    
    # Standardize abbreviations to full names
    std = {
        'R. hipposideros': 'Rhinolophus hipposideros',
        'R. blasii': 'Rhinolophus blasii',
        'R. euryale': 'Rhinolophus euryale',
        'R. simulator': 'Rhinolophus simulator',
        'R. aegyptiacus': 'Rousettus aegyptiacus',
        'R. leschenaulti': 'Rousettus leschenaulti',
        'H. gigas': 'Hipposideros gigas',
        'H. pomona': 'Hipposideros pomona',
        'T. melanopogon': 'Taphozous melanopogon',
        'T. mauritianus': 'Taphozous mauritianus',
        'S. kuhlii': 'Scotophilus kuhlii',
        'M. molossus': 'Molossus molossus',
        'G. soricina': 'Glossophaga soricina',
        'C. perspicillata': 'Carollia perspicillata',
        'E. helvum': 'Eidolon helvum',
        'Eidolon helvum': 'Eidolon helvum',
    }
    for pat, full in std.items():
        if h.lower() == pat.lower():
            return full
    # Try partial match
    for pat, full in std.items():
        if pat.lower() in h.lower():
            return full
    return h

def main():
    all_records = {}
    
    # Parse GenBank files for VP6, VP7, NSP4
    for segment, gb_path in GB_FILES.items():
        print(f"\nParsing {segment} GenBank: {gb_path}")
        records = parse_genbank(segment, gb_path)
        print(f"  Found {len(records)} records")
        all_records[segment] = records
    
    # Parse VP4 matched FASTA
    vp4_fasta = "/mnt/user-uploads/rva_cub_study/data/segments/Bat_VP4_matched.fasta"
    print(f"\nParsing VP4 FASTA: {vp4_fasta}")
    vp4_records = parse_vp4_fasta(vp4_fasta)
    print(f"  Found {len(vp4_records)} records")
    all_records["VP4"] = vp4_records
    
    # Parse study sequences
    print(f"\nParsing study sequences: {STUDY_FASTA}")
    study_seqs = parse_study_fasta(STUDY_FASTA)
    for seg, seqs in study_seqs.items():
        print(f"  {seg}: {len(seqs)} study sequences")
        if seg in all_records:
            all_records[seg].extend(seqs)
    
    # Load metadata CSVs and merge constellation group info
    for segment in all_records:
        meta = load_metadata_csv(segment)
        for rec in all_records[segment]:
            acc = rec["accession"]
            # For study sequences, accession is the NODE name
            # Try matching by accession or by full_name containing the accession
            if acc in meta:
                rec["GC"] = meta[acc]["GC"]
                # Use CSV host if GenBank host is Unknown or generic "Bat"
                csv_host = meta[acc]["host_csv"]
                if rec["host"] in ("Unknown", "Bat") and csv_host:
                    rec["host"] = csv_host
            else:
                # Try matching by the first part of NODE name (e.g., B10-BT-OAU-4)
                short_acc = acc.split("_NODE_")[0] if "_NODE_" in acc else acc
                if short_acc in meta:
                    rec["GC"] = meta[short_acc]["GC"]
                else:
                    rec["GC"] = "Unknown"
            
            # Normalize host
            rec["host"] = normalize_host(rec["host"])
            
            # For study sequences, set GC to Red
            if rec["source"] == "Study":
                rec["GC"] = "Red"
                rec["host"] = "Eidolon helvum"
        
        # Manual host overrides for sequences with generic "bat" host in GenBank
        # but identifiable from related sequences/publications
        host_overrides = {
            # MYAS33: same study as MSLH14 (He et al. 2013) → Rhinolophus hipposideros
            "KF649187.1": "Rhinolophus hipposideros",  # VP4
            "KJ020894.1": "Rhinolophus hipposideros",  # VP6
            "KF649188.1": "Rhinolophus hipposideros",  # VP7
            "KJ020890.1": "Rhinolophus hipposideros",  # NSP4
            # BatLi10: same study as BatLi09/BatLy03 (Esona et al.) → Eidolon helvum
            "KX268747.1": "Eidolon helvum",  # VP6
            "KX268748.1": "Eidolon helvum",  # VP7
            "KX268752.1": "Eidolon helvum",  # NSP4
            # KSA402: G25P[25], Saudi Arabia → Eidolon helvum (G25 associated with E. helvum)
            "KX420942.1": "Eidolon helvum",  # VP4
        }
        # GC overrides
        gc_overrides = {
            "KX420942.1": "Yellow",  # G25P[25], Saudi Arabia
        }
        
        for rec in all_records[segment]:
            acc = rec["accession"]
            if acc in host_overrides:
                rec["host"] = host_overrides[acc]
            if acc in gc_overrides:
                rec["GC"] = gc_overrides[acc]
    
    # Save per-segment FASTA and metadata CSV
    for segment, records in all_records.items():
        # Save FASTA
        fasta_path = os.path.join(OUT_DIR, f"{segment}_cds.fasta")
        with open(fasta_path, 'w') as f:
            for rec in records:
                seq = rec["sequence"].upper().replace("-", "")
                f.write(f">{rec['accession']}\n{seq}\n")
        
        # Save metadata CSV
        csv_path = os.path.join(OUT_DIR, f"{segment}_metadata_parsed.csv")
        with open(csv_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=["accession", "full_name", "host", "GC", "source", "segment", "length", "country", "year"])
            writer.writeheader()
            for rec in records:
                writer.writerow({
                    "accession": rec["accession"],
                    "full_name": rec["full_name"],
                    "host": rec["host"],
                    "GC": rec["GC"],
                    "source": rec["source"],
                    "segment": rec["segment"],
                    "length": rec["length"],
                    "country": rec["country"],
                    "year": rec["year"],
                })
        
        # Print summary
        hosts = {}
        gcs = {}
        for rec in records:
            h = rec["host"]
            g = rec["GC"]
            hosts[h] = hosts.get(h, 0) + 1
            gcs[g] = gcs.get(g, 0) + 1
        
        print(f"\n{segment}: {len(records)} sequences")
        print(f"  Hosts: {hosts}")
        print(f"  Groups: {gcs}")
        print(f"  Saved: {fasta_path}, {csv_path}")

if __name__ == "__main__":
    main()
