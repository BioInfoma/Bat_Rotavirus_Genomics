import os
import pandas as pd
from Bio import SeqIO
import re

def extract_segment_number(blast_result):
    """Extract segment number from BLAST result column"""
    if isinstance(blast_result, str):
        match = re.search(r'segment\s*(\d+)', blast_result, re.IGNORECASE)
        if match:
            return int(match.group(1))
    return None

def process_sample(csv_path, fasta_path, segments_dict):
    """Process a single sample and store sequences."""
    
    df = pd.read_csv(csv_path)
    fasta_sequences = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    
    for index, row in df.iterrows():
        contig_name_csv = row['name']
        segment_num = extract_segment_number(row['BLAST result'])
        
        if segment_num and 1 <= segment_num <= 11:
            matching_sequence = None
            
            for fasta_header, sequence in fasta_sequences.items():
                if fasta_header.startswith(contig_name_csv):
                    matching_sequence = sequence
                    break
            
            if matching_sequence:
                segment_key = f"segment{segment_num}"
                if segment_key not in segments_dict:
                    segments_dict[segment_key] = []
                
                segments_dict[segment_key].append(matching_sequence)
                print(f"✓ Matched: {contig_name_csv} -> {segment_key}")
            else:
                print(f"✗ No FASTA match found for: {contig_name_csv}")

def main():
    csv_folder = "RVA_blast_results"
    fasta_folder = "RotavirusA_samplefastas" 
    output_folder = "segments_fasta"
    
    os.makedirs(output_folder, exist_ok=True)
    
    segments_dict = {}
    
    csv_files = [f for f in os.listdir(csv_folder) if f.endswith('.csv')]
    
    print(f"Found {len(csv_files)} CSV files to process...")
    
    for csv_file in csv_files:
        sample_name = csv_file.replace('.csv', '')
        csv_path = os.path.join(csv_folder, csv_file)
        fasta_path = os.path.join(fasta_folder, f"{sample_name}.fasta")
        
        print(f"\nProcessing sample: {sample_name}")
        
        if os.path.exists(fasta_path):
            process_sample(csv_path, fasta_path, segments_dict)
        else:
            print(f"✗ FASTA file not found: {fasta_path}")
    
    print(f"\nWriting output FASTA files to {output_folder}/...")
    
    for segment, sequences in segments_dict.items():
        output_path = os.path.join(output_folder, f"{segment}.fasta")
        
        if sequences:
            SeqIO.write(sequences, output_path, "fasta")
            print(f"✓ Created {output_path} with {len(sequences)} sequences")
        else:
            print(f"✗ No sequences found for {segment}")
    
    print(f"\n=== PROCESSING COMPLETE ===")
    for segment in sorted(segments_dict.keys()):
        count = len(segments_dict[segment])
        print(f"{segment}: {count} sequences")

if __name__ == "__main__":
    main()
