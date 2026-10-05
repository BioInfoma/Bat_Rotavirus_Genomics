#!/bin/bash
# 00c_prep_trinity_inputs.sh
# Download RNA-Seq reads for 5 bat species for Trinity transcriptome assembly
# Downloads 1 run per species (sufficient for reliable RSCU)
#
# Usage: bash 00c_prep_trinity_inputs.sh

set -e

export PATH="/workspace/sratoolkit.3.2.1-ubuntu64/bin:$PATH"
DL_DIR="/workspace/trinity_data"
mkdir -p "$DL_DIR"

# ---- Species 1: R. simulator (not on ENA, use NCBI SRA Toolkit) ----
echo "=== Downloading R. simulator (SRR38065611) via NCBI ==="
cd "$DL_DIR"
mkdir -p Rsimulator
cd Rsimulator
if [ ! -f SRR38065611_1.fastq.gz ]; then
    prefetch SRR38065611 -o SRR38065611.sra 2>&1
    fasterq-dump SRR38065611.sra --split-files --threads 4 2>&1
    gzip -f SRR38065611_1.fastq SRR38065611_2.fastq 2>/dev/null || true
    rm -f SRR38065611.sra
fi
echo "R. simulator done: $(ls -la *.fastq.gz 2>/dev/null | wc -l) files"

# ---- Species 2: T. melanopogon (ENA, smallest run) ----
echo "=== Downloading T. melanopogon (SRR6796937) via ENA ==="
cd "$DL_DIR"
mkdir -p Tmelanopogon
cd Tmelanopogon
if [ ! -f SRR6796937_1.fastq.gz ]; then
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR679/007/SRR6796937/SRR6796937_1.fastq.gz" -O SRR6796937_1.fastq.gz
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR679/007/SRR6796937/SRR6796937_2.fastq.gz" -O SRR6796937_2.fastq.gz
fi
echo "T. melanopogon done: $(ls -la *.fastq.gz 2>/dev/null | wc -l) files"

# ---- Species 3: R. leschenaulti (ENA) ----
echo "=== Downloading R. leschenaulti (SRR835441) via ENA ==="
cd "$DL_DIR"
mkdir -p Rleschenaulti
cd Rleschenaulti
if [ ! -f SRR835441_1.fastq.gz ]; then
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR835/SRR835441/SRR835441_1.fastq.gz" -O SRR835441_1.fastq.gz
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR835/SRR835441/SRR835441_2.fastq.gz" -O SRR835441_2.fastq.gz
fi
echo "R. leschenaulti done: $(ls -la *.fastq.gz 2>/dev/null | wc -l) files"

# ---- Species 4: G. soricina (ENA, smallest run) ----
echo "=== Downloading G. soricina (SRR13300776) via ENA ==="
cd "$DL_DIR"
mkdir -p Gsoricina
cd Gsoricina
if [ ! -f SRR13300776_1.fastq.gz ]; then
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR133/007/SRR13300776/SRR13300776_1.fastq.gz" -O SRR13300776_1.fastq.gz
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR133/007/SRR13300776/SRR13300776_2.fastq.gz" -O SRR13300776_2.fastq.gz
fi
echo "G. soricina done: $(ls -la *.fastq.gz 2>/dev/null | wc -l) files"

# ---- Species 5: C. perspicillata (ENA, smallest run) ----
# Note: C. perspicillata files on ENA are interleaved (single file)
# Need to split into _1 and _2 for Trinity
echo "=== Downloading C. perspicillata (SRR30428459) via ENA ==="
cd "$DL_DIR"
mkdir -p Cperspicillata
cd Cperspicillata
if [ ! -f SRR30428459_1.fastq.gz ]; then
    wget -q "https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR304/009/SRR30428459/SRR30428459.fastq.gz" -O SRR30428459.fastq.gz
    # Split interleaved FASTQ into R1 and R2
    echo "Splitting interleaved FASTQ..."
    zcat SRR30428459.fastq.gz | paste - - - - - - - - | \
        awk -F'\t' '{print $1"\n"$2"\n"$3"\n"$4 > "SRR30428459_1.fastq"; print $5"\n"$6"\n"$7"\n"$8 > "SRR30428459_2.fastq"}'
    gzip -f SRR30428459_1.fastq SRR30428459_2.fastq
    rm -f SRR30428459.fastq.gz
fi
echo "C. perspicillata done: $(ls -la *.fastq.gz 2>/dev/null | wc -l) files"

echo ""
echo "=== ALL DOWNLOADS COMPLETE ==="
du -sh "$DL_DIR"/*
