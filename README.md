# Genomic Characterisation and Phylogenetic Analysis of Rotavirus A in Nigerian Bats

This repository contains the analysis scripts, code, and documentation for the bioinformatic and phylogenetic workflows described in our manuscript on the evolutionary dynamics and reassortment of Rotavirus A (RVA) in Nigerian bat populations.

## Overview

The study involved metagenomic sequencing of bat rectal swabs to identify and characterise RVA genomes. The analysis is broadly split into two parts:
1. **Upstream Analysis:** Read filtering, de novo assembly, and RVA contig identification.
2. **Downstream Analysis:** Genotyping, sequence alignment, phylogenetic inference, and codon usage analysis.

*Note: Raw sequencing data and contigs are be made available via public repositories.*

---

## 1. Upstream Analysis Pipeline

Raw metagenomic sequencing reads were processed using the **Virome Paired-End Reads (ViPER)** pipeline. 
The pipeline code is publicly available here: [Matthijnssenslab/ViPER](https://github.com/Matthijnssenslab/ViPER).

Briefly, the upstream workflow included:
- **Quality Control & Trimming:** Trimmomatic was used to remove sequencing adapters and low-quality bases.
- **Host Depletion:** Host-derived reads were filtered out using Bowtie2.
- **De Novo Assembly:** Filtered reads were assembled using metaSPAdes and MEGAHIT.
- **Contig Annotation:** Contigs were annotated via DIAMOND (sensitive mode) and confirmed using BLASTn against reference databases. 
- **Read Mapping:** Trimmed reads were mapped back to assembled RVA contigs using Bowtie2 to determine genome coverage and assembly depth.

---

## 2. Downstream Analysis & Phylogenetics

The scripts for the downstream characterisation, phylogenetic analysis, and codon usage assessment are provided in this repository.

### Prerequisites & Dependencies
To run the phylogenetic pipeline, ensure the following tools are installed:
- [MAFFT](https://mafft.cbrc.jp/alignment/software/) (v7.x)
- [AliView](https://ormbunkar.se/aliview/) (for manual alignment inspection/trimming)
- [BMGE] (Block Mapping and Gathering with Entropy)
- [IQ-TREE](http://www.iqtree.org/) (v2.x)
- **R** (>= 4.0.0) with packages: `ggtree`, `ggplot2`, `seqinr`, etc. (Check individual scripts for specific library requirements).
- **Python** (>= 3.x) with `Biopython`.

### Repository Structure

- `Genomic_characteriszation_and_and Phylogenetic_scripts/`
  - `00_phylogenetic_pipeline_commands.sh`: A shell script documenting the exact command-line arguments used for alignment (MAFFT), automated trimming (BMGE), and Maximum Likelihood tree construction (IQ-TREE).
  - `extract_rotavirus_segments_csv.py`: Script to parse and extract segment metadata.
  - `Genetic_distance_segments2_11.R`: R script to calculate pairwise genetic distances across segments.
  - `Tree_Plotting_Scripts/`: R scripts utilizing `ggtree` to visualize and annotate the IQ-TREE outputs for each of the 11 genome segments (e.g., `Plot_Tree_VP1.R`, `Plot_Composite_Panels.R`).

- `Codon_Usage_scripts/`
  - Python scripts (`00a_download_genome_cds.py`, `01a_parse_genbank.py`, etc.) for downloading sequences, parsing GenBank files, and preparing transcriptome tables.
  - R scripts (`02_rscu_analysis.R`, `03_gc_enc_analysis.R`, `04_cai_analysis.R`, `05_tanglegram_ari.R`, etc.) for calculating Relative Synonymous Codon Usage (RSCU), Effective Number of Codons (ENC), GC3s bias, Codon Adaptation Index (CAI), and plotting UPGMA clustering and tanglegrams.

### Reproducibility Steps

1. **Genotyping and ORF Prediction**
   - Assembled contigs were assigned to segments using the [RIVM rotavirus genotyping tool](https://mpf.rivm.nl/mpf/typingtool/rotavirusa/).
   - Open Reading Frames (ORFs) were predicted via NCBI ORFfinder, retaining only those >30% of the expected gene length.

2. **Phylogenetic Inference**
   - See `Genomic_characteriszation_and_and Phylogenetic_scripts/00_phylogenetic_pipeline_commands.sh` for the exact line-by-line commands.
   - Alignments were performed with `mafft --auto`.
   - Alignments were trimmed either manually in AliView (for conserved segments) or automatically with `bmge -t DNA` (for divergent segments like NSP1 and NSP5).
   - ML trees were inferred using `iqtree -m MFP -B 1000 -alrt 1000`.

3. **Tree Visualization**
   - Run the respective segment plotting scripts located in `Tree_Plotting_Scripts/` (e.g., `Rscript Plot_Tree_VP6.R`) to generate the annotated phylogenetic trees.

4. **Codon Usage Analysis**
   - Navigate to `Codon_Usage_scripts/`.
   - Execute the R scripts sequentially (from `02_rscu_analysis.R` through `10_enc_deviation.R`) to reproduce the nucleotide composition metrics, RSCU heatmaps, and UPGMA clustering visualizations.

---

## Citation
If you use this code or data, please cite the associated manuscript.
