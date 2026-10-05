#!/bin/bash
# -----------------------------------------------------------------------------
# Rotavirus A (RVA) Phylogenetic Pipeline Commands
# -----------------------------------------------------------------------------
# This script documents the line-by-line commands used for the downstream
# phylogenetic analysis of Rotavirus A contigs as described in the manuscript.
#
# Prerequisites:
# - MAFFT (v7.x)
# - AliView (used for manual inspection/trimming)
# - BMGE (Block Mapping and Gathering with Entropy)
# - IQ-TREE (v2.x)
# -----------------------------------------------------------------------------

# Note: The input sequences should be unaligned FASTA files containing both 
# the study sequences and globally representative reference sequences for 
# each of the 11 genome segments.

SEGMENT="VP1" # Example segment. Iterate for VP1-VP7, NSP1-NSP5.
INPUT_FASTA="${SEGMENT}_unaligned.fasta"
ALIGNED_FASTA="${SEGMENT}_aligned.fasta"
TRIMMED_FASTA="${SEGMENT}_trimmed.fasta"

# -----------------------------------------------------------------------------
# Step 1: Multiple Sequence Alignment
# -----------------------------------------------------------------------------
# We used MAFFT with the --auto parameter, which automatically selects the 
# most appropriate alignment strategy based on sequence length and complexity.

echo "Running MAFFT alignment for ${SEGMENT}..."
mafft --auto $INPUT_FASTA > $ALIGNED_FASTA

# -----------------------------------------------------------------------------
# Step 2: Alignment Trimming
# -----------------------------------------------------------------------------
# Alignments were trimmed to the start and end of the mature sequence. 
# While this was done manually using AliView for conserved segments, highly 
# divergent segments (like NSP1 and NSP5) were trimmed using BMGE.
#
# Below is the BMGE command used for automated trimming:

echo "Running BMGE trimming for ${SEGMENT}..."
# -t DNA specifies nucleotide sequences
bmge -i $ALIGNED_FASTA -t DNA -o $TRIMMED_FASTA

# -----------------------------------------------------------------------------
# Step 3: Model Selection and Tree Construction
# -----------------------------------------------------------------------------
# Phylogenetic trees were constructed using Maximum Likelihood in IQ-TREE.
# -m MFP : ModelFinder (automatically selects the best-fit substitution model)
# -B 1000 : 1,000 ultrafast bootstrap replicates (UFBoot2)
# -alrt 1000 : 1,000 SH-aLRT replicates to assess branch support

echo "Building phylogenetic tree for ${SEGMENT} with IQ-TREE..."
iqtree -s $TRIMMED_FASTA -m MFP -B 1000 -alrt 1000

echo "Pipeline complete for ${SEGMENT}. Output tree is in ${TRIMMED_FASTA}.treefile"
