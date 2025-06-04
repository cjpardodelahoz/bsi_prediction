#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/ggcaller_align.out
#SBATCH --error=log/enterococcus_diversity/ggcaller_align.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Load conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ggcaller

# List of genomes
GENOME_LIST="analyses/enterococcus_diversity/phylogenetics/mags_plus_lebreton/genome_list.txt"

# Output directory
OUTPUT_DIR="analyses/enterococcus_diversity/phylogenetics/mags_plus_lebreton/ggcaller"
mkdir -p $OUTPUT_DIR

# Run ggcaller
ggcaller --refs ${GENOME_LIST} \
         --out ${OUTPUT_DIR} \
         --alignment core \
         --aligner def \
         --core-threshold 0.99 \
         --ignore-pseduogenes \
         --clean-mode strict \
         --threads 16 \
         --annotation none