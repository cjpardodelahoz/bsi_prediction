#!/bin/bash

#SBATCH --output=log/enterococcus_diversity_II/enterococcus_isabl1_checkm2.out
#SBATCH --error=log/enterococcus_diversity_II/enterococcus_isabl1_checkm2.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --partition=cpu

# Activate the CheckM2 environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm2

# Define paths
MAGS_DIR="analyses/enterococcus_diversity_II/genomes/mags/isabl1"
OUTPUT_DIR="analyses/enterococcus_diversity_II/genomes/checkm2/isabl1"

checkm2 predict \
    --input ${MAGS_DIR} \
    --output-directory ${OUTPUT_DIR} \
    --remove_intermediates \
    --threads 48 \