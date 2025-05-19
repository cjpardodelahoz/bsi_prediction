#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/enterococcus_checkm2.out
#SBATCH --error=log/enterococcus_diversity/enterococcus_checkm2.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Activate the CheckM2 environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm2

# Define paths
MAGS_DIR="analyses/enterococcus_diversity/genomes/mags"
OUTPUT_DIR="analyses/enterococcus_diversity/genomes/checkm2"

checkm2 predict \
    --input ${MAGS_DIR} \
    --output-directory ${OUTPUT_DIR} \
    --remove_intermediates \
    --threads 12 \