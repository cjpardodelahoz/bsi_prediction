#!/bin/bash

#SBATCH --output=log/bacteroides_pul/checkm2_isabl1.out
#SBATCH --error=log/bacteroides_pul/checkm2_isabl1.err
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --partition=cpu

# Activate the CheckM2 environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate checkm2

# Define paths
MAGS_DIR="analyses/bacteroides_pul/mags"
OUTPUT_DIR="analyses/bacteroides_pul/checkm2"

checkm2 predict \
    --input ${MAGS_DIR} \
    --output-directory ${OUTPUT_DIR} \
    --remove_intermediates \
    --threads 48 \