#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/drep_98.out
#SBATCH --error=log/enterococcus_diversity/drep_98.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Activate the dRep environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate drep2

# Define paths
MAGS_DIR="analyses/enterococcus_diversity/genomes/mags"
#CHECKM_RESULTS="analyses/enterococcus_diversity/genomes/checkm2/quality_report.csv"
OUTPUT_DIR="analyses/enterococcus_diversity/genomes/drep"

# Create the output directory if it doesn't exist
mkdir -p ${OUTPUT_DIR}

# Run dRep with a 98% ANI threshold and precomputed CheckM2 results
dRep dereplicate ${OUTPUT_DIR} \
    -g ${MAGS_DIR}/*.fna \
    -p 16 \
    --S_algorithm ANImf \
    --completeness 70 \
    --contamination 5 \
    --cov_thresh 0.1 \
    --S_ani 0.98
#    --genomeInfo ${CHECKM_RESULTS}

    