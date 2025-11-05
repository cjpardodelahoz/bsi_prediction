#!/bin/bash

#SBATCH --output=log/enterococcus_diversity_II/mlst.out
#SBATCH --error=log/enterococcus_diversity_II/mlst.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=cpu

# Activate the conda environment with mlst installed
source $(conda info --base)/etc/profile.d/conda.sh
conda activate mlst

# Define paths
GENOMES_DIR="analyses/enterococcus_diversity_II/genomes/mags/isabl1"
OUTPUT_FILE="analyses/enterococcus_diversity_II/genomes/mlst/mlst_2002_isabl1.csv"
mkdir -p $(dirname ${OUTPUT_FILE})

# Run mlst for all .fna files in the directory
mlst ${GENOMES_DIR}/*.fna --csv --threads 4 > ${OUTPUT_FILE}