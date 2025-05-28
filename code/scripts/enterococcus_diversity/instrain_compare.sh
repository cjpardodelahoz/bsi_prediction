#!/bin/bash

#SBATCH --output=log/enterococcus_diversity/instrain_compare.out
#SBATCH --error=log/enterococcus_diversity/instrain_compare.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=1024G
#SBATCH --partition=cpu_highmem

# Load conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate instrain

# Output directory for inStrain comparison results
OUTPUT_DIR="analyses/enterococcus_diversity/strains/instrain_compare"

# inStrain profiles
PROFILES=$(find analyses/enterococcus_diversity/strains/instrain -mindepth 2 -maxdepth 2 -type d -name instrain_profile)

# Genome list
GENOME_LIST=$(ls analyses/yan_sd_2022/binning/drep98/dereplicated_genomes/*.fna)

# Run inStrain compare
inStrain compare -i ${PROFILES} \
    -o ${OUTPUT_DIR} \
    --stb $(echo ${GENOME_LIST}) \
    -p 16
