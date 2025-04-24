#!/bin/bash

#SBATCH --array=1-395%10  # Adjust the range (1-395) based on the number of libraries
#SBATCH --mem=64G         # Memory per task
#SBATCH -c 16             # Number of CPU cores per task
#SBATCH --time=120:00:00   # Maximum runtime
#SBATCH --error=log/bacteroides_variation/yan_sd_2022_metassembly_%A_%a.err
#SBATCH --output=log/bacteroides_variation/yan_sd_2022_metassembly_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Set paths
READS_DIR="data/reads/yan_sd_2022"
FASTP_DIR="analyses/yan_sd_2022/fastp"
ASSEMBLY_DIR="analyses/yan_sd_2022/assembly/metaplasmid"

# Get the sample name for the current task
sample_name=$(ls ${READS_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directories
mkdir -p ${ASSEMBLY_DIR}/${sample_name}

# Step 1: Assemble the metagenome using metaplasmidSPAdes
conda activate spades
spades.py --metaplasmid \
          -1 ${FASTP_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz \
          -2 ${FASTP_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz \
          -o ${ASSEMBLY_DIR}/${sample_name} \
          -k 21,33,55,75,95 \
          -t 16 \
          -m 64
conda deactivate