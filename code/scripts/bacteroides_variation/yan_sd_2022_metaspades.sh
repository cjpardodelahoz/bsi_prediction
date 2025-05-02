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

# # Set paths
READS_DIR="data/reads/yan_sd_2022"
FASTP_DIR="analyses/yan_sd_2022/fastp"
ASSEMBLY_DIR="analyses/yan_sd_2022/assembly/metaspades"

# Get the sample name for the current task
sample_name=$(ls ${READS_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directories
mkdir -p ${FASTP_DIR}/${sample_name} ${ASSEMBLY_DIR}/${sample_name}

# Input files
R1="${READS_DIR}/${sample_name}/${sample_name}_R1_all.fastq.gz"
R2="${READS_DIR}/${sample_name}/${sample_name}_R2_all.fastq.gz"

# Step 1: Run fastp for quality control and trimming
conda activate fastp
fastp -i ${R1} -I ${R2} \
      -o ${FASTP_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz \
      -O ${FASTP_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz \
      --html ${FASTP_DIR}/${sample_name}/${sample_name}_fastp.html \
      --json ${FASTP_DIR}/${sample_name}/${sample_name}_fastp.json \
      --cut_front \
      --cut_right \
      --cut_mean_quality 20 \
      --length_required 30 \
      --thread 4
conda deactivate

# Step 2: Assemble the metagenome using metaSPAdes
conda activate spades
spades.py --meta \
          -1 ${FASTP_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz \
          -2 ${FASTP_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz \
          -o ${ASSEMBLY_DIR}/${sample_name} \
          -k 21,33,55,75,95 \
          -t 16 \
          -m 64
conda deactivate

# Step 3: Remove error-corrected reads and fastp-trimmed reads to save space
rm -r ${ASSEMBLY_DIR}/${sample_name}/K*/* 
rm -r ${ASSEMBLY_DIR}/${sample_name}/corrected/*  # Remove error-corrected reads
rm -r ${FASTP_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz
rm -r ${FASTP_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz