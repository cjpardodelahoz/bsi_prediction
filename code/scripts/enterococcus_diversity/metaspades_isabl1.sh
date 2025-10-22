#!/bin/bash

#SBATCH --array=1-2077%50           # 2077 samples, 20 concurrent jobs
#SBATCH --mem=48G                  # Memory per task
#SBATCH -c 16                       # Number of CPU cores per task
#SBATCH --time=120:00:00            # Maximum runtime
#SBATCH --error=log/enterococcus_diversity/metaspades_isabl1_%A_%a.err
#SBATCH --output=log/enterococcus_diversity/metaspades_isabl1_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Set output directories
FASTP_DIR="analyses/enterococcus_diversity/metagenomes/fastp"
ASSEMBLY_DIR="analyses/enterococcus_diversity/metagenomes/assembly/metaspades"

# Get sample name (with brackets and quotes) for current task
sample_bracketed=$(sed -n "${SLURM_ARRAY_TASK_ID}p" misc_files/enterococcus_diversity/all_target_samples_hostdepleted.txt)

# Remove brackets and single quotes: e.g., ['275C'] -> 275C
sample=$(echo $sample_bracketed | sed "s/\[\('\([^']*\)'\)\]/\2/")

# Get first record for this sample from the paths table
record=$(awk -v s="$sample_bracketed" '$2==s{print; exit}' misc_files/enterococcus_diversity/all_target_samples_paths.txt)

# Extract R1 and R2 paths (columns 3 and 4)
R1=$(echo "$record" | cut -f3)
R2=$(echo "$record" | cut -f4)

# Replace storage path
R1=$(echo "$R1" | sed 's|/data/brinkvd/|/data1/collab004/|')
R2=$(echo "$R2" | sed 's|/data/brinkvd/|/data1/collab004/|')

# Create output directories
mkdir -p ${FASTP_DIR}/${sample} ${ASSEMBLY_DIR}/${sample}

# Step 1: Run fastp for quality control and trimming
conda activate fastp
fastp -i ${R1} -I ${R2} \
      -o ${FASTP_DIR}/${sample}/${sample}_R1_trimmed.fastq.gz \
      -O ${FASTP_DIR}/${sample}/${sample}_R2_trimmed.fastq.gz \
      --html ${FASTP_DIR}/${sample}/${sample}_fastp.html \
      --json ${FASTP_DIR}/${sample}/${sample}_fastp.json \
      --cut_front \
      --cut_right \
      --cut_mean_quality 20 \
      --length_required 30 \
      --thread 4
conda deactivate

# Step 2: Assemble the metagenome using metaSPAdes
conda activate spades
spades.py --meta \
          -1 ${FASTP_DIR}/${sample}/${sample}_R1_trimmed.fastq.gz \
          -2 ${FASTP_DIR}/${sample}/${sample}_R2_trimmed.fastq.gz \
          -o ${ASSEMBLY_DIR}/${sample} \
          -k 21,33,55,75,95 \
          -t 16 \
          -m 48
conda deactivate

# Step 3: Remove error-corrected reads and fastp-trimmed reads to save space
rm -r ${ASSEMBLY_DIR}/${sample}/K*/* 
rm -r ${ASSEMBLY_DIR}/${sample}/corrected/* 
rm -r ${FASTP_DIR}/${sample}/${sample}_R1_trimmed.fastq.gz
rm -r ${FASTP_DIR}/${sample}/${sample}_R2_trimmed.fastq.gz