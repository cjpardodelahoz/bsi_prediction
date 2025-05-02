#!/bin/bash

#SBATCH --array=1-395%10  # Adjust the range (1-395) based on the number of libraries
#SBATCH --mem-per-cpu=8G         # Memory per task
#SBATCH -c 4             # Number of CPU cores per task
#SBATCH --time=120:00:00  # Maximum runtime
#SBATCH --error=log/bacteroides_variation/classify_assemble_%A_%a.err
#SBATCH --output=log/bacteroides_variation/classify_assemble_%A_%a.out
#SBATCH --partition=cpu

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Set paths
READS_DIR="data/reads/yan_sd_2022"
CLASSIFIED_DIR="analyses/bacteroides_variation/genomes/kraken_binning/kraken_reports"
BACTEROIDES_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_raw_reads"
TRIMMED_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_trimmed_reads"
ASSEMBLY_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_metaspades"
KRAKEN_DB="/data1/xavierj/carlos/dbs/kraken/18042025"

# Get the sample name for the current task
sample_name=$(ls ${READS_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directories
rm -r ${CLASSIFIED_DIR}/${sample_name} ${BACTEROIDES_DIR}/${sample_name} ${TRIMMED_DIR}/${sample_name} ${ASSEMBLY_DIR}/${sample_name}
mkdir -p ${CLASSIFIED_DIR}/${sample_name} ${BACTEROIDES_DIR}/${sample_name} ${TRIMMED_DIR}/${sample_name} ${ASSEMBLY_DIR}/${sample_name}

# Input files
R1="${READS_DIR}/${sample_name}/${sample_name}_R1_all.fastq.gz"
R2="${READS_DIR}/${sample_name}/${sample_name}_R2_all.fastq.gz"

# Step 1: Classify reads using Kraken2
conda activate kraken2
kraken2 --db ${KRAKEN_DB} \
        --paired ${R1} ${R2} \
        --threads 4 \
        --report ${CLASSIFIED_DIR}/${sample_name}/${sample_name}_kraken_report.txt \
        --output ${CLASSIFIED_DIR}/${sample_name}/${sample_name}_kraken_output.txt

# Step 2: Extract reads classified as Bacteroides
extract_kraken_reads.py -k ${CLASSIFIED_DIR}/${sample_name}/${sample_name}_kraken_output.txt \
                        -r ${CLASSIFIED_DIR}/${sample_name}/${sample_name}_kraken_report.txt \
                        -s ${R1} -s2 ${R2} \
                        -o ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R1.fastq \
                        -o2 ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R2.fastq \
                        --fastq-output \
                        --include-children --taxid 816  # TaxID for Bacteroides
conda deactivate

# Step 3: Clean up extracted Bacteroides reads with fastp
conda activate fastp
fastp -i ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R1.fastq \
      -I ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R2.fastq \
      -o ${TRIMMED_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz \
      -O ${TRIMMED_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz \
      --html ${TRIMMED_DIR}/${sample_name}/${sample_name}_fastp.html \
      --json ${TRIMMED_DIR}/${sample_name}/${sample_name}_fastp.json \
      --cut_front \
      --cut_right \
      --cut_mean_quality 20 \
      --length_required 30 \
      --thread 4
conda deactivate

# Delete untrimmed Bacteroides reads to save space
rm -f ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R1.fastq
rm -f ${BACTEROIDES_DIR}/${sample_name}/${sample_name}_R2.fastq

# Step 4: Assemble the cleaned Bacteroides reads using metaSPAdes
conda activate spades
spades.py --meta \
          -1 ${TRIMMED_DIR}/${sample_name}/${sample_name}_R1_trimmed.fastq.gz \
          -2 ${TRIMMED_DIR}/${sample_name}/${sample_name}_R2_trimmed.fastq.gz \
          -o ${ASSEMBLY_DIR}/${sample_name} \
          -k 21,33,55,75,95 \
          -t 4 \
          -m 32
conda deactivate

echo "Pipeline completed for sample: ${sample_name}"

# Remove error-corrected reads and Kmer directories to save space
rm -f ${ASSEMBLY_DIR}/${sample_name}/corrected/*.fastq.gz
rm -rf ${ASSEMBLY_DIR}/${sample_name}/K*/*