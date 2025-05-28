#!/bin/bash

#SBATCH --array=1-216%15  # This is the number of samples with E. faecium detection
#SBATCH --output=log/enterococcus_diversity/floria_phasing_%A_%a.out
#SBATCH --error=log/enterococcus_diversity/floria_phasing_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Prepare to activate conda environments
source $(conda info --base)/etc/profile.d/conda.sh

# Define paths
READS_DIR="data/reads/yan_sd_2022"
REF_DIR="analyses/enterococcus_diversity/strains/floria"
OUTPUT_DIR="analyses/enterococcus_diversity/strains/floria"
SAMPLES_WITH_EFAECIUM=$(awk -F'\t' '$2 == "s_1044K_1_meta_vamb_single" && $4 == "Positive" {print $1}' analyses/enterococcus_diversity/detection/detection_results.tsv)
SAMPLE_NAME=$(echo ${SAMPLES_WITH_EFAECIUM} | cut -d " " -f ${SLURM_ARRAY_TASK_ID})

# Input files
R1="${READS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R1_all.fastq.gz"
R2="${READS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R2_all.fastq.gz"
REF_FASTA="analyses/enterococcus_diversity/strains/ref_contigs/ref_contigs.fna"

# Output files
BAM_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}/mapped_reads/${SAMPLE_NAME}.sorted.bam"
VCF_FILE="${OUTPUT_DIR}/${SAMPLE_NAME}/variants/${SAMPLE_NAME}.vcf"

# Floria output path
FLORIA_OUTPUT="${OUTPUT_DIR}/${SAMPLE_NAME}/floria_output"

# E. faecium reference contigs
REF_CONTIGS=$(grep '^>' analyses/enterococcus_diversity/genomes/drep/dereplicated_genomes/s_1044K_1_meta_vamb_single.fna | sed 's/^>//' | tr '\n' ' ')

# Create output directories
mkdir -p ${OUTPUT_DIR}/${SAMPLE_NAME}/mapped_reads
mkdir -p ${OUTPUT_DIR}/${SAMPLE_NAME}/variants

# Step 1: Index the reference FASTA file
conda activate bowtie2
if [[ ! -f ${REF_FASTA}.fai ]]; then
    samtools faidx ${REF_FASTA}
fi

# Step 2: Build Bowtie2 index for the reference contigs
if [[ ! -f ${REF_FASTA}.1.bt2 ]]; then
    bowtie2-build ${REF_FASTA} ${REF_FASTA}
fi

# Step 3: Map reads to the reference contigs using Bowtie2
bowtie2 -x ${REF_FASTA} -1 ${R1} -2 ${R2} -p 16 | samtools view -bS - | samtools sort -o ${BAM_FILE}
samtools index ${BAM_FILE}

# Step 4: Call variants using FreeBayes
conda activate floria
freebayes -f ${REF_FASTA} -F 0.01 -C 1 --pooled-continuous --bam ${BAM_FILE} > ${VCF_FILE}

# Step 5: Run floria
floria -b ${BAM_FILE} \
    -v ${VCF_FILE} \
    -o ${FLORIA_OUTPUT} \
    -t 16 \
    -r ${REF_FASTA} \
    --contigs ${REF_CONTIGS}

# Step 6: Clean up
rm ${BAM_FILE}
rm ${BAM_FILE}.bai