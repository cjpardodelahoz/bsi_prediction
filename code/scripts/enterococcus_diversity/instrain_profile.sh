#!/bin/bash

#SBATCH --array=1-216%15  # Adjust to the number of positive samples
#SBATCH --output=log/enterococcus_diversity/instrain_profile_%A_%a.out
#SBATCH --error=log/enterococcus_diversity/instrain_profile_%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

source $(conda info --base)/etc/profile.d/conda.sh

# Paths
READS_DIR="data/reads/yan_sd_2022"
REF_FASTA="analyses/enterococcus_diversity/strains/ref_contigs/ref_contigs.fna"
OUTPUT_DIR="analyses/enterococcus_diversity/strains/instrain"

# Get sample list and select sample for this task
SAMPLES_WITH_EFAECIUM=$(awk -F'\t' '$2 == "s_1044K_1_meta_vamb_single" && $4 == "Positive" {print $1}' analyses/enterococcus_diversity/detection/detection_results.tsv)
SAMPLE_NAME=$(echo ${SAMPLES_WITH_EFAECIUM} | cut -d " " -f ${SLURM_ARRAY_TASK_ID})

# Input files
R1="${READS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R1_all.fastq.gz"
R2="${READS_DIR}/${SAMPLE_NAME}/${SAMPLE_NAME}_R2_all.fastq.gz"

# Output files
SAMPLE_OUTDIR="${OUTPUT_DIR}/${SAMPLE_NAME}"
MAPPING_DIR="${SAMPLE_OUTDIR}/mapping"
BAM_FILE="${MAPPING_DIR}/${SAMPLE_NAME}.sorted.bam"
mkdir -p ${MAPPING_DIR}

# E. faecium reference contigs
REF_CONTIGS="analyses/enterococcus_diversity/strains/ref_contigs/ref_contig_list.txt"

# Genome list
GENOME_LIST=$(ls analyses/yan_sd_2022/binning/drep98/dereplicated_genomes/*.fna)

# Step 1: Index reference if needed
conda activate bowtie2
if [[ ! -f ${REF_FASTA}.1.bt2 ]]; then
    bowtie2-build ${REF_FASTA} ${REF_FASTA}
fi

# Step 2: Map reads and sort BAM
bowtie2 -x ${REF_FASTA} -1 ${R1} -2 ${R2} -p 16 | samtools view -bS - | samtools sort -o ${BAM_FILE}
samtools index ${BAM_FILE}

# Step 3: Run inStrain profile
conda activate instrain
inStrain profile ${BAM_FILE} \
    ${REF_FASTA} \
    -o ${SAMPLE_OUTDIR}/instrain_profile \
    --stb $(echo ${GENOME_LIST}) \
    --scaffolds_to_profile ${REF_CONTIGS} \
    -p 16

# Step 4: Clean up BAM files if desired
rm ${BAM_FILE} ${BAM_FILE}.bai