#!/bin/bash
#SBATCH --output=log/bacteroides_pul/get_pulabund_%A_%a.out
#SBATCH --error=log/bacteroides_pul/get_pulabund_%A_%a.err
#SBATCH --array=1-2053%40
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH -p cpu
#SBATCH --time=24:00:00

# Activate conda
source $(conda info --base)/etc/profile.d/conda.sh

# Set sample and read paths
# Remove brackets and single quotes: e.g., ['275C'] -> 275C
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" misc_files/bacteroides_pul/isabl1_succesful_assemblies.txt)
# Get sample name (with brackets and quotes) for current task
sample_bracketed="['${SAMPLE}']"
# Get first record for this sample from the paths table
record=$(awk -v s="$sample_bracketed" '$2==s{print; exit}' misc_files/enterococcus_diversity/all_target_samples_paths.txt)
# Extract R1 and R2 paths (columns 3 and 4)
R1=$(echo "$record" | cut -f3)
R2=$(echo "$record" | cut -f4)
# Replace storage path
R1=$(echo "$R1" | sed 's|/data/brinkvd/|/data1/collab004/|')
R2=$(echo "$R2" | sed 's|/data/brinkvd/|/data1/collab004/|')

# Set directories
ASSEMBLY_DIR=analyses/bacteroides_pul/metagenomes/assembly/metaspades
SAM_DIR="analyses/bacteroides_pul/pul_prediction/samfiles"
ABUND_DIR="analyses/bacteroides_pul/pul_prediction/abund"
SAMPLE_ABUND_DIR="$ABUND_DIR/$SAMPLE"
DBCAN_OUT="analyses/bacteroides_pul/pul_prediction/dbcan/${SAMPLE}"
mkdir -p "$SAM_DIR" "$ABUND_DIR" "$SAMPLE_ABUND_DIR"

# Contigs and FFN files
CONTIGS="$ASSEMBLY_DIR/$SAMPLE/contigs_1000.fasta"
#FFN="$DBCAN_OUT/$SAMPLE.ffn"

# 1. Index contigs
conda activate bwa
bwa index "$CONTIGS"

# 2. Map reads to CDS and contigs
bwa mem -t 8 -o "$SAM_DIR/${SAMPLE}.sam" "$CONTIGS" "${R1}" "${R2}"

# 3. Sort SAM files and convert to indexed BAM
conda activate bowtie2
samtools sort -@ 8 -o "$SAM_DIR/${SAMPLE}.bam" "$SAM_DIR/${SAMPLE}.sam"
samtools index "$SAM_DIR/${SAMPLE}.bam"

# 4. Calculate read coverage for all proteins using dbcan_utils
conda activate dbcan
dbcan_utils cal_coverage \
    -g "$DBCAN_OUT/uniInput.gff" \
    -i "$SAM_DIR/${SAMPLE}.bam" \
    -o "$SAMPLE_ABUND_DIR/${SAMPLE}.depth.txt" \
    -t 8

# 5. Remove SAM and BAM files to save space
rm "$SAM_DIR/${SAMPLE}.sam" "$SAM_DIR/${SAMPLE}.bam" "$SAM_DIR/${SAMPLE}.bam.bai"

# 6. Estimate abundance of CAZymes, CGCs, and PULs using dbcan_utils
cd "$SAMPLE_ABUND_DIR"
dbcan_utils fam_abund \
    -bt "${SAMPLE}.depth.txt" \
    -i "../../../../../$DBCAN_OUT" \
    -o "${SAMPLE}.fam_abund.txt" \
    -a TPM
dbcan_utils fam_substrate_abund \
    -bt "${SAMPLE}.depth.txt" \
    -i "../../../../../$DBCAN_OUT" \
    -o "$SAMPLE_ABUND_DIR/${SAMPLE}.fam_substrate_abund.txt" \
    -a TPM
dbcan_utils CGC_abund \
    -bt "${SAMPLE}.depth.txt" \
    -i "../../../../../$DBCAN_OUT" \
    --output "$SAMPLE_ABUND_DIR/" \
    -a TPM
dbcan_utils CGC_substrate_abund \
    -bt "${SAMPLE}.depth.txt" \
    -i "../../../../../$DBCAN_OUT" \
    -o "$SAMPLE_ABUND_DIR/${SAMPLE}.CGC_substrate_abund.txt" \
    -a TPM