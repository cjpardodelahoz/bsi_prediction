#!/bin/bash

#SBATCH --array=1-395%10  # Adjust the range (1-395) based on the number of samples
#SBATCH --output=log/bacteroides_variation/bacteroides_comebin_single_%A_%a.out
#SBATCH --error=log/bacteroides_variation/bacteroides_comebin_single_%A_%a.err
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --partition=cpu

# Activate COMEBin environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate comebin_env

# Define paths
ASSEMBLY_DIR="analyses/bacteroides_variation/genomes/kraken_binning/bacteroides_metaspades"
READS_DIR="data/reads/yan_sd_2022"
OUTPUT_DIR="analyses/bacteroides_variation/genomes/kraken_binning/comebin/sinlge"

# Get the sample name for the current task
sample_name=$(ls ${ASSEMBLY_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Define paths for the current sample
ASSEMBLY="${ASSEMBLY_DIR}/${sample_name}/contigs.fasta"
R1="${READS_DIR}/${sample_name}/${sample_name}_R1_all.fastq.gz"
R2="${READS_DIR}/${sample_name}/${sample_name}_R2_all.fastq.gz"
SAMPLE_BAM_DIR="${OUTPUT_DIR}/${sample_name}/bamfiles"

# Create output directories
mkdir -p ${SAMPLE_BAM_DIR}

# Step 1: Map reads to the assembly using Bowtie2
bowtie2-build ${ASSEMBLY} ${ASSEMBLY}_index
bowtie2 -x ${ASSEMBLY}_index -1 ${R1} -2 ${R2} -S ${SAMPLE_BAM_DIR}/${sample_name}.sam -p 16
samtools view -bS ${SAMPLE_BAM_DIR}/${sample_name}.sam > ${SAMPLE_BAM_DIR}/${sample_name}.bam
samtools sort ${SAMPLE_BAM_DIR}/${sample_name}.bam -o ${SAMPLE_BAM_DIR}/${sample_name}_sorted.bam
samtools index ${SAMPLE_BAM_DIR}/${sample_name}_sorted.bam
rm ${SAMPLE_BAM_DIR}/${sample_name}.sam ${SAMPLE_BAM_DIR}/${sample_name}.bam  # Remove intermediate files

# Step 2: Run COMEBin
bash run_comebin.sh -a ${ASSEMBLY} \
    -o ${OUTPUT_DIR}/${sample_name} \
    -p ${SAMPLE_BAM_DIR} \
    -t 48

# Step 3: Clean up BAM files to save space
rm -r ${SAMPLE_BAM_DIR}

echo "COMEBin binning completed for sample ${sample_name}!"