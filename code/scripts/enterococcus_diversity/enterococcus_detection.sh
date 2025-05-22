#!/bin/bash

#SBATCH --array=1  # Adjust the range based on the number of libraries
#SBATCH --output=log/enterococcus_diversity/enterococcus_detection_%A_%a.out
#SBATCH --error=log/enterococcus_diversity/enterococcus_detection_%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=cpu

# Activate the conda environment with Bowtie2 and bedtools installed
source $(conda info --base)/etc/profile.d/conda.sh
conda activate bowtie2

# Define paths
READS_DIR="data/reads/yan_sd_2022"
MAGS_DIR="analyses/enterococcus_diversity/detection/ref_mags"
OUTPUT_DIR="analyses/enterococcus_diversity/detection"
COVERAGE_DIR="${OUTPUT_DIR}/coverage"
RESULTS_FILE="${OUTPUT_DIR}/detection_results.tsv"

# Get the sample name for the current task
sample_name=$(ls ${READS_DIR} | sed -n ${SLURM_ARRAY_TASK_ID}p)

# Create output directories
mkdir -p ${COVERAGE_DIR}

# Input files
R1="${READS_DIR}/${sample_name}/${sample_name}_R1_all.fastq.gz"
R2="${READS_DIR}/${sample_name}/${sample_name}_R2_all.fastq.gz"

# Initialize results file header (only for the first task)
if [[ ${SLURM_ARRAY_TASK_ID} -eq 1 ]]; then
    echo -e "sample\tmag\tcoverage (%)\tdetection" > ${RESULTS_FILE}
fi

# Iterate over all MAGs in the directory
for MAG_FILE in ${MAGS_DIR}/*.fna; do
    MAG_NAME=$(basename ${MAG_FILE} .fna)

    # Step 1: Build Bowtie2 index for the MAG (only needs to be done once per MAG)
    if [[ ! -f ${MAG_FILE}.1.bt2 ]]; then
        bowtie2-build ${MAG_FILE} ${MAG_FILE}
    fi

    # Step 2: Map reads to the MAG
    conda activate bowtie2
    bowtie2 -x ${MAG_FILE} -1 ${R1} -2 ${R2} -p 16 | samtools view -bS - | samtools sort -o ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}.sorted.bam
    samtools index ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}.sorted.bam

    # Step 3: Calculate coverage using bedtools
    conda activate bedtools
    bedtools genomecov -ibam ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}.sorted.bam -bga > ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}_coverage.bed

    # Step 4: Determine detection based on >50% coverage
    awk -v sample=${sample_name} -v mag=${MAG_NAME} '
    BEGIN {
        total_length = 0;
        covered_length = 0;
    }
    {
        if ($4 > 0) covered_length += $3 - $2;
        total_length += $3 - $2;
    }
    END {
        coverage = (covered_length / total_length) * 100;
        detection = (coverage > 50) ? "Positive" : "Negative";
        print sample "\t" mag "\t" coverage "\t" detection;
    }' ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}_coverage.bed >> ${RESULTS_FILE}

    # Clean up intermediate files to save space
    rm -f ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}.sorted.bam
    rm -f ${COVERAGE_DIR}/${sample_name}_${MAG_NAME}.sorted.bam.bai
done