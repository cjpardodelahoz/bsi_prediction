#!/bin/bash

# Command-line application to search sequences against multiple HMM profiles and output e-values in TSV format

# Function to display the help menu
function show_help() {
    echo "Usage: $0 -m <HMM_FILES> -i <INPUT_FASTA> -o <OUTPUT_TSV> -d <DATA_TYPE> [-c <CPUS>]"
    echo
    echo "Options:"
    echo "  -m, --hmms            Comma-separated list of HMM profile files (e.g., profile1.hmm,profile2.hmm)"
    echo "  -i, --input           Input FASTA file containing sequences to search"
    echo "  -o, --output          Output TSV file with e-values (rows=sequences, columns=HMMs)"
    echo "  -d, --datatype        Type of sequence data: 'dna' or 'aa' (protein)"
    echo "  -c, --cpus            Number of CPU cores to use (optional, default: 4)"
    echo "  -h, --help            Display this help menu"
    echo
    echo "Description:"
    echo "  This script searches sequences against multiple HMM profiles using HMMER."
    echo "  For DNA sequences, it uses 'nhmmer'. For protein sequences, it uses 'hmmsearch'."
    echo "  Output is a TSV with one row per sequence and one column per HMM profile."
    echo "  E-values are reported for hits; missing hits are filled with value 1."
    echo
    echo "Requirements:"
    echo "  - nhmmer and hmmsearch (from HMMER suite) must be in PATH"
    echo "  - seqkit must be in PATH"
}

# Parse command-line arguments
CPUS=4  # Default number of CPU cores
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -m|--hmms) HMM_FILES="$2"; shift ;;
        -i|--input) INPUT_FASTA="$2"; shift ;;
        -o|--output) OUTPUT_TSV="$2"; shift ;;
        -d|--datatype) DATA_TYPE="$2"; shift ;;
        -c|--cpus) CPUS="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown option: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Check if all required arguments are provided
if [[ -z "$HMM_FILES" || -z "$INPUT_FASTA" || -z "$OUTPUT_TSV" || -z "$DATA_TYPE" ]]; then
    echo "Error: Missing required arguments."
    show_help
    exit 1
fi

# Validate data type
DATA_TYPE=$(echo "$DATA_TYPE" | tr '[:upper:]' '[:lower:]')
if [[ "$DATA_TYPE" != "dna" && "$DATA_TYPE" != "aa" ]]; then
    echo "Error: Data type must be 'dna' or 'aa'."
    exit 1
fi

# Check if required binaries are in the PATH
if [[ "$DATA_TYPE" == "dna" ]]; then
    if ! command -v nhmmer &> /dev/null; then
        echo "Error: 'nhmmer' is not installed or not in the PATH."
        exit 1
    fi
    HMMER_CMD="nhmmer"
else
    if ! command -v hmmsearch &> /dev/null; then
        echo "Error: 'hmmsearch' is not installed or not in the PATH."
        exit 1
    fi
    HMMER_CMD="hmmsearch"
fi

if ! command -v seqkit &> /dev/null; then
    echo "Error: 'seqkit' is not installed or not in the PATH."
    exit 1
fi

# Check if input FASTA exists
if [[ ! -f "$INPUT_FASTA" ]]; then
    echo "Error: Input FASTA file '$INPUT_FASTA' not found."
    exit 1
fi

# Convert comma-separated HMM files to array
IFS=',' read -ra HMM_ARRAY <<< "$HMM_FILES"

# Validate that all HMM files exist
for hmm_file in "${HMM_ARRAY[@]}"; do
    if [[ ! -f "$hmm_file" ]]; then
        echo "Error: HMM file '$hmm_file' not found."
        exit 1
    fi
done

# Create temporary directory for intermediate files
TEMP_DIR=$(mktemp -d)
trap "rm -rf $TEMP_DIR" EXIT

echo "Extracting sequence headers from input FASTA..."
# Extract all sequence headers from input FASTA (only first field before space)
seqkit seq -n "$INPUT_FASTA" | awk '{print $1}' > "$TEMP_DIR/all_sequences.txt"

# Initialize associative array to store e-values
# We'll use a file-based approach for simplicity
echo "Running HMMER searches..."

# Process each HMM file
hmm_index=0
declare -a HMM_NAMES
for hmm_file in "${HMM_ARRAY[@]}"; do
    # Get basename of HMM file for column name
    hmm_basename=$(basename "$hmm_file" .hmm)
    HMM_NAMES[$hmm_index]="$hmm_basename"
    
    echo "  Searching with $hmm_basename..."
    
    # Run HMMER search
    tblout_file="$TEMP_DIR/${hmm_basename}.tblout"
    
    if [[ "$DATA_TYPE" == "dna" ]]; then
        nhmmer --tblout "$tblout_file" --cpu "$CPUS" "$hmm_file" "$INPUT_FASTA" > /dev/null 2>&1
    else
        hmmsearch --tblout "$tblout_file" --cpu "$CPUS" "$hmm_file" "$INPUT_FASTA" > /dev/null 2>&1
    fi
    
    # Parse tblout file and extract sequence IDs and e-values
    # For each sequence, keep only the best (lowest) e-value if multiple hits
    grep -v "^#" "$tblout_file" | awk '{print $1, $5}' | \
    awk '{if (!seen[$1] || $2 < best[$1]) {best[$1] = $2; seen[$1] = 1}} END {for (seq in best) print seq, best[seq]}' \
    > "$TEMP_DIR/${hmm_basename}.evalues"
    
    hmm_index=$((hmm_index + 1))
done

echo "Compiling results into TSV format..."

# Write TSV header
{
    echo -n "sequence"
    for hmm_name in "${HMM_NAMES[@]}"; do
        echo -n $'\t'"$hmm_name"
    done
    echo
} > "$OUTPUT_TSV"

# For each sequence, compile e-values from all HMM searches
while read -r seq_id; do
    echo -n "$seq_id" >> "$OUTPUT_TSV"
    
    for hmm_name in "${HMM_NAMES[@]}"; do
        # Look up e-value for this sequence in this HMM's results
        evalue=$(grep "^${seq_id} " "$TEMP_DIR/${hmm_name}.evalues" | awk '{print $2}')
        
        # If no hit found, use 1 as default
        if [[ -z "$evalue" ]]; then
            evalue=1
        fi
        
        echo -n $'\t'"$evalue" >> "$OUTPUT_TSV"
    done
    
    echo >> "$OUTPUT_TSV"
done < "$TEMP_DIR/all_sequences.txt"

echo "Results saved to $OUTPUT_TSV"
echo "Total sequences: $(wc -l < "$TEMP_DIR/all_sequences.txt")"
echo "Total HMM profiles: ${#HMM_NAMES[@]}"
