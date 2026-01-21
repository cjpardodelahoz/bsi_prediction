import os
import argparse
import pandas as pd

# Mapping of output types to file names and column labels
OUTPUT_TYPES = {
    "fam_substrate": {
        "file": "fam_substrate_abund.out",
        "feature_col": "Substrate",
        "abund_col": "Abundance"
    },
    "fam": {
        "file": "fam_abund.out",
        "feature_col": "Family",
        "abund_col": "Abundance"
    },
    "subfam": {
        "file": "subfam_abund.out",
        "feature_col": "Subfamily",
        "abund_col": "Abundance"
    },
    "PUL": {
        "file": "CGC_substrate_PUL_homology.out",
        "feature_col": "Substrate",
        "abund_col": "Abundance(sum of CGC)"
    },
    "EC": {
        "file": "EC_abund.out",
        "feature_col": "EC",
        "abund_col": "Abundance"
    }
}

def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a feature table from DBcan quantification output files for specified samples."
    )
    parser.add_argument(
        "--samples",
        required=True,
        help="Path to a text file with sample names (one per line)."
    )
    parser.add_argument(
        "--output",
        default="feature_table.tsv",
        help="Output TSV file for the feature table (default: feature_table.tsv)."
    )
    parser.add_argument(
        "--base_dir",
        default="analyses/bacteroides_pul/pul_prediction/abund/",
        help="Base directory containing sample subdirectories (default: analyses/bacteroides_pul/pul_prediction/abund/)."
    )
    parser.add_argument(
        "--output_type",
        required=True,
        choices=OUTPUT_TYPES.keys(),
        help="Type of DBcan output to process. Choices: " + ", ".join(OUTPUT_TYPES.keys())
    )
    return parser.parse_args()

def main():
    args = parse_args()
    output_info = OUTPUT_TYPES[args.output_type]
    file_name = output_info["file"]
    feature_col = output_info["feature_col"]
    abund_col = output_info["abund_col"]

    with open(args.samples) as f:
        sample_names = [line.strip() for line in f if line.strip()]

    feature_dfs = []
    for sample in sample_names:
        file_path = os.path.join(args.base_dir, sample, file_name)
        if not os.path.isfile(file_path):
            print(f"Warning: {file_path} not found, skipping sample {sample}")
            continue
        try:
            df = pd.read_csv(file_path, sep="\t", usecols=[feature_col, abund_col], header=0)
        except ValueError:
            print(f"Warning: Could not find columns {feature_col} and {abund_col} in {file_path}, skipping sample {sample}")
            continue
        df = df.set_index(feature_col).T
        df.index = [sample]
        feature_dfs.append(df)

    if not feature_dfs:
        print("No valid sample files found. Exiting.")
        return

    feature_table = pd.concat(feature_dfs, sort=True).fillna(0)
    feature_table = feature_table.reset_index().rename(columns={"index": "sample"})
    feature_table.to_csv(args.output, sep="\t", index=False)
    print(f"Feature table saved as {args.output}")

if __name__ == "__main__":
    main()