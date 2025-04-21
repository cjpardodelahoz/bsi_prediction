import json
import pandas as pd
import argparse
import sys

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract genome metadata from an NCBI JSONL file and save it as a CSV file."
    )
    parser.add_argument(
        "input_file",
        type=str,
        help="Path to the input NCBI JSONL file containing genome metadata."
    )
    parser.add_argument(
        "output_file",
        type=str,
        help="Path to the output CSV file where extracted metadata will be saved."
    )
    return parser.parse_args()

def process_jsonl(input_file, output_file):
    """Process the JSONL file and save metadata to a CSV file."""
    # List to store extracted data
    data = []

    try:
        # Open the JSONL file and process each line
        with open(input_file, "r") as f:
            for line in f:
                # Parse the JSON object
                entry = json.loads(line)
                
                # Extract relevant fields
                data.append({
                    "accession": entry.get("accession"),
                    "taxid": entry.get("organism", {}).get("taxId"),
                    "organism_name": entry.get("organism", {}).get("organismName"),
                    "species": entry.get("averageNucleotideIdentity", {}).get("submittedSpecies"),
                    "strain": entry.get("organism", {}).get("infraspecificNames", {}).get("strain"),
                    "typelabel": entry.get("typeMaterial", {}).get("typeLabel"),
                    "assembly_level": entry.get("assemblyInfo", {}).get("assemblyLevel"),
                    "assembly_name": entry.get("assemblyInfo", {}).get("assemblyName"),
                    "bioproject": entry.get("assemblyInfo", {}).get("bioprojectAccession"),
                    "biosample": entry.get("assemblyInfo", {}).get("biosample", {}).get("accession"),
                    "refseq_category": entry.get("assemblyInfo", {}).get("refseqCategory"),
                    "total_sequence_length": entry.get("assemblyStats", {}).get("totalSequenceLength"),
                    "gc_percent": entry.get("assemblyStats", {}).get("gcPercent"),
                    "num_contigs": entry.get("assemblyStats", {}).get("numberOfContigs"),
                    "submitter": entry.get("assemblyInfo", {}).get("submitter"),
                    "release_date": entry.get("assemblyInfo", {}).get("releaseDate"),
                    "source_database": entry.get("sourceDatabase"),
                    "checkm_completeness": entry.get("checkmInfo", {}).get("completeness"),
                    "checkm_contamination": entry.get("checkmInfo", {}).get("contamination"),
                    "best_ani": entry.get("averageNucleotideIdentity", {}).get("bestAniMatch", {}).get("ani"),
                    "best_ani_organism": entry.get("averageNucleotideIdentity", {}).get("bestAniMatch", {}).get("organismName"),
                    "best_ani_accession": entry.get("averageNucleotideIdentity", {}).get("bestAniMatch", {}).get("assembly"),
                    "ani_match_status": entry.get("averageNucleotideIdentity", {}).get("matchStatus"),
                    "habitat": next(
                        (attr.get("value") for attr in entry.get("assemblyInfo", {}).get("biosample", {}).get("attributes", [])
                         if attr.get("name") == "habitat"), None
                    ),
                    "host_taxid": next(
                        (attr.get("value") for attr in entry.get("assemblyInfo", {}).get("biosample", {}).get("attributes", [])
                         if attr.get("name") == "host_taxid"), None
                    )
                })
    except FileNotFoundError:
        print(f"Error: The file '{input_file}' was not found.")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"Error: Failed to parse JSON in the file '{input_file}'.")
        print(f"Details: {e}")
        sys.exit(1)

    # Convert the list of dictionaries to a Pandas DataFrame
    df = pd.DataFrame(data)

    # Save the DataFrame to a CSV file
    try:
        df.to_csv(output_file, index=False)
        print(f"Metadata has been successfully saved to {output_file}")
    except Exception as e:
        print(f"Error: Failed to save the output file '{output_file}'.")
        print(f"Details: {e}")
        sys.exit(1)

if __name__ == "__main__":
    # Parse command-line arguments
    args = parse_arguments()

    # Process the JSONL file and save metadata to a CSV file
    process_jsonl(args.input_file, args.output_file)