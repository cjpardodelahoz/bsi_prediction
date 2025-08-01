#!/usr/bin/env python3

import argparse
import os
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed

def extract_from_one_gbk(args):
    gbk_file, genes = args
    basename = os.path.splitext(os.path.basename(gbk_file))[0]
    gene_seqs = {gene: [] for gene in genes}
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                gene_name = feature.qualifiers.get("gene", [""])[0]
                locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                for gene in genes:
                    if gene == gene_name or gene == locus_tag:
                        aa_seq = feature.qualifiers.get("translation", [""])[0]
                        if aa_seq:
                            gene_seqs[gene].append((basename, aa_seq))
    return gene_seqs

def extract_gene_aa_parallel(gbk_files, genes, outdir, nproc=4):
    os.makedirs(outdir, exist_ok=True)
    all_gene_seqs = {gene: [] for gene in genes}
    with ProcessPoolExecutor(max_workers=nproc) as executor:
        futures = [executor.submit(extract_from_one_gbk, (gbk_file, genes)) for gbk_file in gbk_files]
        for future in as_completed(futures):
            gene_seqs = future.result()
            for gene in genes:
                all_gene_seqs[gene].extend(gene_seqs[gene])
    # Write output
    for gene, seqs in all_gene_seqs.items():
        with open(os.path.join(outdir, f"{gene}.faa"), "w") as f:
            for basename, aa_seq in seqs:
                f.write(f">{basename}\n{aa_seq}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract amino acid sequences of specific genes from one or more .gbk files.\n"
                    "Writes one output FASTA per gene, with headers as the basename of the source file.\n"
                    "Requires Biopython (tested with version 1.78).",
        epilog="Example:\n"
               "  python extract_gene_aa.py file1.gbk file2.gbk -g vanA vanB -o output_dir\n"
               "This will create vanA.faa and vanB.faa in output_dir with sequences from all input files.",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("gbk_files", nargs="+", help="Input .gbk files (one or more)")
    parser.add_argument("-g", "--genes", nargs="+", required=True, help="Gene name(s) or locus_tag(s) to extract")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: current directory)")
    parser.add_argument("-p", "--nproc", type=int, default=4, help="Number of parallel processes (default: 4)")
    args = parser.parse_args()
    extract_gene_aa_parallel(args.gbk_files, args.genes, args.outdir, args.nproc)