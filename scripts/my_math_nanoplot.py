#!/usr/bin/env python
import gzip
from Bio import SeqIO
import pandas as pd
from math import log
import sys

# Precompute error probabilities for PHRED 0–128
def errs_tab(n):
    return [10 ** (-q / 10) for q in range(n + 1)]

# Calculate NanoPlot-style average quality for a read
def ave_qual(quals, tab=errs_tab(128)):
    if quals:
        return -10 * log(sum([tab[q] for q in quals]) / len(quals), 10)
    else:
        return None

# Extract per-read info: ID, mean quality, length
def extract_from_fastq(fq_file):
    if fq_file.endswith(".gz"):
        handle = gzip.open(fq_file, "rt")
    else:
        handle = open(fq_file, "r")
    
    for rec in SeqIO.parse(handle, "fastq"):
        yield {
            "id": rec.id,
            "mean_qual": ave_qual(rec.letter_annotations["phred_quality"]),
            "length": len(rec)
        }
    handle.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fastq_quality_table.py <fastq_file>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    
    # Collect all reads
    records = list(extract_from_fastq(fastq_file))
    
    # Convert to DataFrame
    df = pd.DataFrame(records)
    
    # Print table
    print(df)
    
    # Optional: write to CSV
    df.to_csv(f"{fastq_file}.read_stats.tsv", sep="\t", index=False)
    
    # Print overall mean quality
    print("\nProcessed {} reads".format(len(df)))
    print("Mean quality over all reads: {:.2f}".format(df["mean_qual"].mean()))

