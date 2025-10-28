#!/usr/bin/env python
import gzip
from Bio import SeqIO
import pandas as pd
import sys

def phred33_to_qscore(char):
    return ord(char) - 33

def extract_fastplong_style(fq_file):
    if fq_file.endswith(".gz"):
        handle = gzip.open(fq_file, "rt")
    else:
        handle = open(fq_file, "r")
    
    for rec in SeqIO.parse(handle, "fastq"):
        qscores = rec.letter_annotations["phred_quality"]
        mean_q = sum(qscores) / len(qscores)
        yield {
            "id": rec.id,
            "mean_qual": mean_q,
            "length": len(rec)
        }
    
    handle.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python fastplong_style.py <fastq_file>")
        sys.exit(1)
    
    fastq_file = sys.argv[1]
    
    records = list(extract_fastplong_style(fastq_file))
    df = pd.DataFrame(records)
    
    print(df)
    df.to_csv(f"{fastq_file}.fastplong_stats.tsv", sep="\t", index=False)
    
    print("\nProcessed {} reads".format(len(df)))
    print("Mean quality over all reads: {:.2f}".format(df["mean_qual"].mean()))

