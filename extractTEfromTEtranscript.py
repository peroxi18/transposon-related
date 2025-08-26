#!/usr/bin/env python3
import argparse
import os
import sys
import gzip
import pandas as pd
import pysam
import pybedtools


def qualities_to_str(read):
    """Convert Phred quality scores to FASTQ ASCII string."""
    if read.query_qualities is None:
        return None
    return "".join(chr(q + 33) for q in read.query_qualities)


def main(teRef, sampleBam, outputFastq):
    # --- Validate input files ---
    if not os.path.isfile(teRef):
        sys.exit(f"Error: TE reference file not found -> {teRef}")
    if not os.path.isfile(sampleBam):
        sys.exit(f"Error: Sample BAM file not found -> {sampleBam}")
    if not outputFastq.endswith(".fastq.gz"):
        sys.exit("Error: Output file must end with .fastq.gz")

    # --- Load BAM file with pysam ---
    bamfile = pysam.AlignmentFile(sampleBam, "rb") # bam should be seperated into two ends

    # --- Convert BAM to BED and intersect with TE reference ---
    teRef = pybedtools.BedTool(teRef) 
    sampleBed = pybedtools.BedTool(sampleBam).bam_to_bed(split=True)
    sampleTE = sampleBed.intersect(teRef, wa=False, wb=False, wo=True, split=True)

    # --- Convert to DataFrame and clean ---
    df = sampleTE.to_dataframe()
    df.columns = ["rchr", "rpos", "rend", "rname", "score", "rstrand", "sub_chr", "sub_start", "sub_end", "sub_length"]
    # remove read-pair suffix like /1 /2
    df["rname"] = df["rname"].str.replace(r"/[12]$", "", regex=True)

    # --- Keep the longest overlapping substring per read ---
    df = df.loc[df.groupby(["rname", "sub_chr", "sub_start", "sub_end"])["sub_length"].idxmax()].reset_index(drop=True)

    # --- Extract substrings ---
    seen = set()
    seq_l = []
    for row in df.itertuples(index=False):
        for read in bamfile.fetch(row.rchr, row.rpos, row.rend):
            if read.query_name == row.rname:
                key = (row.rname, row.rchr, row.rpos, row.sub_chr, row.sub_start)
                if key in seen:
                    continue
                seen.add(key)

                substring_start = max(row.sub_start, read.pos)
                substring_end = min(row.sub_end, read.reference_end)

                seq = read.query_sequence[substring_start - read.pos:row.sub_length + substring_start - read.pos]
                qual = qualities_to_str(read)[substring_start - read.pos:row.sub_length + substring_start - read.pos]

                seq_l.append({
                    "seq_name": f"{row.rname}_{row.sub_chr}:{substring_start}-{substring_end}",
                    "seq": seq,
                    "qual": qual
                })

    # --- Write to FASTQ ---
    with gzip.open(outputFastq, "wt") as f:
        for row in seq_l:
            f.write(f"@{row['seq_name']}\n{row['seq']}\n+\n{row['qual']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        Extract TE-overlapping substrings from a BAM file and write them to a compressed FASTQ.
        1. Bam should be seperated into two ends then merged into one bam file.
        2. TE reference should be in BED format. (0-based coordinates)
        3. Sample BAM file must be indexed.
        4. module load python3/3.10.5
        5. module load bedtools

        Usage:
          extractTEfromTEtranscript.py <TE_reference.bed> <sample.bam> <output.fastq.gz>

        Example:
          extractTEfromTEtranscript.py TE_ref.bed sample.bam output.fastq.gz
        """
    )
    parser.add_argument("teRef", help="TE reference file in BED format")
    parser.add_argument("sampleBam", help="Sample BAM file (must be indexed)")
    parser.add_argument("outputFastq", help="Output FASTQ file (.fastq.gz)")

    args = parser.parse_args()
    main(args.teRef, args.sampleBam, args.outputFastq)
