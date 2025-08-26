#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import pybedtools

def main(starReport, sampleBam, outputFile):
    # --- Validate input files ---
    if not os.path.isfile(starReport):
        sys.exit(f"Error: STAR report file not found -> {starReport}")
    if not starReport.endswith("Log.final.out"):
        sys.exit("Error: STAR report file must be Log.final.out")
    if not os.path.isfile(sampleBam):
        sys.exit(f"Error: Sample BAM file not found -> {sampleBam}")

    # get the library size from the STAR report
    starReport = pd.read_csv(starReport, sep="\t", header=None)
    lib_size = int(starReport.loc[4,1])

    # --- Load BAM file with pybedtools and convert to df---
    df = pybedtools.BedTool(sampleBam).bam_to_bed().to_dataframe()

    df["postive"] = [1 if x == "+" else 0 for x in df["strand"] ]
    df["negative"] = [1 if x == "-" else 0 for x in df["strand"] ]

    # Create bins for the "start" column
    bins = range(0, df["start"].max() + 25, 25)
    df["start_bin"] = pd.cut(df["start"], bins=bins, right=False)

    # Group by the bins and sum the "positive" and "negative" columns
    grouped_df = df.groupby("start_bin",observed=False)[["postive", "negative"]].sum()

    grouped_df["pos"] = grouped_df["postive"]/lib_size*1000000
    grouped_df["neg"] = -grouped_df["negative"]/lib_size*1000000

    grouped_df["Position"] = [interval.left + 1 for interval in grouped_df.index.categories]
    grouped_df["Chromosome"] = "L1ME1_dfam"
    grouped_df[["18_23.pos","18_23.neg","24_35.pos","24_35.neg"]] = 0

    grouped_df[['Chromosome','Position', '18_23.pos', '18_23.neg', '24_35.pos', '24_35.neg', 'pos', 'neg']].to_csv(outputFile, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""
        This tool converts a BAM file into a TSV file for generating TE coverage plots.
        Input: A single-end BAM file aligned to a TE consensus or extracted from a specific region. If the BAM file is paired-end, please first separate forward and reverse reads and then merged into a single BAM file before using this tool.
        Processing: Read counts are normalized to RPM (Reads Per Million) using the library size reported by STAR.
        Output: A TSV file compatible with https://bitbucket.org/bucab/retrotransposon_smrnaseq_viz/src/master/ for TE coverage visualization.
        
        Usage:
          fromBam2cpBed.py <star_report> <sample.bam> <output.file>

        Example:
          fromBam2cpBed.py Log.final.out sample.bam TE@BC_Y1_L1ME1
        """
    )
    parser.add_argument("starReport", help="Star alignment report file (Log.final.out)")
    parser.add_argument("sampleBam", help="a specific region or a TE consensus BAM file")
    parser.add_argument("outputFile")

    args = parser.parse_args()
    main(args.starReport, args.sampleBam, args.outputFile)


