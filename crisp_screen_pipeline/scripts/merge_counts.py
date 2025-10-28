#!/usr/bin/env python3
"""
Merge count tables from multiple samples into a single matrix
"""

import pandas as pd
import sys
from pathlib import Path


def merge_count_tables(count_files, output_file):
    """
    Merge multiple count files into a single count matrix
    
    Args:
        count_files: List of paths to individual count files
        output_file: Path to output merged count file
    """
    
    # Read all count files
    count_dfs = {}
    
    for count_file in count_files:
        # Extract sample name from filename
        sample_name = Path(count_file).stem.replace("_counts", "")
        
        # Read count file
        df = pd.read_csv(count_file, sep="\t")
        
        # Check if required columns exist
        if "sgRNA_ID" not in df.columns or "count" not in df.columns:
            print(f"Warning: {count_file} missing required columns", file=sys.stderr)
            continue
            
        # Store counts for this sample
        count_dfs[sample_name] = df.set_index("sgRNA_ID")["count"]
    
    # Merge all count dataframes
    merged_df = pd.DataFrame(count_dfs)
    
    # Add gene column if available from first file
    first_file = pd.read_csv(count_files[0], sep="\t")
    if "gene" in first_file.columns:
        gene_map = first_file.set_index("sgRNA_ID")["gene"]
        merged_df.insert(0, "gene", merged_df.index.map(gene_map))
    
    # Reset index to make sgRNA_ID a column
    merged_df = merged_df.reset_index()
    merged_df = merged_df.rename(columns={"index": "sgRNA_ID"})
    
    # Sort by gene and sgRNA_ID
    if "gene" in merged_df.columns:
        merged_df = merged_df.sort_values(["gene", "sgRNA_ID"])
    else:
        merged_df = merged_df.sort_values("sgRNA_ID")
    
    # Fill NaN values with 0
    merged_df = merged_df.fillna(0)
    
    # Convert counts to integers
    count_cols = [col for col in merged_df.columns if col not in ["sgRNA_ID", "gene"]]
    merged_df[count_cols] = merged_df[count_cols].astype(int)
    
    # Write merged counts
    merged_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"Merged {len(count_files)} count files")
    print(f"Total sgRNAs: {len(merged_df)}")
    print(f"Samples: {', '.join(count_cols)}")
    
    return merged_df


def main():
    """Main function to run when used in Snakemake"""
    
    # Get input and output from Snakemake
    if hasattr(snakemake, "input"):
        count_files = snakemake.input.counts
        output_file = snakemake.output.merged
    else:
        # Standalone usage
        import argparse
        parser = argparse.ArgumentParser(description="Merge CRISPR count tables")
        parser.add_argument("--inputs", nargs="+", required=True, help="Input count files")
        parser.add_argument("--output", required=True, help="Output merged count file")
        args = parser.parse_args()
        
        count_files = args.inputs
        output_file = args.output
    
    # Merge count tables
    merge_count_tables(count_files, output_file)


if __name__ == "__main__":
    main()
