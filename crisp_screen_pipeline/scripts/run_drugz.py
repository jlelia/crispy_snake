#!/usr/bin/env python3
"""
Run DrugZ analysis for CRISPR screen data
DrugZ is a method for identifying essential genes from CRISPR screens
"""

import pandas as pd
import numpy as np
from scipy import stats
import sys
from pathlib import Path


def calculate_drugz_score(treatment_counts, control_counts, pseudocount=5):
    """
    Calculate DrugZ score for comparing treatment vs control
    
    Args:
        treatment_counts: Array of treatment replicate counts
        control_counts: Array of control replicate counts
        pseudocount: Pseudocount to add to avoid log(0)
    
    Returns:
        normZ score for the comparison
    """
    
    # Add pseudocount
    treatment_counts = np.array(treatment_counts) + pseudocount
    control_counts = np.array(control_counts) + pseudocount
    
    # Calculate log2 fold changes
    log2fc_values = []
    for t_count in treatment_counts:
        for c_count in control_counts:
            log2fc = np.log2(t_count / c_count)
            log2fc_values.append(log2fc)
    
    # Calculate mean and standard deviation
    mean_log2fc = np.mean(log2fc_values)
    std_log2fc = np.std(log2fc_values, ddof=1) if len(log2fc_values) > 1 else 1.0
    
    # Avoid division by zero
    if std_log2fc == 0:
        std_log2fc = 1e-6
    
    # Calculate Z-score
    z_score = mean_log2fc / std_log2fc
    
    return mean_log2fc, std_log2fc, z_score


def normalize_counts(df, count_columns):
    """
    Normalize count columns to reads per million (RPM)
    
    Args:
        df: DataFrame with counts
        count_columns: List of column names with counts
    
    Returns:
        DataFrame with normalized counts
    """
    
    df_norm = df.copy()
    
    for col in count_columns:
        total_counts = df[col].sum()
        if total_counts > 0:
            df_norm[col] = (df[col] / total_counts) * 1e6
        else:
            df_norm[col] = 0
    
    return df_norm


def run_drugz(count_file, treatment_samples, control_samples, 
              output_file, pseudocount=5, min_observations=1):
    """
    Run DrugZ analysis
    
    Args:
        count_file: Path to merged count table
        treatment_samples: List of treatment sample names
        control_samples: List of control sample names
        output_file: Path to output results
        pseudocount: Pseudocount for log calculations
        min_observations: Minimum number of observations required
    """
    
    # Read count table
    counts_df = pd.read_csv(count_file, sep="\t")
    
    # Check if gene column exists
    if "gene" not in counts_df.columns:
        print("Error: gene column not found in count table", file=sys.stderr)
        sys.exit(1)
    
    # Verify sample columns exist
    missing_samples = []
    for sample in treatment_samples + control_samples:
        if sample not in counts_df.columns:
            missing_samples.append(sample)
    
    if missing_samples:
        print(f"Error: Missing samples in count table: {missing_samples}", file=sys.stderr)
        sys.exit(1)
    
    # Normalize counts to RPM
    all_samples = treatment_samples + control_samples
    counts_norm = normalize_counts(counts_df, all_samples)
    
    # Calculate scores per gene
    results = []
    
    for gene in counts_df["gene"].unique():
        gene_data = counts_norm[counts_norm["gene"] == gene]
        
        # Skip if not enough observations
        if len(gene_data) < min_observations:
            continue
        
        # Get treatment and control counts
        treatment_counts = gene_data[treatment_samples].values.flatten()
        control_counts = gene_data[control_samples].values.flatten()
        
        # Calculate DrugZ score
        mean_log2fc, std_log2fc, z_score = calculate_drugz_score(
            treatment_counts, control_counts, pseudocount
        )
        
        # Calculate normalized Z score (normZ)
        num_sgrnas = len(gene_data)
        normZ = z_score * np.sqrt(num_sgrnas)
        
        # Calculate p-value
        pvalue = stats.norm.sf(abs(normZ)) * 2  # two-tailed
        
        results.append({
            "gene": gene,
            "numObs": num_sgrnas,
            "log2FC": mean_log2fc,
            "stdLog2FC": std_log2fc,
            "zscore": z_score,
            "normZ": normZ,
            "pvalue": pvalue
        })
    
    # Create results dataframe
    results_df = pd.DataFrame(results)
    
    # Calculate FDR using Benjamini-Hochberg
    from scipy.stats import false_discovery_control
    if len(results_df) > 0:
        results_df["fdr"] = false_discovery_control(results_df["pvalue"].values)
    else:
        results_df["fdr"] = []
    
    # Sort by normZ score (most depleted first)
    results_df = results_df.sort_values("normZ")
    
    # Create output directory if it doesn't exist
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    
    # Write results
    results_df.to_csv(output_file, sep="\t", index=False)
    
    print(f"DrugZ analysis complete")
    print(f"Total genes analyzed: {len(results_df)}")
    print(f"Significant genes (FDR < 0.05): {sum(results_df['fdr'] < 0.05)}")
    
    return results_df


def main():
    """Main function to run when used in Snakemake"""
    
    # Get parameters from Snakemake
    if hasattr(snakemake, "input"):
        count_file = snakemake.input.counts
        output_file = snakemake.output.results
        treatment_samples = snakemake.params.treatment
        control_samples = snakemake.params.control
        pseudocount = snakemake.params.pseudocount
        min_obs = snakemake.params.min_obs
    else:
        # Standalone usage
        import argparse
        parser = argparse.ArgumentParser(description="Run DrugZ analysis")
        parser.add_argument("--counts", required=True, help="Merged count table")
        parser.add_argument("--treatment", nargs="+", required=True, help="Treatment samples")
        parser.add_argument("--control", nargs="+", required=True, help="Control samples")
        parser.add_argument("--output", required=True, help="Output file")
        parser.add_argument("--pseudocount", type=int, default=5, help="Pseudocount")
        parser.add_argument("--min-obs", type=int, default=1, help="Minimum observations")
        args = parser.parse_args()
        
        count_file = args.counts
        output_file = args.output
        treatment_samples = args.treatment
        control_samples = args.control
        pseudocount = args.pseudocount
        min_obs = args.min_obs
    
    # Run DrugZ analysis
    run_drugz(count_file, treatment_samples, control_samples,
              output_file, pseudocount, min_obs)


if __name__ == "__main__":
    main()
