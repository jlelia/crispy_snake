#!/usr/bin/env python3
"""
Count sgRNA sequences from FASTQ files by matching against library
"""

import gzip
import argparse
import sys
from pathlib import Path
from collections import defaultdict


def load_library(library_file):
    """
    Load sgRNA library file
    
    Args:
        library_file: Path to library file (tab-separated)
    
    Returns:
        Dictionary mapping sequence to (sgRNA_ID, gene)
    """
    
    library = {}
    gene_map = {}
    
    with open(library_file, 'r') as f:
        # Skip header
        header = f.readline().strip().split('\t')
        
        # Find column indices
        try:
            id_idx = header.index('sgRNA_ID')
            seq_idx = header.index('sequence')
            gene_idx = header.index('gene')
        except ValueError as e:
            print(f"Error: Required column not found in library file: {e}", file=sys.stderr)
            sys.exit(1)
        
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) > max(id_idx, seq_idx, gene_idx):
                sgrna_id = fields[id_idx]
                sequence = fields[seq_idx].upper()
                gene = fields[gene_idx]
                
                library[sequence] = sgrna_id
                gene_map[sgrna_id] = gene
    
    return library, gene_map


def open_fastq(fastq_file):
    """Open FASTQ file (handles gzipped and plain text)"""
    
    if fastq_file.endswith('.gz'):
        return gzip.open(fastq_file, 'rt')
    else:
        return open(fastq_file, 'r')


def count_sgrnas(fastq_file, library, gene_map, mismatches=0):
    """
    Count sgRNA sequences in FASTQ file
    
    Args:
        fastq_file: Path to FASTQ file
        library: Dictionary mapping sequences to sgRNA IDs
        gene_map: Dictionary mapping sgRNA IDs to genes
        mismatches: Number of allowed mismatches (0 for exact match)
    
    Returns:
        Dictionary of sgRNA counts
    """
    
    counts = defaultdict(int)
    total_reads = 0
    matched_reads = 0
    
    with open_fastq(fastq_file) as f:
        while True:
            # Read FASTQ entry (4 lines)
            header = f.readline()
            if not header:
                break
            
            sequence = f.readline().strip().upper()
            plus = f.readline()
            quality = f.readline()
            
            total_reads += 1
            
            # Try exact match first
            if sequence in library:
                sgrna_id = library[sequence]
                counts[sgrna_id] += 1
                matched_reads += 1
            elif mismatches > 0:
                # Try fuzzy matching (simple implementation)
                for lib_seq, sgrna_id in library.items():
                    if len(sequence) == len(lib_seq):
                        diff = sum(c1 != c2 for c1, c2 in zip(sequence, lib_seq))
                        if diff <= mismatches:
                            counts[sgrna_id] += 1
                            matched_reads += 1
                            break
    
    print(f"Total reads: {total_reads:,}", file=sys.stderr)
    print(f"Matched reads: {matched_reads:,} ({100*matched_reads/total_reads:.2f}%)", file=sys.stderr)
    print(f"Unique sgRNAs detected: {len(counts)}", file=sys.stderr)
    
    return counts


def write_counts(counts, gene_map, output_file):
    """
    Write count table to file
    
    Args:
        counts: Dictionary of sgRNA counts
        gene_map: Dictionary mapping sgRNA IDs to genes
        output_file: Path to output file
    """
    
    with open(output_file, 'w') as f:
        # Write header
        f.write("sgRNA_ID\tgene\tcount\n")
        
        # Write counts (sorted by sgRNA ID)
        for sgrna_id in sorted(counts.keys()):
            gene = gene_map.get(sgrna_id, "Unknown")
            count = counts[sgrna_id]
            f.write(f"{sgrna_id}\t{gene}\t{count}\n")
        
        # Also write zero counts for sgRNAs not detected
        for sgrna_id in sorted(gene_map.keys()):
            if sgrna_id not in counts:
                gene = gene_map[sgrna_id]
                f.write(f"{sgrna_id}\t{gene}\t0\n")


def main():
    """Main function"""
    
    parser = argparse.ArgumentParser(
        description="Count sgRNA sequences from FASTQ files"
    )
    parser.add_argument(
        "--fastq", required=True,
        help="Input FASTQ file (can be gzipped)"
    )
    parser.add_argument(
        "--library", required=True,
        help="sgRNA library file (tab-separated)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output count file"
    )
    parser.add_argument(
        "--mismatches", type=int, default=0,
        help="Number of allowed mismatches (default: 0)"
    )
    
    args = parser.parse_args()
    
    # Load library
    print("Loading sgRNA library...", file=sys.stderr)
    library, gene_map = load_library(args.library)
    print(f"Loaded {len(library)} sgRNAs targeting {len(set(gene_map.values()))} genes", 
          file=sys.stderr)
    
    # Count sgRNAs
    print(f"Counting sgRNAs in {args.fastq}...", file=sys.stderr)
    counts = count_sgrnas(args.fastq, library, gene_map, args.mismatches)
    
    # Write output
    print(f"Writing counts to {args.output}...", file=sys.stderr)
    write_counts(counts, gene_map, args.output)
    
    print("Done!", file=sys.stderr)


if __name__ == "__main__":
    main()
