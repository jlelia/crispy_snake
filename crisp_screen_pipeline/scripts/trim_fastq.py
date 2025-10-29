#!/usr/bin/env python3
"""
Trim FASTQ reads to a specified length (typically 20bp for sgRNA guides)
"""

import gzip
import argparse
import sys
from pathlib import Path


def open_fastq(fastq_file, mode='r'):
    """
    Open FASTQ file (handles gzipped and plain text)
    
    Args:
        fastq_file: Path to FASTQ file
        mode: File mode ('r' for read, 'w' for write)
    
    Returns:
        File handle
    """
    if fastq_file.endswith('.gz'):
        if mode == 'r':
            return gzip.open(fastq_file, 'rt')
        else:
            return gzip.open(fastq_file, 'wt')
    else:
        return open(fastq_file, mode)


def trim_fastq(input_file, output_file, length=20):
    """
    Trim FASTQ reads to specified length
    
    Args:
        input_file: Path to input FASTQ file
        output_file: Path to output FASTQ file
        length: Number of bases to keep from the start of each read
    
    Returns:
        Tuple of (total_reads, written_reads, skipped_reads)
    """
    total_reads = 0
    written_reads = 0
    skipped_reads = 0
    
    with open_fastq(input_file, 'r') as fin, open_fastq(output_file, 'w') as fout:
        while True:
            # Read FASTQ entry (4 lines)
            header = fin.readline()
            if not header:
                break
            
            sequence = fin.readline().strip()
            plus = fin.readline()
            quality = fin.readline().strip()
            
            total_reads += 1
            
            # Skip reads that are shorter than requested length
            if len(sequence) < length:
                skipped_reads += 1
                continue
            
            # Trim sequence and quality to requested length
            trimmed_seq = sequence[:length]
            trimmed_qual = quality[:length]
            
            # Write trimmed record
            fout.write(header)
            fout.write(trimmed_seq + '\n')
            fout.write(plus)
            fout.write(trimmed_qual + '\n')
            
            written_reads += 1
    
    return total_reads, written_reads, skipped_reads


def main():
    """Main function"""
    
    parser = argparse.ArgumentParser(
        description="Trim FASTQ reads to specified length"
    )
    parser.add_argument(
        "--input", required=True,
        help="Input FASTQ file (can be gzipped)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output FASTQ file (will be gzipped if .gz extension)"
    )
    parser.add_argument(
        "--length", type=int, default=20,
        help="Number of bases to keep from start of each read (default: 20)"
    )
    
    args = parser.parse_args()
    
    # Validate length
    if args.length <= 0:
        print(f"Error: Length must be positive (got {args.length})", file=sys.stderr)
        sys.exit(1)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Trim FASTQ
    print(f"Trimming reads in {args.input} to {args.length}bp...", file=sys.stderr)
    total, written, skipped = trim_fastq(args.input, args.output, args.length)
    
    # Print statistics
    print(f"Total reads: {total:,}", file=sys.stderr)
    print(f"Written reads: {written:,} ({100*written/total:.2f}%)" if total > 0 else "Written reads: 0", file=sys.stderr)
    print(f"Skipped reads (< {args.length}bp): {skipped:,} ({100*skipped/total:.2f}%)" if total > 0 else f"Skipped reads (< {args.length}bp): 0", file=sys.stderr)
    
    print(f"Trimmed FASTQ written to {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
