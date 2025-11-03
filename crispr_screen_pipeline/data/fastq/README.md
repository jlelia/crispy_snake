# FASTQ Files

Place your FASTQ files in this directory.

## File Format

- Files should be gzipped FASTQ format (`.fastq.gz` or `.fq.gz`)
- Alternatively, you can specify absolute paths in `config.yaml`

## Example Files

```
control_rep1.fastq.gz
control_rep2.fastq.gz
treated_rep1.fastq.gz
treated_rep2.fastq.gz
```

## Quality Requirements

- Reads should contain sgRNA sequences
- Typical read length: 50-150 bp
- Recommended quality score: Q20+
