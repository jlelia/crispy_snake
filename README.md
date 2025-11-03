# CRISPR Screen Analysis Pipeline

A comprehensive Snakemake pipeline for analyzing CRISPR knockout screening data. This pipeline processes raw FASTQ files through quality control, sgRNA counting, and statistical analysis using MAGeCK and DrugZ.

## Features

- **Quality Control**: FastQC analysis of raw sequencing reads
- **sgRNA Counting**: Extraction and quantification of sgRNA sequences
- **Statistical Analysis**: 
  - MAGeCK for robust identification of essential genes
  - DrugZ for sensitive detection of gene depletion/enrichment
- **Automated Reporting**: HTML reports with analysis summaries
- **Reproducible**: Conda environments ensure consistent software versions

## Directory Structure

```
crispr_screen_pipeline/
├── config.yaml              # Main configuration file
├── Snakefile               # Workflow definition
├── envs/                   # Conda environment specifications
│   ├── env_qc.yaml        # FastQC environment
│   ├── env_mageck.yaml    # MAGeCK environment
│   └── env_python.yaml    # Python analysis environment
├── scripts/                # Analysis scripts
│   ├── merge_counts.py    # Merge count tables
│   └── run_drugz.py       # DrugZ analysis
├── data/                   # Input data directory
│   ├── fastq/             # FASTQ files
│   └── library/           # sgRNA library files
└── results/                # Output directory
    ├── qc/                # Quality control reports
    ├── counts/            # Count tables
    ├── mageck/            # MAGeCK results
    ├── drugz/             # DrugZ results
    └── reports/           # Summary reports
```

## Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) (>= 7.0)
- [Conda](https://docs.conda.io/en/latest/) or [Mamba](https://mamba.readthedocs.io/)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/crispy_snake.git
cd crispy_snake/crispr_screen_pipeline
```

2. Install Snakemake (if not already installed):
```bash
conda install -c conda-forge -c bioconda snakemake
```

## Naming Conventions (Required)

**All names in the configuration file must follow strict naming rules to ensure compatibility with analysis tools.**

### Allowed Characters

Names (sample names, comparison names, treatment/control labels) may only contain:
- Letters (A-Z, a-z)
- Numbers (0-9)
- Dots (.)
- Underscores (_)
- Dashes (-)

**Spaces and other special characters are NOT allowed.**

### Valid Examples
```yaml
✓ control_rep1
✓ treated_rep2
✓ MMC_vs_DMSO
✓ treatment-vs-control
✓ exp.2024.batch1
✓ sample_A1
```

### Invalid Examples
```yaml
✗ control rep1          # contains space
✗ treated(rep2)         # contains parentheses
✗ MMC vs DMSO           # contains spaces
✗ treatment/control     # contains forward slash
✗ sample@A1             # contains @ symbol
```

### Why These Rules?

The pipeline uses various bioinformatics tools (like MAGeCK RRA) that parse filenames and can break when encountering spaces or special characters. By enforcing these naming conventions early, we prevent pipeline failures and ensure reliable analysis.

**The pipeline will validate all names at startup and fail immediately with a clear error message if any invalid characters are detected.**

## Configuration

Edit `config.yaml` to specify your experimental setup:

### 1. Define Samples

```yaml
samples:
  - name: "control_rep1"
    fastq: "data/fastq/control_rep1.fastq.gz"
  - name: "control_rep2"
    fastq: "data/fastq/control_rep2.fastq.gz"
  - name: "treated_rep1"
    fastq: "data/fastq/treated_rep1.fastq.gz"
  - name: "treated_rep2"
    fastq: "data/fastq/treated_rep2.fastq.gz"
```

### 2. Specify sgRNA Library

```yaml
library:
  file: "data/library/library.txt"
```

Library file format (tab-separated):
```
sgRNA_ID    sequence    gene
sgRNA_1     ATCGATCGATCGATCG    GENE1
sgRNA_2     GCTAGCTAGCTAGCTA    GENE1
sgRNA_3     TTAATTAATTAATTAA    GENE2
```

### 3. Define Comparisons

```yaml
comparisons:
  - name: "treatment_vs_control"
    treatment: ["treated_rep1", "treated_rep2"]
    control: ["control_rep1", "control_rep2"]
```

### 4. Adjust Analysis Parameters (Optional)

```yaml
params:
  qc:
    min_quality: 20
    min_length: 18
  counting:
    mismatches: 0
  mageck:
    fdr_threshold: 0.25
    norm_method: "median"
  drugz:
    pseudocount: 5
    min_observations: 1
```

## Usage

### Basic Run

Execute the entire pipeline:

```bash
snakemake --use-conda --cores 8
```

### Dry Run

Preview what will be executed:

```bash
snakemake --use-conda --cores 8 -n
```

### Generate Workflow Diagram

```bash
snakemake --dag | dot -Tpdf > workflow.pdf
```

### Run Specific Steps

Quality control only:
```bash
snakemake --use-conda --cores 8 results/qc/{sample}_fastqc.html
```

Count generation:
```bash
snakemake --use-conda --cores 8 results/counts/merged_counts.txt
```

MAGeCK analysis:
```bash
snakemake --use-conda --cores 8 results/mageck/{comparison}.gene_summary.txt
```

## Output Files

### Quality Control
- `results/qc/{sample}_fastqc.html`: FastQC report for each sample
- `results/qc/{sample}_fastqc.zip`: FastQC data archive

### Counts
- `results/counts/{sample}_counts.txt`: Individual sample counts
- `results/counts/merged_counts.txt`: Combined count matrix

### MAGeCK Results
- `results/mageck/{comparison}.gene_summary.txt`: Gene-level statistics
- `results/mageck/{comparison}.sgrna_summary.txt`: sgRNA-level statistics

Columns in gene_summary.txt:
- `id`: Gene identifier
- `pos|score`: Positive selection score
- `pos|p-value`: P-value for positive selection
- `pos|fdr`: FDR for positive selection
- `neg|score`: Negative selection score
- `neg|p-value`: P-value for negative selection
- `neg|fdr`: FDR for negative selection

### DrugZ Results
- `results/drugz/{comparison}_drugz.txt`: DrugZ analysis results

Columns:
- `gene`: Gene identifier
- `numObs`: Number of sgRNAs per gene
- `log2FC`: Mean log2 fold change
- `normZ`: Normalized Z-score
- `pvalue`: P-value
- `fdr`: False discovery rate

### Reports
- `results/reports/analysis_summary.html`: Comprehensive analysis summary

## Analysis Methods

### MAGeCK (Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout)

MAGeCK uses a negative binomial model to identify genes where sgRNAs are systematically enriched or depleted. It's particularly robust for screens with multiple replicates.

**Key features:**
- Robust to outliers
- Accounts for variability between sgRNAs
- Provides both positive and negative selection statistics

### DrugZ

DrugZ uses a normalized Z-score approach to identify essential genes. It's particularly sensitive for detecting gene depletion in survival screens.

**Key features:**
- Simple, interpretable statistics
- Good sensitivity for depletion screens
- Works well with variable sgRNA numbers per gene

## Troubleshooting

### Issue: Pipeline fails with "No such file or directory"
**Solution**: Ensure all paths in `config.yaml` are correct and files exist.

### Issue: Conda environment creation fails
**Solution**: Update conda/mamba: `conda update -n base conda`

### Issue: Out of memory errors
**Solution**: Reduce the number of cores or increase available RAM.

### Issue: FastQC fails
**Solution**: Check FASTQ file format and compression. Ensure files are not corrupted.

## Performance Tips

1. **Use Mamba**: Faster than conda for environment creation
   ```bash
   conda install -c conda-forge mamba
   snakemake --use-conda --conda-frontend mamba --cores 8
   ```

2. **Parallel Execution**: Use more cores for faster processing
   ```bash
   snakemake --use-conda --cores 16
   ```

3. **Cluster Execution**: Submit to HPC cluster
   ```bash
   snakemake --use-conda --cluster "sbatch -c {threads} -t 4:00:00" --jobs 10
   ```

## Citation

If you use this pipeline, please cite:

- **Snakemake**: Mölder, F., et al. (2021). Sustainable data analysis with Snakemake. F1000Research, 10:33.
- **MAGeCK**: Li, W., et al. (2014). MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15:554.
- **DrugZ**: Colic, M., et al. (2019). Identifying chemogenetic interactions from CRISPR screens with drugZ. Genome Medicine, 11:52.
- **FastQC**: Andrews, S. (2010). FastQC: A quality control tool for high throughput sequence data.

## License

This pipeline is provided under the MIT License. See LICENSE file for details.

## Support

For issues, questions, or contributions, please open an issue on the GitHub repository.

## Authors

Developed for CRISPR screening analysis workflows.
