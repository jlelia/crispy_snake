# sgRNA Library Files

Place your sgRNA library file(s) in this directory.

## File Format

Tab-separated file with the following columns:

```
sgRNA_ID    sequence    gene
```

### Example:

```
sgRNA_BRCA1_1    ATCGATCGATCGATCG    BRCA1
sgRNA_BRCA1_2    GCTAGCTAGCTAGCTA    BRCA1
sgRNA_BRCA2_1    TTAATTAATTAATTAA    BRCA2
sgRNA_CTRL_1     AAAAAAAAAAAAAAAA    NonTargeting
```

## Requirements

- **sgRNA_ID**: Unique identifier for each sgRNA
- **sequence**: 20-23 bp sgRNA sequence (without PAM)
- **gene**: Target gene name (use "NonTargeting" for negative controls)

## Common Libraries

Popular genome-wide sgRNA libraries include:
- Brunello (Broad Institute)
- GeCKO v2 (Feng Zhang lab)
- Toronto KnockOut Library (TKO)
- Yusa v1/v1.1

Ensure your library matches your experimental design and organism.
