# PRS / ProtRS notes

A polygenic risk score is an individual-level score calculated from genotype dosage and variant weights.

```text
PRS_i = sum_j dosage_ij × weight_j
```

Two files are conceptually different:

```text
Reference scoring file: SNP, effect allele, weight
Individual score file:  participant ID, calculated PRS
```

## Practical checklist

1. Choose a scoring file with clear ancestry, trait definition, genome build, and allele coding.
2. Harmonize variants between the scoring file and genotype data.
3. Remove ambiguous variants when allele frequency cannot resolve strand issues.
4. Calculate PRS using PLINK, PRSice, LDpred, PRS-CS, or another method.
5. Standardize the score in the target cohort.
6. Evaluate performance using discrimination, calibration, and risk reclassification.
7. Compare PRS with clinical risk factors rather than reporting PRS alone.

## Public resource

- [PRS Catalog](https://www.pgscatalog.org/)
