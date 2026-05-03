# Mendelian randomization practical notes

This page records small but important details that often determine whether a Mendelian randomization analysis is interpretable.

## 1. Instrument selection

Common filters:

- Remove weak instruments, often defined as `F_stat < 10`.
- Remove variants in complex LD regions when appropriate, especially the MHC/HLA region.
- LD-clump instruments before harmonization unless the method explicitly models correlated instruments.
- Check whether the exposure and outcome samples overlap.

MHC/HLA region examples:

```text
GRCh37: chr6:28477897-33448354
GRCh38: chr6:28510120-33480577
```

The Michigan Imputation Server uses a broader HLA-region definition:

```text
rs57232568 to rs116206961
GRCh37: 29,000,554-33,999,992
GRCh38: 29,032,777-34,032,215
```

## 2. Allele harmonization

The `allele.qc`-style workflow is useful for aligning `BETA` and `EAF` between exposure and outcome datasets. One practical detail: after harmonization, the output table may still keep the original `EA` and `NEA` column labels, so always check the direction of `BETA` explicitly.

Recommended checks:

```text
SNP match
EA/NEA match or strand flip
BETA direction
EAF consistency
Palindromic SNP handling
Missing N or EAF
```

## 3. Common file and software pitfalls

1. `fread()` and `fwrite()` are much faster than `read.table()` and `write.table()`, but check quote and missing-value behavior.
2. If `EAF` is exactly `0`, some workflows may output an empty field; check that the output still has the expected number of columns.
3. `gcta --cojo-slct` may generate `.ldr.cojo` files with an extra trailing tab.
4. Some input GWAS files may not have `EAF` or `N`; code should handle one or both files missing these columns.
5. Long chains of `%>%` can make debugging harder. Use intermediate objects when harmonization or filtering becomes complex.
6. After `group_by()`, remember to `ungroup()` before later operations.
7. In PLINK `.bim` files, chromosome X can be written as `X`, but GCTA usually expects `23`.
8. Always confirm whether coordinates are GRCh37 or GRCh38 before applying MHC/HLA filters.

## 4. Robustness checks

MR results can be sensitive to instrument definition, LD clumping parameters, palindromic SNP handling, and pleiotropic outliers. A result should not be treated as robust until it has passed several sensitivity checks.

Useful checks include:

```text
IVW estimate
Weighted median estimate
MR-Egger intercept
Leave-one-out analysis
Heterogeneity test
Outlier or pleiotropy checks
Colocalization, if cis instruments or locus-level evidence are used
```

A useful cautionary example is the discussion around `cisMRcML` robustness: <https://github.com/ZhaotongL/cisMRcML/issues/6>
