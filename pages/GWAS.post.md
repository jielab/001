# GWAS summary statistics QC checklist

This page records a practical checklist for cleaning and preparing GWAS summary statistics before downstream analyses such as PheWeb, FUMA, LDSC, MR, colocalization, or locus-level visualization.

## 1. Recommended column names

A clean GWAS file is easier to reuse if the core columns are standardized.

```text
SNP  CHR  POS  CHRPOS  EA  NEA  EAF  N  BETA  SE  Z  P
```

Minimal recommended columns:

```text
SNP  CHR  POS  EA  NEA  EAF  N  BETA  SE  P
```

where:

- `SNP`: rsID or variant ID
- `CHR`, `POS`: chromosome and base-pair position
- `EA`, `NEA`: effect allele and non-effect allele
- `EAF`: effect allele frequency
- `N`: sample size
- `BETA`, `SE`, `P`: effect estimate, standard error, and P value

## 2. Check whether every row has the same number of columns

This is the “first button” of GWAS QC. If rows have different numbers of fields, many tools will silently fail or misread columns.

```bash
zcat GWAS.gz | awk '{print NF}' | sort -nu
```

If the output has more than one value, inspect the file before continuing. For tab-delimited files, a quick repair for empty fields is:

```bash
sed -e 's/^\t/NA\t/' -e 's/\t\t/\tNA\t/g' -e 's/\t\t/\tNA\t/g' -e 's/\t$/\tNA/' input.tsv > input.fixed.tsv
```

## 3. Standardize column names

Example command for remapping arbitrary column names to a standard format. Replace the shell variables with the column numbers in your file.

```bash
header_names $file

cat $file | awk \
	-v snp_col=$SNP_col \
	-v chr_col=$CHR_col \
	-v pos_col=$POS_col \
	-v ea_col=$EA_col \
	-v nea_col=$NEA_col \
	-v eaf_col=$EAF_col \
	-v n_col=$N_col \
	-v beta_col=$BETA_col \
	-v se_col=$SE_col \
	-v p_col=$P_col \
	-v log10p_col=$LOG10P_col \
'{
	if(NR==1) print "SNP CHR POS EA NEA EAF N BETA SE P";
	else{
		if(p_col=="") pval=10^(-$log10p_col); else pval=$p_col;
		print $snp_col,$chr_col,$pos_col,$ea_col,$nea_col,$eaf_col,$n_col,$beta_col,$se_col,pval
	}
}' | grep -vwE "nan|inf|NA" | sed 's/  */\t/g' | sort -k2,2n -k3,3n > chr$chr.$phe
```

## 4. Basic GWAS corrections

Recommended checks:

1. Sort by chromosome and position.
2. Convert allele columns to upper case.
3. Avoid extremely small P values being parsed as zero by downstream tools.
4. Reconstruct missing `BETA`, `SE`, or `P` if only one of them is missing.

```bash
sort -k1,1V -k2,2n input.tsv > input.sorted.tsv
```

For very small P values, some tools may fail when `P < 1e-312`. A practical workaround is to use `Z` if available, or cap extremely small P values at a safe floor such as `1e-300`.

Useful formulas:

```r
# If beta and p are available
se <- abs(beta / qnorm(p / 2))

# If se and p are available
beta <- se * qnorm(p / 2)

# If beta and se are available
p <- 2 * pnorm(-abs(beta / se))

# If 95% CI is available
se <- (CI_upper - CI_lower) / (1.96 * 2)
```

## 5. Add rsID, CHR, POS, REF, or ALT when needed

If a GWAS file is missing rsID or genomic coordinates, a reference `.pvar`/variant map can be used to supplement the missing fields.

Examples:

```bash
# Input has SNP, add CHR/POS
python /mnt/d/scripts/f/0add_rsid.py \
	-i afib.txt \
	--sep $'\t' \
	--snp SNP \
	--chr chr_name \
	--pos hm_pos \
	--ref effect_allele \
	--alt other_allele \
	-d /mnt/d/data/ukb/gen/imp/ukb.imp.pvar.gz \
	-o out.tsv.gz
```

Notes:

- If the input only has `SNP`, add `--chr` and `--pos`.
- If the input only has `CHR` and `POS`, add `--snp`.
- `--ref` and `--alt` are optional but useful when allele-aware matching is needed.

## 6. Merge lifted coordinates or external annotations

Example:

```bash
python scripts/library/join_file.py \
	-i "$dat,TAB,0 $dat.lifted.3col,TAB,2" \
	-o $dat.NEW.tmp

sed -i 's/  */\t/g' $dat.NEW.tmp
awk '$NF=="NA"' $dat.NEW.tmp | wc -l
cut -f 1-10,12 $dat.NEW.tmp | sed '1 s/POS/POS.b38/' > $dat.NEW.txt
```

## 7. Slim large GWAS files

Some web servers cannot accept very large summary-statistic files. For example, FUMA may reject files larger than its upload limit. One practical approach is to keep common variants and round selected numeric columns.

```bash
zcat $dat.gz | awk '
function r(x){return sprintf("%.4f", x)}
{
	if(NR==1) print;
	else if($6 > 0.005 && $6 <= 0.995){
		$6=r($6); $8=r($8); $9=r($9); print
	}
}' | sed 's/  */\t/g' | bgzip > $dat.lean.gz
```

## 8. Index GWAS files

For tabix-indexed summary statistics:

```bash
tabix -f -S 1 -s 1 -b 2 -e 2 GWAS.gz
```

Typical assumptions:

```text
-S 1: skip one header line
-s 1: chromosome column
-b 2: start position column
-e 2: end position column
```
