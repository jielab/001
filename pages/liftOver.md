# liftOver notes: GWAS/COJO table and VCF conversion between GRCh38 and GRCh37

This note records a practical workflow for using `bcftools +liftover` to convert variant coordinates between genome builds.

Main use cases:

```text
1. LiftOver a GWAS/COJO summary table with CHR/POS/REF/ALT columns from GRCh38 to GRCh37.
2. LiftOver whole VCF/VCF.GZ files between GRCh37 and GRCh38.
3. Avoid using rsID as the primary join key; use a temporary row ID instead.
```

Important note:

```text
bcftools +liftover is VCF-aware. It is more suitable than UCSC liftOver for variants because it can handle REF/ALT, reverse-complement alleles, indels, and multi-allelic variants.
For GWAS/COJO tables, first convert the table into a temporary minimal VCF, run bcftools +liftover, and then merge the lifted coordinates back to the original table by row ID.
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s0: Install bcftools +liftover plugin
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

```bash
mkdir -p /mnt/d/software/bcftools_plugins
cd /mnt/d/software/bcftools_plugins

wget -c https://software.broadinstitute.org/software/score/score_1.22-20250819.zip
unzip -o score_1.22-20250819.zip

echo 'export BCFTOOLS_PLUGINS=/mnt/d/software/bcftools_plugins' >> ~/.bashrc
source ~/.bashrc

bcftools +liftover -h
```

If `bcftools +liftover -h` does not work, check whether the plugin `.so` files are directly under `$BCFTOOLS_PLUGINS`.

```bash
echo $BCFTOOLS_PLUGINS
find /mnt/d/software/bcftools_plugins -name "liftover.so"
bcftools +liftover -vv
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s1: Download chain files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Save chain files to `/mnt/d/files/liftOver`.

```bash
mkdir -p /mnt/d/files/liftOver
cd /mnt/d/files/liftOver

wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
wget -c https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s2: Download reference FASTA files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Save FASTA files to `/mnt/d/data/refGen/fasta`.

```bash
mkdir -p /mnt/d/data/refGen/fasta
cd /mnt/d/data/refGen/fasta

wget -O human_g1k_v37.fasta.gz https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
wget -O GCA_GRCh38.fna.gz https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz
wget -O chm13v2.0.fa.gz https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/chm13v2.0.fa.gz
```

Check whether the downloaded files are valid gzip files.

```bash
gzip -t human_g1k_v37.fasta.gz GCA_GRCh38.fna.gz chm13v2.0.fa.gz
```

Decompress and create FASTA indexes.

```bash
gunzip -c human_g1k_v37.fasta.gz > human_g1k_v37.fasta
gunzip -c GCA_GRCh38.fna.gz > GCA_GRCh38.fna
gunzip -c chm13v2.0.fa.gz > chm13v2.0.fa

samtools faidx human_g1k_v37.fasta
samtools faidx GCA_GRCh38.fna
samtools faidx chm13v2.0.fa
```

Recommended path variables:

```bash
fa37=/mnt/d/data/refGen/fasta/human_g1k_v37.fasta
fa38=/mnt/d/data/refGen/fasta/GCA_GRCh38.fna

chain38to37=/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz
chain37to38=/mnt/d/files/liftOver/hg19ToHg38.over.chain.gz
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s3: Create script for GWAS/COJO table liftover, GRCh38 → GRCh37
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This script expects a tab-delimited GWAS/COJO table with chromosome, position, reference allele, and alternate allele columns.

Example column names:

```text
CHR POS REF ALT
```

or, if the GWAS file uses effect/non-effect alleles and they are confirmed to represent REF/ALT:

```text
CHR POS NEA EA
```

Important warning:

```text
EA/NEA are not always equivalent to ALT/REF. bcftools +liftover requires true REF and ALT alleles. If REF does not match the source reference FASTA, many records will be rejected or incorrectly represented.
```

Create `lift38to37.sh`.

```bash
cat > lift38to37.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

dat=$1
chr=${2:-CHR}
pos=${3:-POS}
ref=${4:-REF}
alt=${5:-ALT}

fa38=${fa38:-/mnt/d/data/refGen/fasta/GCA_GRCh38.fna}
fa37=${fa37:-/mnt/d/data/refGen/fasta/human_g1k_v37.fasta}
chain=${chain:-/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz}

awk -v chr=$chr -v pos=$pos -v ref=$ref -v alt=$alt '
BEGIN{
	FS=OFS="\t"
	print "##fileformat=VCFv4.2"
	print "#CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO"
}
NR==1{
	for(i=1;i<=NF;i++) h[$i]=i
	next
}
{
	c=$h[chr]
	gsub(/^chr/,"",c)
	if(c==23)c="X"
	if(c==24)c="Y"
	print "chr"c,$h[pos],"RID"NR,$h[ref],$h[alt],".",".","."
}' "$dat" | bcftools view -Oz -o "$dat.tmp38.vcf.gz"

tabix -f -p vcf "$dat.tmp38.vcf.gz"

bcftools +liftover --no-version -Ou "$dat.tmp38.vcf.gz" -- \
	-s "$fa38" \
	-f "$fa37" \
	-c "$chain" \
	--reject "$dat.reject37.bcf" \
	--reject-type b \
	| bcftools sort -Oz -o "$dat.tmp37.vcf.gz" --write-index

bcftools query -f '%ID\t%CHROM\t%POS\t%REF\t%ALT\n' "$dat.tmp37.vcf.gz" \
	| sed 's/^RID//; s/\tchr/\t/' > "$dat.pos37.tsv"

Rscript -e 'library(data.table); dat=commandArgs(T)[1]; x=fread(dat); x[,RID:=.I+1L]; y=fread(paste0(dat,".pos37.tsv"),col.names=c("RID","CHR.37","POS.37","REF.37","ALT.37")); z=merge(x,y,by="RID",all.x=TRUE,sort=FALSE); z[,RID:=NULL]; fwrite(z,paste0(dat,".lift37.tsv.gz"),sep="\t")' "$dat"
EOF

chmod +x lift38to37.sh
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s4: Run GWAS/COJO table liftover
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If the table has true `CHR POS REF ALT` columns:

```bash
./lift38to37.sh XYZ.txt CHR POS REF ALT
```

If the table uses `NEA` as REF and `EA` as ALT, and you have confirmed this is correct:

```bash
./lift38to37.sh XYZ.txt CHR POS NEA EA
```

Output:

```text
XYZ.txt.lift37.tsv.gz
```

The output table keeps the original columns and adds:

```text
CHR.37 POS.37 REF.37 ALT.37
```

Rejected variants are saved as:

```text
XYZ.txt.reject37.bcf
```

Check reject count:

```bash
bcftools view -H XYZ.txt.reject37.bcf | wc -l
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s5: Batch liftover for multiple GWAS/COJO files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For multiple plain text files:

```bash
ls -1 *.txt | parallel -j 4 './lift38to37.sh {} CHR POS REF ALT'
```

For multiple gzipped files:

```bash
ls -1 *.gz | parallel -j 4 'zcat {} > {.}.txt && ./lift38to37.sh {.}.txt CHR POS REF ALT'
```

If the columns are `CHR POS NEA EA`:

```bash
ls -1 *.gz | parallel -j 4 'zcat {} > {.}.txt && ./lift38to37.sh {.}.txt CHR POS NEA EA'
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s6: LiftOver whole VCF.GZ files, GRCh37 → GRCh38
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create output folders.

```bash
mkdir -p GRCh38 reject tmp
```

Set path variables.

```bash
fa37=/mnt/d/data/refGen/fasta/human_g1k_v37.fasta
fa38=/mnt/d/data/refGen/fasta/GCA_GRCh38.fna
chain=/mnt/d/files/liftOver/hg19ToHg38.over.chain.gz

export fa37 fa38 chain
```

Run liftover for all VCF files in the current folder.

```bash
ls -1 *.vcf.gz | parallel -j 2 '
	f={}
	b=$(basename "$f" .vcf.gz)
	echo "$f"

	bcftools +liftover --no-version -Ou "$f" -- \
		-s "$fa37" \
		-f "$fa38" \
		-c "$chain" \
		--reject "reject/${b}.reject38.bcf" \
		--reject-type b \
		| bcftools sort -m 4G -T "tmp/$b" -Oz -o "GRCh38/${b}.GRCh38.vcf.gz" --write-index
'
```

Notes:

```text
-j 2 means two VCF files are processed at the same time.
bcftools sort -m 4G means each sort process can use around 4 GB memory.
Do not set -j too high for very large VCF files because sorting is I/O- and memory-heavy.
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s7: LiftOver whole VCF.GZ files, GRCh38 → GRCh37
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create output folders.

```bash
mkdir -p GRCh37 reject tmp
```

Set path variables.

```bash
fa38=/mnt/d/data/refGen/fasta/GCA_GRCh38.fna
fa37=/mnt/d/data/refGen/fasta/human_g1k_v37.fasta
chain=/mnt/d/files/liftOver/hg38ToHg19.over.chain.gz

export fa38 fa37 chain
```

Run liftover for all VCF files in the current folder.

```bash
ls -1 *.vcf.gz | parallel -j 2 '
	f={}
	b=$(basename "$f" .vcf.gz)
	echo "$f"

	bcftools +liftover --no-version -Ou "$f" -- \
		-s "$fa38" \
		-f "$fa37" \
		-c "$chain" \
		--reject "reject/${b}.reject37.bcf" \
		--reject-type b \
		| bcftools sort -m 4G -T "tmp/$b" -Oz -o "GRCh37/${b}.GRCh37.vcf.gz" --write-index
'
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s8: Re-bgzip and tabix VCF files if needed
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If tabix reports:

```text
not compressed with bgzip
```

then the `.vcf.gz` file was compressed by ordinary gzip, not bgzip. Recompress it with bgzip.

Single file:

```bash
f=chr13.noRB.vcf.gz

mv "$f" "$f.old.gz"
gzip -dc "$f.old.gz" | bgzip -@ 8 -c > "$f"
tabix -f -p vcf "$f"
```

Batch version:

```bash
ls -1 *.vcf.gz | parallel -j 4 'echo {}; mv {} {}.old.gz; gzip -dc {}.old.gz | bgzip -@ 4 -c > {} && tabix -f -p vcf {}'
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s9: Useful checks
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check FASTA chromosome names.

```bash
grep "^>" /mnt/d/data/refGen/fasta/human_g1k_v37.fasta | head
grep "^>" /mnt/d/data/refGen/fasta/GCA_GRCh38.fna | head
```

Check VCF chromosome names.

```bash
bcftools view -h input.vcf.gz | grep "^##contig" | head
bcftools view -H input.vcf.gz | head
```

Check whether a VCF is bgzip-compressed and indexed.

```bash
file input.vcf.gz
tabix -l input.vcf.gz | head
```

Check liftover reject variants.

```bash
bcftools view -H input.reject37.bcf | wc -l
bcftools view -H input.reject38.bcf | wc -l
```

Check the first few lifted variants.

```bash
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' output.vcf.gz | head
```


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 s10: Practical recommendation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For GWAS/COJO plus 1000G/archaic VCF analysis:

```text
If possible, keep 1000G and archaic VCF files in their original GRCh37/hg19 coordinates.
Convert only GWAS/COJO lead SNP tables or summary statistics to GRCh37.
This is usually safer and faster than lifting over the entire 1000G or archaic VCF dataset.
```

For variant IDs:

```text
Do not use rsID as the only merge key.
Use CHR:POS:REF:ALT plus genome build as the main key.
Use rsID only as an annotation.
For long-term cross-reference representation, consider SPDI or GA4GH VRS-style identifiers.
```
