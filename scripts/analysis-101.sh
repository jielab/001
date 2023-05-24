
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download and format UKB 35 traits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for dat in `ls *array.gz | sed 's/\.array.gz//g'`; do
	echo "#!/bin/bash
	zcat $dat.array.gz | awk '{if (\$8<1e-300) \$8=1e-300; if (NR==1) print \"SNP CHRPOS CHR POS EA NEA EAF N BETA SE Z P\"; else print \$3, \$1\":\"\$2,\$1,\$2, toupper(\$4),toupper(\$5),\"NA\",500000, \$6,\$7,\$6/\$7,\$8}' > $dat.gwas
	zcat $dat.rsid.gz | awk 'NR >1 && $15>=0.4 {if ($8<1e-300) $8=1e-300; print \$NF,\$1\":\"\$2,\$1,\$2, \$5,\$4,$14,500000, \$6,\$7,\$6/\$7,\$8}' >> $dat.gwas
	gzip -f $dat.gwas
	" > $dat.cmd
	chmod 777 $dat.cmd; ./$dat.cmd 
done
for dat in `ls ukbb_*.top.snps | sed 's/.top.snps//g'`; do
	fgrep -vw NA $dat.top.snps > tmp/$dat.top.snps
	zcat $dat.gwas.gz | fgrep -wf tmp/$dat.top.snps | awk -v f=$dat '{if (NR==1) f2="Trait"; else f2=f; print f2, $0}' > tmp/$dat.top.txt
done
awk 'NR==1 || FNR>1' tmp/*.top.txt > tmp/ukbb.top.txt


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download and format GWAS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir=/mnt/d
outdir=$dir/data/gwas/pheweb2
mrlist=$dir/files/mr.list.txt
cut -f 1,4 $mrlist | sed '1d' | fgrep -vw "NA" | while IFS=$'\t' read -r pheno_code file_url; do
	file_name=$(basename $file_url)
	file_ext=${file_name##*.}
	if [[ -f $pheno_code.gz ]]; then
		echo $pheno_code.gz already downloaded
		continue
	fi
	echo -e "#!/bin/bash
	wget $file_url
	mv $file_name $pheno_code.$file_ext
	if [[ $file_ext == zip ]]; then
		unzip -p $pheno_code.$file_ext | gzip -f > $pheno_code.gz
	fi
	if [[ $file_ext == txt || $file_ext == tsv ]]; then
		mv $pheno_code.$file_ext $pheno_code
		gzip -f $pheno_code
	fi
	" > $pheno_code.cmd
	chmod 777 $pheno_code.cmd; ./$pheno_code.cmd 
done
conda activate ldsc
ldscdir=/mnt/d/data/ldsc
for dat in `ls *gz | sed 's/\.gz$//g'`; do
	head_row=`zcat $dat.gz | head -1 | sed 's/\t/ /g'`;
	snp_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'snp|variant_id|markername'`; snp_col=`echo "$"$snp_str | sed 's/:.*//'`
	chr_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'chr|chrom|chromosome'`; chr_col=`echo "$"$chr_str | sed 's/:.*//'`
	pos_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'pos|bp|base_pair_location'`; pos_col=`echo "$"$pos_str | sed 's/:.*//'`
	ea_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'ea|alt|effect_allele|allele1'`; ea_col=`echo "$"$ea_str | sed 's/:.*//'`
	nea_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'nea|ref|noneffect_allele|other_allele|allele0'`; nea_col=`echo "$"$nea_str | sed 's/:.*//'`
	eaf_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'eaf|af|effect_allele_frequency|a1freq'`; eaf_col=`echo "$"$eaf_str | sed 's/:.*//'`
	n_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'n'`; n_col=`echo "$"$n_str | sed 's/:.*//'`
	beta_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'beta'`; beta_col=`echo "$"$beta_str | sed 's/:.*//'`
	se_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'se|standard_error'`; se_col=`echo "$"$se_str | sed 's/:.*//'`
	p_str=`echo $head_row | sed 's/\t/ /g' | tr ' ' '\n' | grep -Einw 'p|pval|p.value|p_value|p_linreg'`; p_col=`echo "$"$p_str | sed 's/:.*//'`
	echo $snp_str $chr_str $pos_str $ea_str $nea_str $eaf_str $n_str $beta_str $se_str $p_str
	echo columns $snp_col $chr_col $pos_col $ea_col $nea_col $eaf_col $n_col $beta_col $se_col $p_col 
	if [[ $chr_str == "" ]]; then chr_col="\"NA\""; fi
	if [[ $pos_str == "" ]]; then pos_col="\"NA\""; fi
	if [[ $eaf_str == "" ]]; then eaf_col=0.5; fi
	if [[ $n_str == "" ]]; then n_col=100000; fi
	echo "#!/bin/bash
	zcat $dat.gz | awk '{if ($p_col<1e-300) $p_col=1e-300; if (NR==1) print \"SNP CHRPOS CHR POS EA NEA EAF N BETA SE Z P\"; else if ($beta_col*$se_col !=0) print $snp_col,$chr_col\":\"$pos_col,$chr_col,$pos_col, toupper($ea_col), toupper($nea_col),$eaf_col,$n_col, $beta_col,$se_col,$beta_col/$se_col,$p_col}' | sort -k 3,3V -k 4,4n | gzip -f > $dat.gwas.gz
	zcat $dat.gwas.gz | awk 'NR==1 || \$NF <=1e-03' > $dat.gwas.p3
	zcat $dat.gwas.gz | awk 'NR==1 || \$NF<5e-8 {b=sprintf(\"%.0f\",\$4/1e6); print \$1,\$3,\$4,\$NF,b}' | sort -k 2,2n -k 5,5n -k 4,4g | awk '{chunk=\$2\".\"\$5; if (arr[chunk] !=\"Y\") print \$1; arr[chunk] =\"Y\"}' > $dat.top1.snps
	munge_sumstats.py --chunksize 10000 --sumstats $dat.gwas.gz --merge-alleles $ldscdir/w_hm3.snplist --out $dat --snp SNP --a1 EA --a2 NEA --ignore EAF --signed-sumstats Z,0 --p P --N-col N 
	zcat $dat.gwas.gz | fgrep -wf /mnt/d/data/hm3/hm3.b37.snps | awk '{print \$1,\$5,\$6,\$8,\$11}' | sed '1 s/EA NEA/A1 A2/' > $dat.hdl.txt
	" > $dat.cmd
	chmod 777 $dat.cmd; ./$dat.cmd 
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LDSC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conda activate ldsc
ldscdir=/mnt/d/data/ldsc
ldrefdir=$ldscdir/eur_w_ld_chr # $ldscdir/weights_hm3_no_hla
trait=covid_resf; traits="covid_resf covid_icu covid_r5 covid_r6"
traits_c=`echo $traits | sed 's/ /,/g'`; traits_cp=`echo $traits_c |  sed 's/^/P./; s/,/,P./g'`
ldsc.py --h2 $trait.sumstats.gz --ref-ld-chr $ldrefdir/ --w-ld-chr $ldrefdir/
#ldsc.py --h2 $trait.sumstats.gz --ref-ld-chr $ldrefdir/CNS.,baseline. --w-ld-chr $ldrefdir/weights. --overlap-annot --frqfile-chr $datadir/1000G.mac5eur. --out $trait --print-coefficients
echo $traits | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % ldsc.py --rg % --out $trait.rg --ref-ld-chr $ldrefdir/ --w-ld-chr $ldrefdir/
awk '$1=='Summary' {printf NR}' $trait.rg.log | xargs -n1 -I % awk -v s=% 'FNR >=s' $trait.rg.log | sed 's/.sumstats.gz//g' > $trait.rg.txt
