#!/bin/bash

conda activate ldsc # 采用 python2，而不是 python3
gwasdir=/mnt/d/data/gwas/main/ldsc
refdir=/mnt/d/data/ldsc
Arr1=("SNP" "CHR" "POS" "EA" "NEA" "EAF" "N" "BETA" "SE" "P")
Arr2=("snp|rsid|variant_id" "chr|chrom|chromosome" "pos|bp|base_pair" "alt|ea|eff.allele|effect_allele|a1|allele1" "ref|nea|ref.allele|other_allele|a2|allele0" "eaf|a1freq|effect_allele_freq" "n|Neff" "beta" "se|standard_error" "p|pval|p_bolt_lmm")

dats=`cd $gwasdir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	if [ -f $dat.sumstats.gz ]; then echo $dat already run; continue; fi
	echo $dat
	datf=$gwasdir/$dat.gz
	head_row=`zcat $datf | head -1 | sed 's/\t/ /g'`
	snp=""; ea=""; nea=""; n=""; beta=""; p=""
	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/\w*://'`
		eval ${Arr1[$i]}_col=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:\w*//'`
	done
	echo dat $dat, snp $SNP $SNP_col, ea $EA $EA_col, nea $NEA $NEA_col, n $N $N_col, beta $BETA $BETA_col, p $P $P_col
	if [[ -z "$N" ]]; then n_str="--N 100000"; else n_str="--N-col $N"; fi
	bad_p=`zcat $datf | awk -v p=$P_col 'NR>1 && $p<=0' | wc -l` # awk和muange_sumstats.py都会把 <1e-312当作0。 直接用Z，不跑muange_sumstats.py，就能避开这个问题。
	if [[ bad_p -gt 0 ]]; then 
		zcat $datf | awk -v p=$P_col '{if($p<=1e-300) $p=1e-300; print $0}' | sed 's/ /\t/g' | gzip -f > $dat.gz
		datf=$dat.gz
	fi
	python /mnt/d/software/ldsc/munge_sumstats.py --chunksize 10000 --sumstats $datf --merge-alleles $refdir/w_hm3.snplist --out $dat --snp $SNP --a1 $EA --a2 $NEA $n_str --signed-sumstats $BETA,0 --p $P --ignore SNPID,OR 
done
for dat in $dats; do
	echo $dat $dats | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % \
	/mnt/d/software/ldsc/ldsc.py --rg % --out $dat.rg --ref-ld-chr $refdir/ --w-ld-chr $refdir/
	beginl=`awk '$1=="Summary" {printf NR}' $dat.rg.log` 
	awk -v s=$beginl 'FNR >s' $dat.rg.log | head -n -3 | sed 's/.sumstats.gz//g'  > $dat.rg.txt
done
awk 'NR==1 || FNR>1' *.rg.txt | sed 's/  */ /g; s/^ //' > all.rg.res


## 下面是R代码
library(corrplot); library(reshape); library(corrplot)
dat0 <- read.table("D:/analysis/ldsc/all.rg.res", header=T, as.is=T)
dat <- reshape::cast(dat0, p1 ~ p2, value="rg") 
dat <- dat[, -1]; rownames(dat) <- dat$p1; dat$p1 <- NULL
write.table(dat, "dat.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, append=FALSE)
corrplot(dat, method="number", is.corr = FALSE) # , tl.col = "black", tl.srt = 45