#!/bin/bash
conda activate ldsc # 采用 python2，而不是 python3
dir=/mnt/d
gwasdir=$dir/data/gwas/main/clean
refdir=$dir/data/ldsc
source $dir/scripts/f/f.sh


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# munge_sumstats 数据标准化
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dats=`cd $gwasdir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	if [[ -f $dat.sumstats.gz ]]; then echo $dat already run; continue; fi
	echo $dat
	for var in SNP SNP_col EA EA_col NEA NEA_col N N_col BETA BETA_col P P_col; do eval $var=''; done
	datf=$gwasdir/$dat.gz; header_names $datf 🏮
	if [[ -z "$N" ]]; then n_str="--N 100000"; else n_str="--N-col $N"; fi
	bad_p=`zcat $datf | awk -v p=$P_col 'NR>1 && $p<=0' | wc -l` 
	if [[ bad_p -gt 0 ]]; then 
		zcat $datf | awk -v p=$P_col '{if($p<=1e-300) $p=1e-300; print $0}' | sed 's/ /\t/g' | gzip -f > $dat.gz; datf=$dat.gz
	fi
	python /mnt/d/software/ldsc/munge_sumstats.py --chunksize 10000 --sumstats $datf --merge-alleles $refdir/w_hm3.snplist --out $dat --snp $SNP --a1 $EA --a2 $NEA $n_str --signed-sumstats $BETA,0 --p $P --ignore SNPID,OR 
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对多个GWAS计算Rg
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某些loci分析local的h2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://github.com/huilisabrina/HEELS
dir=/work/sph-huangj
module load python/anaconda3/2020.7
source activate heels
python $dir/software/heels/run_HEELS.py --output_fp jie --sumstats_fp ukb_332k_simul_pheno.txt \
	--ld_fp ukb_332k_chr22_LD.npz --ld_snp_fp ukb_332k_chr22_maf01_geno1.bim \
	--init_values 0.1,0.9 --N 332340 --tol 1e-4 --calc_var --stream-stdout

