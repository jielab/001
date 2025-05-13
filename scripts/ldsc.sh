#!/bin/bash

dir0=/mnt/d
gwasdir=$dir0/data/gwas/main/clean
refdir=$dir0/data/ldref/ldsc
source $dir0/scripts/f/phe.f.sh

conda activate ldsc # 🏮 采用 python2，而不是 python3

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# munge_sumstats 数据标准化
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dats=`cd $gwasdir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	if [[ -f $dat.sumstats.gz ]]; then echo $dat already run; continue; fi
	echo $dat
	for var in SNP SNP_col EA EA_col NEA NEA_col N N_col BETA BETA_col P P_col; do eval $var=''; done
	datf=$gwasdir/$dat.gz; header_names $datf # 🏮
	if [[ -z "$N" ]]; then n_str="--N 100000"; else n_str="--N-col $N"; fi
	bad_p=`zcat $datf | awk -v p=$P_col 'NR>1 && $p<=0' | wc -l` 
	if [[ bad_p -gt 0 ]]; then 
		zcat $datf | awk -v p=$P_col '{if($p<=1e-300) $p=1e-300; print $0}' | sed 's/ /\t/g' | gzip -f > $dat.gz; datf=$dat.gz
	fi
	python2 $dir0/software/ldsc/munge_sumstats.py --chunksize 10000 --sumstats $datf --merge-alleles $refdir/w_hm3.snplist --out $dat --snp $SNP --a1 $EA --a2 $NEA $n_str --signed-sumstats $BETA,0 --p $P --ignore SNPID,OR 
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对多个GWAS计算Rg
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dats=($dats) # 🏮
for ((i=0; i<${#dats[@]}-1; i++)); do
    first=${dats[i]}
    other=("${dats[@]:i+1}") # 这样可避免一个pair计算两次，导致后续R代码出现问题
    echo $first ${other[*]} | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -I % \
	python2 /mnt/d/software/ldsc/ldsc.py --rg % --out $first.rg --ref-ld-chr $refdir/ --w-ld-chr $refdir/
	beginl=`awk '$1=="Summary" {printf NR}' $first.rg.log` 
	awk -v s=$beginl 'FNR >s' $first.rg.log | head -n -3 | sed 's/.sumstats.gz//g'  > $first.rg.txt
done
awk 'NR==1 || FNR>1' *.rg.txt | sed 's/  */ /g; s/^ //' > all.rg.res

## 下面是R代码
library(tidyverse); library(reshape); library(corrplot)
dat0 <- read.table("D:/analysis/ldsc/all.rg.res", header=T, as.is=T) %>% 
	filter(grepl("bald|asd|vte|cad|^st", p1), grepl("bald|asd|vte|cad|^st", p2)) %>%
	dplyr::select(p2, p1, rg) %>% mutate(rg=ifelse(rg > 1, 1, ifelse(rg < -1, -1, rg)))	%>% drop_na() 
	# dat.tmp1 <- dat0[, c("p1", "p2", "rg")] %>% dplyr::rename(p1.tmp=p1, p2.tmp=p2)
	# dat.tmp2 <- dat0[, c("p2", "p1", "rg")] %>% dplyr::rename(p1.tmp=p2, p2.tmp=p1)
	# dat <- rbind(dat.tmp1, dat.tmp2) %>% dplyr::rename(p1=p1.tmp, p2=p2.tmp) %>% drop_na() %>% distinct(p1, p2, .keep_all = TRUE) 
dat <- dat0 %>% pivot_wider(names_from=p1, values_from=rg) %>% as.data.frame()
	write.table(dat, "dat.txt", row.names=TRUE, col.names=TRUE, quote=FALSE, append=FALSE)
	rownames(dat) <- dat$p2; dat$p2 <- NULL; dat[1:20, 1:20] 
	common <- intersect(rownames(dat), colnames(dat)); dat <- dat[common, common]
	mat <- as.matrix(dat); mat[upper.tri(mat)] <- t(mat)[upper.tri(mat)]; diag(mat) <- 1; mat[1:10, 1:10]
	col <- colorRampPalette(c("blue", "white", "red"))(200)
	corrplot(mat, type = "lower", is.corr = FALSE, diag = FALSE, col = col, tl.pos = "lt", tl.cex = 0.8, tl.col = "black", tl.srt = 90)
	corrplot(mat, type = "upper", method = "number", number.cex = 0.6, format = "%.2f", tl.pos = "n", add = TRUE, diag = FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某些loci分析local的h2
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LD block: https://github.com/jmacdon/LDblocks_GRCh38
# FTO chr16:53703962-54114467，属于LD block 1114 [chr16:51034025-53815773]。 
conda create -n heels python=3.12.2 # https://github.com/huilisabrina/HEELS
	conda activate heels
	pip install numpy pandas scipy pandas-plink joblib scikit-learn
	conda install -c conda-forge scikit-learn-intelex 
# module load python/anaconda3/2020.7 # HPC
# source activate heels
	
dat=bmi; mkdir -p $dat
chr=22; pos1=30164130; pos2=37164130
datf=$gwasdir/$dat.gz; header_names $datf #🏮

zcat $datf | awk -v chr=$chr -v pos1=$pos1 -v pos2=$pos2 -v snp_col=$SNP_col -v chr_col=$CHR_col -v pos_col=$POS_col -v ea_col=$EA_col -v nea_col=$NEA_col -v n_col=$N_col -v beta_col=$BETA_col -v se_col=$SE_col -v p_col=$P_col 'NR==1 || ($chr_col==chr && $pos_col >=pos1 && $pos_col <=pos2 && $se_col !=0 && $2 ~ /^[0-9.]+$/) {if (NR==1) print "CHR BP Predictor A1 A2 n Z P"; else print $chr_col,$pos_col,$snp_col, $ea_col,$nea_col,$n_col, $beta_col/$se_col,$p_col}' | sed 's/ /\t/g' | sort -k 1,1n -k 2,2n > $dat.txt
awk 'NR>1 {print $1, $3, 0, $2, $4, $5}' $dat.txt > $dat.bim
python $dir0/software/heels/run_HEELS.py --output_fp $dat --sumstats_fp $dat.txt \
	--ld_fp $dir0/data/ldref/heels/ukb_332k_chr22_LD.npz --ld_snp_fp $dat \
	--init_values 0.1,0.9 --N 300000 --tol 1e-4 --calc_var --stream-stdout

