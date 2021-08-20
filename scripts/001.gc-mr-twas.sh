#!/bin/bash

dir=/mnt/d/data/gwas/pheweb
ldsc_dir=/mnt/d/software_lin/ldsc
traits=`cd $dir; ls *.gcta.txt | sed 's/.gcta.txt//' | tr '\n' ' '`
#traits_x="bilirubin bmi.female bmi bmi.male calcium ck creatinine crp dbp hba1c hdl ht ldl menarche menopause monocyte na plt pulsep rbc sbp tg ua wbc"
#traits_y="arrhythmia asthma breastcancer cad cataract chf copd lungcancer osteoporosis pad ra stroke t2d"


## 第一步：LDSC (https://github.com/bulik/ldsc) ##
conda activate ldsc
for trait in $traits_x $traits_y; do
	# zcat $trait.gz | sed 's/\t\t/\tNA\t/g' | gzip -f > $trait.$dat.chk # 确保没有确实数据
	python2 $ldsc_dir/munge_sumstats.py --sumstats $dir/bbj.$trait.gz --chunksize 10000 --snp rsids --a1 alt --a2 ref --frq maf --p pval --N 200000 --signed-sumstats beta,0 --merge-alleles $ldsc_dir/hm3.snplist --out $trait
	# python2 $ldsc_dir/munge_sumstats.py --sumstats $dir/$trait.gwas.gz --chunksize 10000 --snp SNP --a1 A1 --a2 A2 --frq A1_FREQ --p P --N-col N --signed-sumstats Z,0 --merge-alleles $ldsc_dir/hm3.snplist --out $trait
done
for trait in $traits_x $traits_y; do
	echo process $trait
    echo $trait $traits_x $traits_y | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % /mnt/d/software_lin/ldsc/ldsc.py --rg % --out $trait.rg --ref-ld-chr $ldsc_dir/eur_w_ld_chr/ --w-ld-chr $ldsc_dir/eur_w_ld_chr/
    awk '$1=="Summary" {printf NR}' $trait.rg.log | xargs -n1 -I % awk -v s=% 'FNR >=s' $trait.rg.log | sed 's/.sumstats.gz//g' > $trait.rg.txt
done


## 第二步：GSMR (https://cnsgenomics.com/software/gcta/#GSMR) ** 请用 --gsmr2-beta 而不是 --gsmr-beta ##
for trait in $traits_x $traits_y; do
#	zcat $dir/$trait.gwas.gz | awk '{if (NR==1) print "SNP A1 A2 freq b se p N"; else if (NF==18 && arr[$5] !="Y") print $3,$6,$7,$11, $15,$16,$18,$14; else; arr[$5]="Y"}' > $trait.gcta.txt	
	fi
done
for tx in $traits; do
	echo "$tx $dir/$tx.gcta.txt" > $tx.exposure
	for ty in $traits; do
		if [[ $ty == $tx ]]; then
			echo $tx $ty the same
			continue
		fi
		if [[ ! -f $ty.outcome ]]; then
			echo "$ty $dir/$ty.gcta.txt" > $ty.outcome
		fi
		if [[ ! -f $dir/$tx.top.snps ]]; then
			continue
		fi
		if [[ ! -f $tx.top.txt ]]; then
			cat $dir/$tx.gcta.txt | fgrep -wf $dir/$tx.top.snps > $tx.top.txt
		fi
		if [[ ! -f $ty.for.$tx.top.txt ]]; then
			cat $dir/$ty.gcta.txt | fgrep -wf $dir/$tx.top.snps > $ty.for.$tx.top.txt 
		fi
		echo -e "
		source('/mnt/d/scripts/library/gsmr_plot.r')
		gsmr_data = read_gsmr_data('$tx.$ty.eff_plot.gz')
		pdf('$tx.$ty.gcta.pdf')
		par(mai=c(1,1,0.5,0.5))
		plot_gsmr_effect(gsmr_data, '$tx', '$ty', colors()[75])
		dev.off()
		" > $tx.$ty.gcta-plot.R
		echo -e "
		source('/mnt/d/scripts/library/compareB.f.R')
		pdf('$tx.$ty.compareB.pdf')
		par(mfrow=c(3,2), mai=c(0.5,1,0.5,0.5))
		compareB(
			f1='$tx.top.txt', f1_name='$tx', f1_snp='SNP', f1_ea='A1', f1_nea='A2', f1_eaf='freq', f1_beta='b', f1_se='se', f1_p='p',
			f2='$ty.for.$tx.top.txt', f2_name='$ty', f2_snp='SNP', f2_ea='A1', f2_nea='A2', f2_eaf='freq', f2_beta='b', f2_se='se', f2_p='p'
		)
		dev.off()
		" > $tx.$ty.compareB.R
		echo -e "
		library('MendelianRandomization')
		png('$tx.$ty.mr.png', w=800, res=128)
		dat = read.table('$tx.$ty.merged.txt', header=T)
		XGb = dat\$BETA.aligned_1; XGse = dat\$SE_1; YGb = dat\$BETA.aligned_2; YGse = dat\$SE_2
		mr_forest(mr_input(XGb, XGse, YGb, YGse), snp_estimates=F, methods=c('ivw', 'median', 'wmedian', 'egger', 'maxlik', 'mbe', 'conmix'))
		dev.off()
		" > $tx.$ty.mr.R
		echo run gsmr on $tx vs. $ty
		Rscript $tx.$ty.compareB.R
		Rscript $tx.$ty.mr.R
		gcta64 --bfile /mnt/d/data/hm3/hm3.b37 --gsmr-file $tx.exposure $ty.outcome --gsmr2-beta --gsmr-direction 2 --diff-freq 1 --gwas-thresh 5e-8 --effect-plot --out $tx.$ty
		Rscript $tx.$ty.gcta-plot.R
		awk -v t=$tx 'NR==1 || $1==t' $tx.$ty.gsmr > $tx.$ty.way1.gsmr
		awk -v t=$tx 'NR==1 || $2==t' $tx.$ty.gsmr > $tx.$ty.way2.gsmr
	done
done


fgrep Error *log # make sure no error
awk 'NR==1 || FNR>1' *.way1.gsmr > way1.gsmr
awk 'NR==1 || FNR>1' *.way2.gsmr > way2.gsmr
echo -e "
library(corrplot); library(reshape2)
way1 <- read.table('D:/analysis/gsmr.pheweb/way1.gsmr', header=T, as.is=T)
way2 <- read.table('D:/analysis/gsmr.pheweb/way2.gsmr', header=T, as.is=T)
b1cor <- acast(way1, Outcome ~ Exposure, value.var='bxy'); b1cor[is.na(b1cor)] =0
b2cor <- acast(way2, Exposure ~ Outcome, value.var='bxy'); b2cor[is.na(b2cor)] =0
p1cor <- acast(way1, Outcome ~ Exposure, value.var='p'); p1cor[is.na(p1cor)] =1
p2cor <- acast(way2, Exposure ~ Outcome, value.var='p'); p2cor[is.na(p2cor)] =1
pdf('All-MR.corrplot.pdf')
par(mfrow=c(2,1), mai=c(1,1,0.5,0.5))
corrplot(b1cor, is.corr=F, method='color', type='full', addCoef.col='black', number.cex=0.7, p.mat=p1cor, sig.level=1e-03, insig="pch", pch.col="green", pch.cex=2, tl.col="black", tl.srt=45, outline=T)
corrplot(b2cor, is.corr=F, method='color', type='full', addCoef.col='black', number.cex=0.7, p.mat=p2cor, sig.level=1e-03, insig="pch", pch.col="blue", pch.cex=2, tl.col="black", tl.srt=45, outline=T) # add=T
dev.off()
" > 001.corrplot.R
Rscript 001.corrplot.R


## FUSION TWAS ##
## 安装 TWAS (http://gusevlab.org/projects/fusion/)
dir_tw=/mnt/d/data/twas_data
dir_gt=$dir_tw/GTEx_v7_multi_tissue
dir_ld=$dir_tw/LDREF
fusion=/mnt/d/software_lin/fusion_twas
for trait in RHR T2D; do
for tissue in `ls -d1 $dir_gt/*/ | sed 's/\/$//' | awk -F '/' '{print $NF}' | awk '{printf " "$1}'`; do
for chr in 7; do
    echo now process trait $trait, tissue $tissue, chr $chr
    Rscript $fusion/FUSION.assoc_test.R --sumstats $dir/summary/$trait.sumstats.gz --chr $chr --out $trait.$tissue.chr$chr.txt --weights $dir_gt/$tissue.P01.pos --weights_dir $dir_gt --ref_ld_chr $dir_ld/1000G.EUR.
done
done
done
