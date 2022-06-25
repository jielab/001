
#!/bin/bash
# module avail python3; module load python3/3.8.10
dir=/restricted/projectnb/ukbiobank/jiehuang
dir=/mnt/d


##############################################
## Download and format GWAS
##############################################
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
#	cd $outdir
#	qsub -P ukbiobank -l mem_per_core=1G -pe omp 16 -N $pheno_code.wget -o $pheno_code.LOG -e $pheno_code.ERR < $pheno_code.cmd
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
## for UKB 35 traits
for dat in `ls *array.gz | sed 's/\.array.gz//g'`; do
	echo "#!/bin/bash
	zcat $dat.array.gz | awk '{if (\$8<1e-300) \$8=1e-300; if (NR==1) print \"SNP CHRPOS CHR POS EA NEA EAF N BETA SE Z P\"; else print \$3, \$1\":\"\$2,\$1,\$2, toupper(\$4),toupper(\$5),\"NA\",500000, \$6,\$7,\$6/\$7,\$8}' > $dat.gwas
	zcat $dat.rsid.gz | awk 'NR >1 && $15>=0.4 {if ($8<1e-300) $8=1e-300; print \$NF,\$1\":\"\$2,\$1,\$2, \$5,\$4,$14,500000, \$6,\$7,\$6/\$7,\$8}' >> $dat.gwas
	gzip -f $dat.gwas
	" > $dat.cmd
	chmod 777 $dat.cmd; ./$dat.cmd 
done



##############################################
## Cross-traits meta-analysis
## ##############################################
# MTAG
#conda create -n py2 python=2.7; conda activate py2
conda env create --file environment.yml
# hypocoloc
library(hyprcoloc)
betas <- hyprcoloc::test.betas
ses <- hyprcoloc::test.ses
trait.cor <- hyprcoloc::test.corr
ld.matrix <- hyprcoloc::test.ld
traits <- paste0('T', 1:10)
rsid <- rownames(betas)
sample.overlap = matrix(1, dim(betas)[2], dim(betas)[2]); 
binary.traits = c(1,1,1,rep(0,dim(betas)[2]-3));
# trait.subset = c('T1','T2'), binary.outcomes = binary.traits, 
res <- hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid, snpscores = TRUE, trait.cor = trait.cor, ld.matrix = ld.matrix, sample.overlap = sample.overlap, uniform.priors = FALSE);
cred.sets(res, value = 0.95);
# CM plots of many GWAS
echo "library(tidyverse); library(CMplot)" > plot.R
for t in $traits; do
	echo "$t <- read.table('$t.gwas.p3', header=T) %>% select(SNP, CHR, POS, P) %>% filter(!is.na(SNP)) %>% rename(CHR.$t=CHR, POS.$t=POS, P.$t =P) " >> plot.R
done
echo "
dat <- Reduce(function(x,y) merge(x,y,by='SNP',all=T), list($traits_c)) 
dat1 <- dat %>%
	mutate (chrom=apply( dat[,grepl('CHR',names(dat))], 1, FUN=min, na.rm=T ), position=apply( dat[,grepl('POS',names(dat))], 1, FUN=min, na.rm=T )) %>% 
	select(SNP, chrom, position, $traits_cp)
write.table(dat1, 'gwas.merged.txt', append=F, quote=F, row.names=F, col.names=T)
CMplot(dat1, plot.type='m',multracks=TRUE, threshold=c(5e-8,1e-6), threshold.lty=c(1,2), threshold.lwd=c(1,1), threshold.col=c('black','grey'), amplify=TRUE,bin.size=1e6,
	chr.den.col=c('darkgreen', 'yellow', 'red'), signal.col=c('red','green','blue'), signal.cex=1, file='jpg',memo='',dpi=300,file.output=TRUE,verbose=TRUE)
" >> plot.R



##############################################
## LDSC & HDL
##############################################
conda activate ldsc
ldscdir=/mnt/d/data/ldsc
ldrefdir=$ldscdir/eur_w_ld_chr # $ldscdir/weights_hm3_no_hla
trait=covid_resf; traits="covid_resf covid_icu covid_r5 covid_r6"
traits_c=`echo $traits | sed 's/ /,/g'`; traits_cp=`echo $traits_c |  sed 's/^/P./; s/,/,P./g'`
ldsc.py --h2 $trait.sumstats.gz --ref-ld-chr $ldrefdir/ --w-ld-chr $ldrefdir/
#ldsc.py --h2 $trait.sumstats.gz --ref-ld-chr $ldrefdir/CNS.,baseline. --w-ld-chr $ldrefdir/weights. --overlap-annot --frqfile-chr $datadir/1000G.mac5eur. --out $trait --print-coefficients
echo $traits | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % ldsc.py --rg % --out $trait.rg --ref-ld-chr $ldrefdir/ --w-ld-chr $ldrefdir/
awk '$1=='Summary' {printf NR}' $trait.rg.log | xargs -n1 -I % awk -v s=% 'FNR >=s' $trait.rg.log | sed 's/.sumstats.gz//g' > $trait.rg.txt
#Rscript /mnt/d/software_lin/HDL/HDL.run.R gwas1.df=covid_res.hdl.txt gwas2.df=covid_icu.hdl.txt LD.path=/mnt/d/data/ldsc/UKB_imputed_SVD_eigen99_extraction output.file=jie.Rout
library(HDL)
dat1 <- read.table('covid_res.hdl.txt', header=T)
dat2 <- read.table('covid_icu.hdl.txt', header=T)
LD.path <- 'D:/data/ldsc/UKB_imputed_SVD_eigen99_extraction'
res.HDL <- HDL.rg(dat1, dat2, LD.path)



##############################################
## Mendelian randomization  
##############################################
dir=/mnt/d
gdir=$dir/data/gwas/pheweb
outdir=$dir/analysis/mr
mr_ex_list=$dir/files/mr.ex.list
mr_ou_list=$dir/files/mr.ou.list
cut -f 1,3 $mr_ex_list | sed '1d' | while IFS=$'\t' read -r ex_name ex_ieu; do
	outdir2=$outdir/$ex_name; mkdir -p $outdir2
	ex_top_snps=$gdir/$ex_name.top.snps
	echo -e "library(CMplot); library(tidyverse)
	dat <- read.table('$gdir/$ex_name.gwas.p3', header=T, as.is=T)
	dat\$P[dat\$P<1e-100] <- 1E-100
	top <- read.table('$gdir/$ex_name.top.snps', header=T, as.is=T)
	dat1 <- dat %>% select(SNP, CHR, POS, P)
	CMplot(dat1, plot.type='m', col=c('grey30','grey60'), highlight=top[,1], highlight.col='green', highlight.cex=1, highlight.pch=19, file='jpg', memo='', chr.border=T, dpi=300, file.output=T, verbose=T, width=12, height=8)
	" > $outdir2/$ex_name.mh.R
	echo -e "#!/bin/bash
	Rscript $ex_name.mh.R
	mv Rectangular-Manhattan.P.jpg $ex_name.mh.jpg
	zcat $gdir/$ex_name.gwas.gz | fgrep -wf $gdir/$ex_name.top.snps > $ex_name.top.txt
	" > $outdir2/$ex_name.cmd
	cd $outdir2; chmod 777 $ex_name.cmd; ./$ex_name.cmd 
	cut -f 1,3 $mr_ou_list | sed '1d' | while IFS=$'\t' read -r ou_name ou_ieu; do
		if [[ $ex_name == $ou_name ]]; then continue; fi
		echo -e "
	#	source('$dir/scripts/library/mr2s.f.R')
	#	mr2sfn( top_snps_file='$top_snps_file', ex_name='$ex_name', ex_ieu='$ex_ieu', ex_gwas='$ex_gwas_file', ou_name='$ou_name', ou_ieu='$ou_ieu', ou_gwas='$ou_gwas')
		source('$dir/scripts/library/compareB.f.R')
		compareB( 
			f1='$ex_name.top.txt', f1_name='$ex_name', f1_snp='^snp\$|variant_id', f1_ea='^ea\$|^alt\$|^effect_allele\$', f1_nea='^nea\$|^ref\$|other_allele', f1_eaf='^eaf\$|effect_allele_frequency', f1_beta='^beta\$', f1_se='^se\$|standard_error', f1_p='^p\$|^pval\$|p_value',
			f2='$ex_name.$ou_name.top.txt', f2_name='$ou_name', f2_snp='^snp\$|variant_id', f2_ea='^ea\$|^alt\$|^effect_allele\$', f2_nea='^nea\$|^ref\$|other_allele', f2_eaf='^eaf\$|effect_allele_frequency', f2_beta='^beta\$', f2_se='^se\$|standard_error', f2_p='^p\$|^pval\$|p_value'
		)
		" > $outdir2/$ex_name.$ou_name.mr.R
		ou_gwas_file=$gdir/$ou_name.gwas.gz
		echo -e "#!/bin/bash
		zcat $ou_gwas_file | fgrep -wf $gdir/$ex_name.top.snps > $ex_name.$ou_name.top.txt
		nrow=\`wc -l < $ex_name.$ou_name.top.txt\`
		if [[ \$nrow >1 ]]; then 
			Rscript $ex_name.$ou_name.mr.R
		fi
		" > $outdir2/$ex_name.$ou_name.cmd
		cd $outdir2; chmod 777 $ex_name.$ou_name.cmd; ./$ex_name.$ou_name.cmd 
	done
done
#cat */*.mr2s.res.txt | sed 's/ /-/g; s/|/\t/g' | awk 'NF==6' > 001.mr2s.txt
#pdftk */*.cB.pdf cat output 001.cB.pdf
cat */*.mr.res.txt | awk 'NF ==10' > 001.mr.txt
## run in R
pacman::p_load(data.table, dplyr, tidyverse, magrittr, corrplot, reshape2)
dat <- read.table('001.mr.txt', header=F, as.is=T)
names(dat) <- c('ex_name', 'ou_name', 'cnt', 'p.ivw', 'p.median', 'p.maxlik', 'p.mbe', 'p.conmix', 'pleio.egger', 'p.egger')
#dat1 <- dat %>% group_by(ex_name, ou_name) %>% summarize(Pmin=min(p)) %>% mutate(p= -log10(Pmin))
dat1 <- dat %>% select(ex_name, ou_name, p.ivw) %>% mutate (p=-log10(p.ivw))
pmat <- acast(dat1, ou_name ~ ex_name, value.var='p')
pmat[is.na(pmat)] =0;  pmat=round(pmat,1); pmat[pmat >7] =7; 
pdf('001.mr.pval.pdf', w=20, h=20)
#png('jie.png', w=4000, h=4000, res=512)
corrplot(pmat, is.corr=F, method='shade', bg='black', col=colorRampPalette(c("white","green","gold"))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig="pch", pch.cex=2, tl.srt=45, outline=T)
dev.off()
#dat <- read.table('my-summary.gsmr', header=T, as.is=T)
#beta <- acast(dat, Outcome ~ Exposure, value.var='bxy'); b1cor[is.na(b1cor)] =0
#pval <- acast(dat, Outcome ~ Exposure, value.var='p'); p1cor[is.na(p1cor)] =1
#corrplot(beta, is.corr=F, method='color', type='full', addCoef.col='black', number.cex=0.7, p.mat=pval, sig.level=1e-03, insig="pch", pch.col="green", pch.cex=2, tl.col="black", tl.srt=45, outline=T)
