Arr1=("snp" "chr" "pos" "ea" "nea" "eaf" "n" "beta" "se" "p")
Arr2=("snp|rsid|variant_id" "chr|chrom|chromosome" "pos|bp|base_pair" "ea|alt|a1|eff.allele" "nea|ref|a2|ref.allele" "eaf" "n" "beta" "se|standard_error" "^p|pval")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C0: GWAS download and format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir=/mnt/d
outdir=$dir/data/gwas/pheweb2
files=$dir/files/mr.list.txt
cut -f 1,4 $files | sed '1d' | fgrep -vw "NA" | while IFS=$'\t' read -r pheno_code file_url; do
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
for dat in `ls *gz | sed 's/\.gz$//g'`; do
	head_row=`zcat $dat.gz | head -1 | sed 's/\t/ /g'`;
	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
	done
	echo $snp $chr $pos $ea $nea $eaf $n $beta $se $p
	echo "#!/bin/bash
	zcat $dat.gz | awk '{if ($p_col<1e-300) $p_col=1e-300; if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print $snp,$chr,$pos, toupper($ea), toupper($nea),$eaf,$n, $beta,$se,$p}' | sort -k 2,2V -k 3,3n | gzip -f > $dat.gz
	" > $dat.cmd
	chmod 777 $dat.cmd; ./$dat.cmd 
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1: genetic correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
conda activate ldsc # 采用 python2，而不是 python3
dir=/mnt/d/data/gwas/penguin
ldscdir=/mnt/d/data/ldsc
dats=`cd $dir; ls *.gz | sed 's/\.gz//g'`
for dat in $dats; do
	echo $dat
	datf=$dir/$dat.gz
	head_row=`zcat $datf | head -1 | sed 's/\t/ /g'`;
	munge_sumstats.py --chunksize 10000 --sumstats $dat.tmp.gz --merge-alleles $ldscdir/w_hm3.snplist --out $dat --snp SNP --a1 EA --a2 NEA --ignore SNPID,OR --signed-sumstats BETA,0 --p P --N 100000 
done
for dat in $dats; do
	echo $dat $dats | sed -e 's/ /.sumstats.gz,/g' -e 's/$/.sumstats.gz/' | xargs -n1 -I % ldsc.py --rg % --out $dat.rg --ref-ld-chr $ldscdir --w-ld-chr $ldscdir
	beginl=`awk '$1=="Summary" {printf NR}' $dat.rg.log` 
	awk -v s=$beginl 'FNR >s' $dat.rg.log | head -n -3 | sed 's/.sumstats.gz//g'  > $dat.rg.txt
done
awk 'NR==1 || FNR>1' *.rg.txt > all.rg.res
#下面是R代码
pacman::p_load(tidyverse, reshape2, ggplot2, ggcorrplot)
dat <- read.table('D:/analysis/ldsc/all.rg.res', header=T)
rg <- dat %>% select(p1, p2, rg) %>% acast(p1 ~ p2, value.var='rg'); rg[is.na(rg)] =0;  rg=round(rg,1)
pval <- dat %>% select(p1, p2, p) %>% acast(pval, p1 ~ p2, value.var='p')
plt <- ggcorrplot(rg, lab=T, p.mat=pval, sig.level=5e-4, insig ='blank') 
	plt + theme(axis.title=element_text(size=15, face='bold'), axis.text=element_text(size=12, face='bold'))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2: Mendelian Randomization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exp_dir=/mnt/d/data/gwas/mb
out_dir=/mnt/d/data/gwas/penguin
dats_x=`cd $exp_dir; ls *.gz | sed 's/\.gz//g'`
dats_y=`cd $out_dir; ls *.gz | sed 's/\.gz//g'`
for exp in $dats_x; do
	exp_ieu=NA; exp_file0=$exp_dir/$exp.gz; exp_pheno=$exp
	exp_iv_file=$exp_dir/$exp.top.snp
	if [ ! -e $exp_iv_file ]; then
		exp_iv_file=$exp.top.snp; zcat $exp_file0 | awk 'NR==1 || $10<=1e-6 {b=sprintf("%.0f",$3/1e5); print $4,$2,$3,$10,b}' | sort -k 2,2n -k 5,5n -k 4,4g | awk '{if (arr[$NF] !="Y") print $1; arr[$NF] ="Y"}' > $exp.top.snp
	fi
	exp_iv_line=`wc -l $exp_iv_file | awk '{printf $1}'`
	if [ "$exp_iv_line" -le 1 ]; then
		echo $exp has no IV SNP
		continue 
	fi
	head_row=`zcat $exp_file0 | head -1 | sed 's/\t/ /g'`
	for i in ${!Arr1[@]}; do
		eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
	done
	echo exp $exp, snp $snp, ea $ea, nea $nea, beta $beta, se $se, p $p
	exp_file=$exp.txt; zcat $exp_file0 | fgrep -wf $exp_iv_file | awk -v snp=$snp -v ea=$ea -v nea=$nea -v beta=$beta -v se=$se -v p=$p BEGIN'{print "SNP EA NEA BETA SE P"}{print $snp, toupper($ea), toupper($nea), $beta, $se, $p}' > $exp.txt	
	for out in $dats_y; do
		if [ "$exp" = "$out" ]; then continue; fi
		echo Running: $exp $out 
		out_ieu=NA; out_file0=$out_dir/$out.gz; out_pheno=$out
		head_row=`zcat $out_file0 | head -1 | sed 's/\t/ /g'`
		for i in ${!Arr1[@]}; do
			eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`
		done
		echo out $out, snp $snp, ea $ea, nea $nea, beta $beta, se $se, p $p
		out_file=$exp.list.$out.txt; zcat $out_file0 | fgrep -wf $exp_iv_file | awk -v snp=$snp -v ea=$ea -v nea=$nea -v beta=$beta -v se=$se -v p=$p BEGIN'{print "SNP EA NEA BETA SE P"}{print $snp, toupper($ea), toupper($nea), $beta, $se, $p}' > $exp.list.$out.txt
		echo "
		source('/mnt/d/scripts/library/tsmr.f.R')
		tsmr_fn( exp_iv_file='$exp_iv_file',
			exp_name='$exp', exp_pheno='$exp_pheno', exp_ieu='$exp_ieu', exp_file='$exp_file',  
			out_name='$out', out_pheno='$out_pheno', out_ieu='$out_ieu', out_file='$out_file'
		)
		" > $exp.$out.R
		R CMD BATCH $exp.$out.R
	done
done
cat *.tsmr.out | sed 's/|/\t/g' > ../all.mr.res
#下面是R代码
# pacman::p_load(tidyverse, reshape2, ggplot2, corrplot, ggcorrplot)
# exp_dat <- read_exposure_data(filename="D:/analysis/stoolf.txt", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
# out_dat <- read_outcome_data(filename="D:/analysis/hdl.txt", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
# dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat)
# mr(dat); run_mr_presso(dat)
dat <- read.table('D:/analysis/all.mr.res2', sep='\t', header=F, as.is=T) %>% subset(method %in% c("Wald ratio", "Inverse variance weighted"))
	names(dat) <- c('exp_name', 'exp_cnt', 'out_name', 'out_cnt', 'dat_cnt', 'method', 'beta', 'se', 'p')
beta <- dat %>% select(exp_name, out_name, beta) %>% acast(beta, exp_name ~ out_name, value.var='beta'); beta[is.na(beta)] =0; beta=round(beta,1)
pval <- dat %>% select(exp_name, out_name, p) %>% acast(pval, exp_name ~ out_name, value.var='p'); pval[is.na(pval)] =1
plt <- ggcorrplot(beta[,151:211], lab=T, p.mat=pval[,151:211], sig.level=.25e-4, insig ='blank') 
	plt + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
#corrplot(beta, is.corr=F, method='shade', bg='black', col=colorRampPalette(c('white','green','gold'))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig='pch', pch.cex=2, tl.srt=45, outline=T)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C3: Colocalization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cut -f 1-3 x.education > ref.chrpos; cut -f 1 x.education > ref.snp
for dat in `ls *.joined | sed 's/\.joined//g'`; do
	echo $dat
	#python /mnt/d/scripts/library/join_file.py -i "ref.snp,TAB,0 $dat,TAB,0" -o $dat.joined
	cut -d " " -f 1,5,6,9,10 $dat.joined | sed "1 s/ /\.$dat /g; 1 s/$/\.$dat/" > $dat.5col
done
paste ref.chrpos *.5col > all.coloc.tmp
awk '{print NF}' all.coloc.tmp | sort -nu
num=`ls -l *.5col | wc -l | awk '{printf $1}'`
awk -v num=$num '{for (i=2;i<=num;i++) {if (NR>1 && $(i*5-1) !=$1) print $1, $i, "ERR"; else if (NR>1 && $(i*5) !=$5 && $(i*5+2) !="NA") $(i*5+2) =-$(i*5+2); $(i*5)=""; $(i*5+1)="" }; print $0}' all.coloc.tmp > all.coloc.tmp2; fgrep ERR all.coloc.tmp2
cat all.coloc.tmp2 | awk 'NR==1 || $10<=1e-6 {b=sprintf("%.0f",$3/1e5); print $4,$2,$3,$10,b}' | sort -k 2,2n -k 5,5n -k 4,4g | awk '{if (arr[$NF] !="Y") print $1; arr[$NF] ="Y"}' > all.coloc.txt
#下面是R代码
library(hyprcoloc)
dat <- read.table("D:/analysis/coloc/all.coloc.txt", header=T, as.is=T) %>% na.omit()
betas <- subset(dat, select=grepl("BETA", names(dat))); betas <- as.matrix(betas)
ses <- subset(dat, select=grepl("^SE", names(dat))); ses <- as.matrix(ses)
traits <- grep("x\\.|y\\.", names(dat), value=T)
rsid <- as.matrix(dat[,1])
hyprcoloc(betas, ses, trait.names=traits, trait.subset=c("x.education", "y.dementia", "y.depress"), snp.id=rsid)