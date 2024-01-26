#!/bin/bash -l


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate PRS by PRS-CSx
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

python PRScs.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --out_dir=OUTPUT_DIR [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --chrom=CHROM --beta_std=BETA_STD --write_psi=WRITE_PSI --seed=SEED]
python PRScsx.py --ref_dir=PATH_TO_REFERENCE --bim_prefix=VALIDATION_BIM_PREFIX --sst_file=SUM_STATS_FILE --n_gwas=GWAS_SAMPLE_SIZE --pop=POPULATION --out_dir=OUTPUT_DIR --out_name=OUTPUT_FILE_PREFIX [--a=PARAM_A --b=PARAM_B --phi=PARAM_PHI --n_iter=MCMC_ITERATIONS --n_burnin=MCMC_BURNIN --thin=MCMC_THINNING_FACTOR --chrom=CHROM --meta=META_FLAG --seed=SEED]

--sst_file=path_to_sumstats/sumstats_se.txt 
--sst_file=path_to_sumstats/EUR_sumstats.txt,EAS_sumstats.txt 


python PRScsx.py --ref_dir=path_to_ref --bim_prefix=path_to_bim/test --sst_file=path_to_sumstats/EUR_sumstats.txt,path_to_sumstats/EAS_sumstats.txt \
	--n_gwas=200000,100000 --pop=EUR,EAS --chrom=22 --phi=1e-2 --out_dir=path_to_output --out_name=test


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Calculate PRS by PLINK
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash -l
dir=/work/sph-huangj
gendir=/data/sph-huangj/ukb/imp
gwasdir=$dir/files/posCtrls
for label in `cd $gwasdir; ls *ref | sed 's/\.ref$//g'`; do
	echo RUN $label
	gwas=$gwasdir/$label.ref
	outdir=$dir/tmp/$label
		if [ -d $outdir ]; then rm -r $outdir; fi; mkdir -p $outdir
	head_row=`head -1 $gwas | sed 's/\t/ /g'`
		chr_str=`echo $head_row | tr ' ' '\n' | grep -Einw 'chr|chrom|chromosome'`; chr_col=`echo $chr_str | sed 's/:.*//'`
		snp_str=`echo $head_row | tr ' ' '\n' | grep -Einw 'snp|rsid|variant_id|markername'`; snp_col=`echo $snp_str | sed 's/:.*//'`
		ea_str=`echo $head_row  | tr ' ' '\n' | grep -Einw 'ea|alt|effect_allele|allele1'`; ea_col=`echo $ea_str | sed 's/:.*//'`
		beta_str=`echo $head_row | tr ' ' '\n' | grep -Einw 'beta'`; beta_col=`echo $beta_str | sed 's/:.*//'`
		cols="$snp_col $ea_col $beta_col header"; echo $cols
	for chr in {1..22}; do
		echo "#!/bin/bash
		awk 'NR==1 || \$$chr_col==$chr' $gwas > chr$chr.ref
		nref=\`wc -l chr$chr.ref | awk '{printf \$1}'\`
		if [[ \$nref == 1 ]]; then exit; fi
		plink2 --pfile $gendir/chr$chr --chr $chr --score chr$chr.ref $cols no-mean-imputation ignore-dup-ids cols=+scoresums list-variants --out chr$chr
		" > $outdir/chr$chr.cmd
		cd $outdir
		bsub -q smp -J $label.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
	done
	echo "#!/bin/bash
		paste -d ' ' chr*.sscore > $label.prs.tmp
		awk '{print NF}' $label.prs.tmp | sort -nu
		num=\`ls -l chr*.sscore | wc -l | awk '{printf \$1}'\`
		awk -v num=\$num '{if (NR==1) print \"eid $label.allele_cnt $label.score_sum\"; else {allele_cnt=score_sum=0; for (i=1;i<=num;i++) {if (\$(i*6-5) !=\$1) score_sum=score_sum\"\"i\"ERR,\"; else {allele_cnt=allele_cnt+\$(i*6-3); score_sum=score_sum+\$(i*6-0)}}; print \$1, allele_cnt, score_sum} }' $label.prs.tmp > $label.prs.txt
		fgrep ERR $label.prs.txt
	" > $outdir/step2.cmd
	cd $outdir
	bsub -q smp -J $label.step2 -w "ended($label.chr*)" -o step2.LOG -e step2.ERR < step2.cmd
done
# Merge many PRS files into one
fgrep "error|ERR" */*.log 
paste -d ' ' */*.prs.txt > all.prs.txt
awk '{print NF}' all.prs.txt | sort -nu
num=`ls -l */*prs.txt | wc -l | awk '{printf $1}'`
awk -v num=$num '{for (i=1;i<=num;i++) {if ($(i*3-2) !=$1) print $1, $i, "ERR" }}' all.prs.txt | head
gzip all.prs.txt