#!/bin/bash -l

threads=40
dir=/work/sph-huangj
impdir=$dir/data/ukb/gen/imp
label=abq; traits="a012 b012"; binary=N
#label=walk; traits="walk_brisk"; binary=Y
pheno=$dir/data/ukb/phe/common/ukb.pheno; covar=$pheno

if [[ $binary == Y ]]; then 
	col_str="+a1freqcc,+a1countcc,+totallelecc"; phe_str="--1"; postfix="logistic.hybrid"
else
	col_str="+a1freq,+a1count,+totallelecc"; phe_str="--1"; postfix="linear" ## phe_str="--pheno-quantile-normalize $traits"
fi

if [ ! -f $impdir/sp_grm ]; then
	for par in 2 3; do
	echo "gcta --mpfile $impdir/mpfile.list --make-grm-part 3 $par --sparse-cutoff 0.05 --thread-num $threads --out grm.par$par" > $impdir/grm.par$par.cmd
	cd $impdir
#	bsub -q smp -n $threads -J grm.par$par -o grm.par$par.LOG -e grm.par$par.ERR < grm.par$par.cmd
	done
fi


## 开始运行 GWAS 🏮🏮
for sex in a; do # m f

	echo processing $label $sex
	if [[ $sex == m ]]; then 
		sex_str="--keep-males"; sex_cov=""; sex_post=".$sex"
	elif [[ $sex == f ]]; then 
		sex_str="--keep-females"; sex_cov=""; sex_post=".$sex"
	else 
		sex_str=""; sex_cov=",sex"; sex_post=""
	fi
	outdir=$dir/data/gwas/self/$label$sex_post; mkdir -p $outdir
	for chr in {1..22} X; do
		x_str=""; if [[ $chr == X ]]; then x_str="--xchr-model 2"; fi
		echo "#!/bin/bash
		if [ ! -f $impdir/chr$chr.extract ]; then
			plink2 --memory 12000 --threads $threads --pfile $impdir/chr$chr --freq counts --out $impdir/chr$chr
			awk '\$5>100 {print \$2}' $impdir/chr$chr.acount | sed '1 s/ID/SNP/' > $impdir/chr$chr.extract
		fi
		plink2 --memory 12000 --threads $threads --pfile $impdir/chr$chr --extract $impdir/chr$chr.extract --glm cols=+omitted,+nobs,+beta,$col_str hide-covar allow-no-covars no-x-sex --pheno $pheno --no-psam-pheno --pheno-name $traits $phe_str $sex_str $x_str --covar $covar --covar-name age$sex_cov,PC1-PC4 --out chr$chr
		# gcta --pfile chr$chr --grm-sparse sp_grm --joint-covar --fastGWA-mlm-binary --pheno phenotype.txt --qcovar pc.txt --covar fixed.txt --thread-num 10 --out geno_assoc
			
		Arr1=(\"SNP\" \"CHR\" \"POS\" \"EA\" \"NEA\" \"EAF\" \"N\" \"BETA\" \"SE\" \"P\")
		Arr2=(\"snp|rsid|id|variant_id\" \"chr|chrom|chromosome\" \"pos|bp|base_pair\" \"alt|ea|eff.allele|effect_allele|allele1\" \"ref|nea|ref.allele|other_allele|allele0\" \"eaf|a1_freq|effect_allele_freq\" \"n|obs_ct\" \"beta\" \"se|standard_error\" \"p|pval|p_bolt_lmm\") # avoid multiple matches
		for t in $traits; do 
			sed -i '1 s/^#//; 1 s/?//g' chr$chr.\$t.glm.$postfix
			head_row=\`head -1 chr$chr.\$t.glm.$postfix | sed 's/\t/ /g'\`
			for i in \${!Arr1[@]}; do
				eval \${Arr1[\$i]}=\`echo \$head_row | tr ' ' '\n' | grep -Einw \${Arr2[\$i]} | sed 's/\w*://'\`
				eval \${Arr1[\$i]}_col=\`echo \$head_row | tr ' ' '\n' | grep -Einw \${Arr2[\$i]} | sed 's/:\w*//'\`
			done
			cat chr$chr.\$t.glm.$postfix | awk -v snp=\$SNP_col -v chr=\$CHR_col -v pos=\$POS_col -v ea=\$EA_col -v nea=\$NEA_col -v eaf=\$EAF_col -v n=\$N_col -v beta=\$BETA_col -v se=\$SE_col -v p=\$P_col '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print \$snp,\$chr,\$pos, \$ea,\$nea,\$eaf,\$n, \$beta,\$se,\$p}' | sed 's/ /\t/g' | sort -k 2,2n -k 3,3n | gzip -f > chr$chr.\$t.glm.$postfix.gz
		#	rm  chr$chr.\$t.glm.$postfix
		done
		" > $outdir/chr$chr.cmd 
		cd $outdir # -q smp -n 10 -R "span[ptile=192]"
	#	bsub -q smp -n 40 -J $label$sex_post.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
	done
	
	## 看GWAS是否顺利跑完🏮🏮
	chr_col=2; pos_col=3
	echo "#!/bin/bash
	for dat in \`ls *.gz | sed 's/\\.gz//g'\`; do
		cols=\`zcat \$dat.gz | awk '{print NF}' | sort -nu | wc -l\`
		if [[ \$cols -gt 1 ]]; then echo \"\$dat has \$cols different number of columns\" >> check.cols; fi
		# zcat \$dat.gz | awk -v filename=\$dat -v chr_col=$chr_col -v pos_col=$pos_col -f check_order.awk >> check.order
	done
	" > $outdir/check.cmd
	cd $outdir
	# bsub -q short -n 40 -J check.$label -o check.LOG -e check.ERR < check.cmd
	
	## 生成 P<1e03 的小文件🏮🏮
	echo "#!/bin/bash
	for trait in $traits; do
		ls -v chr*.\$trait.glm.$postfix.gz | xargs zcat | grep -vwE \"nan|inf|NA\" | awk 'NR==1 || (\$NF != \"NA\" && \$NF !=\"P\")' | gzip -f > ../\$trait$sex_post.gz
		zcat ../\$trait$sex_post.gz | awk 'NR==1 || $NF<1e-03' > ../\$trait$sex_post.1e-3
	done
	" > $outdir/step2.cmd
	# bsub -q smp -n 40 -J $label$sex_post.step2 -w "ended($label$sex_post.chr*)" -o step2.LOG -e step2.ERR < step2.cmd #  -w "ended($label$sex_post.chr*)"

done

