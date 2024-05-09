#!/bin/bash -l

threads=10
dir=/work/sph-huangj
impdir=/data/sph-huangj/ukb/imp
label=death60; traits="death60a death60b death60c"; binary=Y
#label=walk; traits="walk_brisk"; binary=Y
pheno=$dir/data/ukb/phe/common/ukb.pheno; covar=$pheno
if [[ $binary == Y ]]; then 
	phe_str="--1"; postfix="logistic"
else
	phe_str="--pheno-quantile-normalize $traits"; postfix="linear"
fi

sex_str=""; sex_cov=",sex"
for sex in a; do
	echo processing $label $sex
	if [[ $sex == m ]]; then sex_str="--keep-males"; sex_cov=""; fi
	if [[ $sex == f ]]; then sex_str="--keep-females"; sex_cov=""; fi
	outdir=$dir/gwas/$label.$sex; mkdir -p $outdir
	for chr in {1..22} X; do
		x_str=""; if [[ $chr == X ]]; then x_str="--xchr-model 2"; fi
		echo "#!/bin/bash
		if [ ! -f $impdir/chr$chr.extract ]; then
			plink2 --memory 12000 --threads $threads --pfile $impdir/chr$chr --freq counts --out $impdir/chr$chr
			awk '\$5>100 {print \$2}' $impdir/chr$chr.acount | sed '1 s/ID/SNP/' > $impdir/chr$chr.extract
		fi
		plink2 --memory 12000 --threads $threads --pfile /data/sph-huangj/ukb/imp/chr$chr --maf 0.001 --glm cols=+omitted,+a1freq,+nobs,+beta hide-covar allow-no-covars no-x-sex --pheno $pheno --no-psam-pheno --pheno-name $traits $phe_str $sex_str $x_str --covar $covar --covar-name age$sex_cov,PC1-PC4 --out chr$chr
		for t in $traits; do 
			if [[ $binary == Y ]]; then
				cat chr$chr.\$t.glm.logistic.hybrid | awk -v OFS='\t' '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print \$3,\$1,\$2, \$7,\$8,\$9,\$12, \$13,\$14,\$16}' | sed '1 s/ /\t/g' | sort -k 2,2n -k 3,3n | gzip -f > chr$chr.\$t.glm.logistic.gz
				rm  chr$chr.\$t.glm.logistic.hybrid
			else
				cat chr$chr.\$t.glm.linear | awk -v OFS='\t' '{if (NR==1) print \"SNP CHR POS EA NEA EAF N BETA SE P\"; else print \$3,\$1,\$2, \$7,\$8,\$9,\$11, \$12,\$13,\$15}' | sed '1 s/ /\t/g' | sort -k 2,2n -k 3,3n | gzip -f > chr$chr.\$t.glm.linear.gz
                                rm  chr$chr.\$t.glm.linear
			fi
		done
		" > $outdir/chr$chr.cmd 
		cd $outdir #ser smp short
		bsub -q smp -n 10 -R "span[ptile=192]" -J $label.$sex.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
	done
	echo "#!/bin/bash
	for trait in $traits; do
		ls -v chr*.\$trait.glm.$postfix.gz | xargs zcat | awk 'NR==1 || \$8 ~ /[0-9]/' | gzip -f > ../\$trait.$sex.gz
	done
	" > $outdir/step2.cmd
	bsub -q ser -n 1 -J $label.$sex.step2 -w "ended($label.$sex.chr*)" -o step2.LOG -e step2.ERR < step2.cmd #  -w "ended($label.$sex.chr*)"
done

