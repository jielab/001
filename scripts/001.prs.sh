#!/bin/bash -l
Arr1=("snp" "chrom" "pos" "ea" "nea" "eaf" "n" "beta" "se" "p")
Arr2=("snp|rsid|variant_id" "chr|chrom|chromosome" "pos|bp|base_pair" "alt|ea|eff.allele|effect_allele|a1|allele1|PGS_EFFECT_ALLELE" "ref|nea|ref.allele|other_allele|a2|allele0|PGS_OTHER_ALLELE" "eaf|a1freq|effect_allele_freq|PGS_WEIGHT" "n" "beta" "se|standard_error|PGS_POSTERIOR_STANDARD_ERROR" "p|pval|p_bolt_lmm")


dir=/work/sph-huangj
gendir=/data/sph-huangj/data/ukb/gen/imp
gwasdir=$dir/files/posCtrls
for label in gallstone; do
	echo RUN $label
	gwas=$gwasdir/$label.ref
	outdir=$dir/analysis/prs/$label; if [ -d "$outdir" ]; then echo $label already run; continue; fi; mkdir $outdir
	head_row=`cat $gwas | head -1 | sed 's/\t/ /g'`
	snp=""; chrom=""; ea=""; beta=""
	for i in ${!Arr1[@]}; do eval ${Arr1[$i]}=`echo $head_row | tr ' ' '\n' | grep -Einw ${Arr2[$i]} | sed 's/:.*//'`; done
	echo snp $snp, chr $chrom, ea $ea, beta $beta
	cols="$snp $ea $beta header"
	echo "#!/bin/bash
		for chr in {1..22} X; do
			cat $gwas | awk -v chr=\$chr 'NR==1 || \$$chrom==chr' > chr\$chr.ref
			nref=\`wc -l chr\$chr.ref | awk '{printf \$1}'\`
			if [[ \$nref == 1 ]]; then continue; fi
			plink2 --pfile $gendir/chr\$chr --chr \$chr --score chr\$chr.ref $cols no-mean-imputation ignore-dup-ids cols=+scoresums list-variants --out chr\$chr
		done
		paste -d ' ' chr*.sscore > $label.prs.tmp
		awk '{print NF}' $label.prs.tmp | sort -nu
		num=\`ls -l chr*.sscore | wc -l | awk '{printf \$1}'\`
		awk -v num=\$num '{if (NR==1) print \"eid $label.allele_cnt $label.score_sum\"; else {allele_cnt=score_sum=0; for (i=1;i<=num;i++) {if (\$(i*6-5) !=\$1) score_sum=score_sum\"\"i\"ERR,\"; else {allele_cnt=allele_cnt+\$(i*6-3); score_sum=score_sum+\$(i*6-0)}}; print \$1, allele_cnt, score_sum} }' $label.prs.tmp > $label.prs.txt
		fgrep ERR $label.prs.txt
	" > $outdir/prs.cmd
	cd $outdir
	bsub -q smp -J $label.prs -o prs.LOG -e prs.ERR < prs.cmd
done

#
echo "
fgrep \"error|ERR\" */*.log 
paste -d ' ' */*.prs.txt > all.prs.txt
awk '{print NF}' all.prs.txt | sort -nu
num=`ls -l */*prs.txt | wc -l | awk '{printf \$1}'`
awk -v num=\$num '{for (i=1;i<=num;i++) {if (\$(i*3-2) !=\$1) print \$1, \$i, \"ERR\" }}' all.prs.txt | head
gzip all.prs.txt
" > prs.step2.cmd

