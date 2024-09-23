#!/bin/bash

dir=/work/sph-huangj
pathX=ppp
pathY=main

for Y in bald12 bald13 bald14 bald134 bald1234; do
	outdir=$dir/analysis/assoc.sum/$Y; mkdir -p $outdir
	for X in `cd $dir/data/gwas/$pathX/clean; ls *gz | sed 's/\.gz//g'`; do 
	# height ABO GLP1R walk_pace x.age_fbirth x.age_fsex x.bmi x.crp x.drink x.education x.smoke x.vitad; do
	# id.1850 id.1853 id.400 id.419 id.432 id.433 id.436
	if [ $X == $Y ]; then continue; fi
	for pathM in NA; do # bb bc inf mb mb2 met ppp
	label="$pathX-$pathM-$pathY"
	if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi

	echo "#!/bin/bash
	module load python/anaconda3/2020.7 # 调用集群conda环境
	source activate # 进入conda的base环境
	conda activate R
	cat $dir/scripts/main/assoc.sum.R | sed '9 s/?/$pathX/; 10 s/?/$pathY/; 11 s/?/$pathM/; 12 s/?/$X/; 13 s/?/$Y/; 14 s/?/$label/g' > $label.R
	R CMD BATCH $label.R
	" > $outdir/$label.cmd
	cd $outdir
	bsub -q short -n 40 -J sum.$label -o $label.LOG -e $label.ERR < $label.cmd
done
done
done

# awk 'FNR==1 && $1 !~ /SKIP/' */*-met-*.log | cut -d ' ' -f 2,4,6-8 | sort -k 5,5g > mr.top.txt
# echo "X M Y iv X2Y.b X2Y.p TOT.b TOT.p ACME.b ACME.p ADE.b ADE.p Prop.b Prop.p" > mrMed.top.txt; 
# cat */*.log | cut -d ' ' -f 2-4,6-8,10-17 | awk 'NF==14 && $9*$11 >0 && $8<0.01 && $10<0.01 && $14<0.01 && $13>=0.1' | sed 's/"//g' | grep -v "^id" | sort -k 13,13gr >> mrMed.top.txt

