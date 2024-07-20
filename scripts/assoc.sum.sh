#!/bin/bash

dir=/work/sph-huangj/analysis
pathX=main
pathY=main

for Y in covid_A2; do # covid_B2 covid_C2 y.cad y.cirrhosis y.dementia y.in y.t2dm y.varicose y.vte; do
	outdir=$dir/assoc.sum/$Y; mkdir -p $outdir
	for X in ABO; do # GLP1R walk_pace x.age_fbirth x.age_fsex x.crp x.vitad; do
	# id.1850 id.1853 id.400 id.419 id.432 id.433 id.436
	if [ $X == $Y ]; then continue; fi
	for M in ppp; do # bb bc in mb mb2 met ppp; do
	if [ $M == $X ]; then continue; fi
	label="$X-$M-$Y"
	if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi

	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	cat /work/sph-huangj/scripts/main/assoc.sum.R | sed '8 s/?/$pathX/; 9 s/?/$pathY/; 10 s/?/$M/; 11 s/?/$X/; 12 s/?/$Y/; 14 s/?/$label/' > $label.R
	R CMD BATCH $label.R
	" > $outdir/$label.cmd
	cd $outdir
	bsub -q short -n 40 -J sum.$label -o $label.LOG -e $label.ERR < $label.cmd
done
done
done

# label=top; echo "X M Y iv X2Y.b X2Y.p TOT.b TOT.p ACME.b ACME.p ADE.b ADE.p Prop.b Prop.p" > $label.txt; cat */*-*.log | cut -d ' ' -f 3-5,7-9,11-18 | awk 'NF==14 && $9*$11 >0 && $10 <0.001 && $14<0.001 && $13>0.1' | sed 's/"//g' | grep -v "^id" | sort -k 13,13gr >> $label.txt
