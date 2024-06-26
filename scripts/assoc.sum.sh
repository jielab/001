#!/bin/bash

dir=/work/sph-huangj/analysis

for Y in y.cad y.vte covid_A2 covid_B2 covid_C2; do # t2dm.EUR y.adhd y.asd y.cad y.dementia y.dep y.ibs y.mi y.ptsd y.scz y.stroke y.t2dm y.varicose y.vte; do
	outdir=$dir/assoc.sum/$Y; mkdir -p $outdir
	for X in ABO height walk_pace x.age_fbirth x.age_fsex x.bmi; do # x.cannabis x.crp x.drink x.education x.hdl x.ldl x.leisure x.lst x.mvpa x.smoke x.vitad lipid_hdl.EUR lipid_ldl.EUR lipid_tc.EUR lipid_tg.EUR; do
	for Z in bb bc met mb ppp; do
	label="$X-$Z-$Y"
	if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi

	echo "
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	cat /work/sph-huangj/scripts/main/assoc.sum.R | sed -e '11 s/?/$X/' -e '12 s/?/$Y/' -e '10 s/?/$Z/' -e '14 s/?/$label/' > $label.R
	R CMD BATCH $label.R
	" > $outdir/$label.cmd
	cd $outdir
	bsub -q short -n 40 -J $label -o $label.LOG -e $label.ERR < $label.cmd
done
done
done

# label=walk; cat */$label*-*.log | awk '$24>0 && $27<0.05 {$1=""; print $0}' | sed 's/^ //; s/"//g; s/| /|\t/g' | sort -k 23,23gr > $label.top.txt

