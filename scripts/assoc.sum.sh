#!/bin/bash

dir=/work/sph-huangj
pathX=ppp
pathY=main
cis=1 # 通常设为0，当进行PPP的cis-pQTL分析时设为0

for Y in bald12; do
	outdir=$dir/analysis/assoc.sum/$Y; mkdir -p $outdir

	for X in `cd $dir/data/gwas/$pathX/clean; ls *gz | sed 's/\.gz//g'`; do 
	echo $X
	if [ $X == $Y ]; then continue; fi

	for pathM in NA; do # bb bc inf mb mb2 met ppp
	label="$X-$pathM-$Y"
	if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi

	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	cat $dir/scripts/main/assoc.sum.R | sed '9 s/?/$pathX/; 10 s/?/$pathY/; 11 s/?/$pathM/; 13 s/?/$X/; 14 s/?/$Y/; 15 s/?/$label/g' > $label.R
	" > $outdir/$label.cmd
	
	if [ $cis==1 ]; then 
		read chr pos_begin pos_end < <(awk -v p=$X '$4==p {print $1, $2, $3}' $dir/files/glist-hg19)
		echo "sed -i '12 s/cis=0/cis=1; chr=$chr; pos_begin=$pos_begin; pos_end=$pos_end; flank=1000000/' $label.R" >> $outdir/$label.cmd
	fi
	echo "R CMD BATCH $label.R" >> $outdir/$label.cmd

	cd $outdir
#	bsub -q short -n 40 -J sum.$label -o $label.LOG -e $label.ERR < $label.cmd
done
done
done

# awk 'FNR==1 && $1 !~ /SKIP/' */*-met-*.log | cut -d ' ' -f 2,4,6-8 | sort -k 5,5g > mr.top.txt
# echo "X M Y iv X2Y.b X2Y.p TOT.b TOT.p ACME.b ACME.p ADE.b ADE.p Prop.b Prop.p" > mrMed.top.txt; 
# cat */*.log | cut -d ' ' -f 2-4,6-8,10-17 | awk 'NF==14 && $9*$11 >0 && $8<0.01 && $10<0.01 && $14<0.01 && $13>=0.1' | sed 's/"//g' | grep -v "^id" | sort -k 13,13gr >> mrMed.top.txt

