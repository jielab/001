#!/bin/bash

dir=/work/sph-huangj
folderX=ppp
folderY=main
folderM="mb met" # bb bc mb met ppp  # 改变之前的代码，避免X-Y为每一个folderM重复跑🏃‍

cis=1 # 通常设为0，当进行PPP的cis-pQTL分析时设为1
cisMRcML=0
coloc=0
mrMed=1 # 通常设为0，会消耗大量计算🏂

for Y in bald12 bald13 bald14 bald134 bald1234; do
	outdir=$dir/analysis/assoc.sum.NEW/$Y; mkdir -p $outdir

	for X in CA5A EPGN PAM PXN FES UXS1; do # ALL
	echo $X
	# if [ -f $outdir/$label.log ]; then echo $label already run; continue; fi
	if [ $X == $Y ]; then continue; fi
	label="$X-$Y"

	if [ $cis == 1 ]; then 
		str=`awk -v p=$X '{if ($4==p) print $1, $2, $3}' $dir/files/ppp.glist-hg19`
		read -r chr pos_begin pos_end <<< $str
		if [ $chr == "" ]; then
			cis_str=0 # 没办法，找不到位置
		else
			cis_str="1; chr=$chr; pos_begin=$pos_begin; pos_end=$pos_end; flank=1000000"
		fi
	else
		cis_str=0
	fi
	
	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	source activate
	conda activate R
	cat $dir/scripts/main/assoc.sum.R | sed '8 s/?/$folderX/; 9 s/?/$folderY/; 10 s/?/$folderM/; 11 s/?/$X/; 12 s/?/$Y/; 13 s/?/$label/; 15 s/?/$cis_str/; 16 s/?/$cisMRcML/; 17 s/?/$coloc/; 18 s/?/$mrMed/' > $label.R
	R CMD BATCH $label.R
	" > $outdir/$label.cmd

	cd $outdir
	bsub -q short -n 40 -J sum.$label -o $label.LOG -e $label.ERR < $label.cmd
done
done

# awk 'FNR==1 && $1 !~ /SKIP/' */*.mr.log | cut -d ' ' -f 2,4,6-8 | sort -k 5,5g > mr.top.txt
# echo "X M Y iv X2Y.b X2Y.p TOT.b TOT.p ACME.b ACME.p ADE.b ADE.p Prop.b Prop.p" > mrMed.top.txt; 
# cat */*mrMed.log | cut -d ' ' -f 2-4,6-8,10-17 | awk 'NF==14 && $9*$11 >0 && $8<0.01 && $10<0.01 && $14<0.01 && $13>=0.1' | sed 's/"//g' | grep -v "^id" | sort -k 13,13gr >> mrMed.top.txt

