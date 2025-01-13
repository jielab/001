#!/bin/bash -l

dir0=/work/sph-huangj
dir=$dir0/data/gwas/ppp
outdir=$dir/clean

for dat in `cd $dir/tmp; ls -d */ | sed 's/\///g'`; do
	if [[ -f $outdir/$dat.gz ]]; then echo $outdir/$dat.gz DONE; continue; fi
	echo $dat
	echo "#!/bin/bash
	module load python/anaconda3/2020.7
	echo -e > $dat
	for chr in {1..22} X; do
		fname=\`ls $dir/tmp/$dat/discovery_chr\${chr}_*\`
		python3 $dir0/scripts/main/library/join_file.py -i \"\$fname,SPACE,2 $dir/map.clean/chr\$chr.map,TAB,0\" -o $dat.chr\$chr.ID	
		if [[ \$chr == 1 ]]; then
			cat $dat.chr\$chr.ID | sed 's/ /\t/g' >> $dat
		else
			cat $dat.chr\$chr.ID | sed '1d' | sed 's/ /\t/g' >> $dat
		fi
		rm $dat.chr\$chr.ID
	done
	cat $dat | awk '{if (NR==1) \$13=\"P\"; else if (\$13 !=\"NA\") \$13=10^-(\$13); \$3=\$9=\$12=\$14=\$15=\$16=\$17=\$20=\"\"; print \$0}' | sort -k 1,1n -k 12,12n | sed '1 s/CHROM GENPOS/CHR POS38/; 1 s/ALLELE0 ALLELE1 A1FREQ/NEA EA EAF/; 1 s/POS19/POS/; s/  */\t/g' | gzip -f > $dat.gz
	rm $dat
	" > $outdir/$dat.cmd
	cd $outdir
#	bsub -q short -n 40 -J $dat -o $dat.LOG -e $dat.ERR < $dat.cmd
done


dir=$dir0/data/gwas/ppp
outdir=$dir/cojo; mkdir -p $outdir
for dat in `cd $dir/clean; ls *.gz | sed 's/\.gz//g'`; do
	read chr pos0 pos1 prot<<< $(grep -w "$dat" $dir/ppp.b19.bed) # GRCH 38
	((pos0=pos0-100000)); if ((pos0 < 0)); then pos0=0; fi; ((pos1=pos1+100000))
	echo DO dat $dat, chr $chr, pos0 $pos0, pos1 $pos1, prot $prot

	echo "#!/bin/bash
	zcat $dir/clean/$dat.gz | awk 'NR==1 || ((\$1==\"$chr\" || (\"$chr\"==\"X\" && \$1==23)) && \$12>=$pos0 && \$12 <=$pos1 && \$10<=5e-08)' | sed 's/^23/X/' > $dat.5e-08
        zcat $dir/clean/$dat.gz | awk 'NR==1 || \$10<=5e-8 {x=\$12/1e+05; y=int(x); b=(x>y?y+1:y); print \$11,\$1,\$12,b,\$10}' | \
        sort -k 2,2n -k 4,4n -k 5,5g | awk '{block=\$2\".\"\$4; if (arr[block] !=\"Y\") print \$1; arr[block] =\"Y\"}' | sed 's/^23/X/' > $dat.NEW.top.snp

        zcat $dir/clean/$dat.gz | awk 'NR==1 || ((\$1==\"$chr\" || (\"$chr\"==\"X\" && \$1==23)) && \$12>=$pos0 && \$12 <=$pos1) {print \$11, \$4, \$3, \$5, \$8, \$9, \$10, \$7}' > $dat.4gcta
        plink2 --pgen $dir0/data/ukb/gen/imp/chr$chr.pgen --psam $dir0/data/ukb/gen/imp/chr$chr.psam --pvar $dir0/data/ukb/gen/imp/chr$chr.pvar --keep $dir0/data/ukb/phe/common/ukb.white --chr $chr --from-bp $pos0 --to-bp $pos1 --make-bed --out $dat # --maf 0.01 --hwe 1e-6
	gcta --bfile $dat --cojo-file $dat.4gcta --cojo-slct --cojo-p 5e-6 --out $dat
	" > $outdir/$dat.cmd
        cd $outdir
	bsub -q short -n 40 -J $dat -o $dat.LOG -e $dat.ERR < $dat.cmd
done
