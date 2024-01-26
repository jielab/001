#!/bin/bash -l

dir=/work/sph-huangj
outdir=/scratch/2023-05-17/sph-huangj
mkdir -p $outdir
echo -e "66137\na93100fe54833270898b8e87f96646fdf52f509de2b87248ed1ea92b560c2ba1" > .ukbkey


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download and format
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
wget -nd biobank.ndph.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_bim.tar
for chr in {1..22} X XY; do # Y MT for typed only
    echo "#!/bin/bash
    echo \"Starting on : \$(date); Running on : \$(hostname); Job ID : \$JOB_ID\"
	gfetch 22418 -c$chr # genotype bed file
	gfetch 22418 -c$chr -m
	gfetch 22438 -c$chr # haplotype
	gfetch 22438 -c$chr -m
	gfetch 22828 -c$chr # imputed genotype
	gfetch 22828 -c$chr -m
	gfetch 23158 -c$chr # WES PLINK bed
	gfetch 23158 -c$chr -m
	" > $outdir/chr$chr.cmd
	cd $outdir # short, medium, large:wq
    bsub -q smp -J ukb.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
done
awk 'FNR!=1 && $1 <0 {print $1,$2}' chr*.psam | sort | uniq > negative.sam.id # 人为加上 5456019
for chr in {1..22} X XY; do
    sam_typ=`ls -1 $dir/ukb/ukb22418_c${chr}_b0_v2_s*.fam | awk '{printf $1}'`
    sam_hap=`ls -1 $dir/ukb/ukb22438_c${chr}_b0_v2_s*.sample | awk '{printf $1}'`
    sam_imp=`ls -1 $dir/ukb/ukb22828_c${chr}_b0_v3_s*.sample | awk '{printf $1}'`
    echo "#!/bin/bash
    echo \"Starting on : \$(date); Running on : \$(hostname); Job ID : \$JOB_ID\"
    cat chr$chr.pvar | awk '{if (arry[\$3]==\"Y\") {\$3=\$3\".DUP\"}; print \$0; arry[\$3]=\"Y\"}' | sed 's/ /\t/g' > chr$chr.NEW.pvar
    plink2 --bed $dir/ukb/ukb22418_c${chr}_b0_v2.bed --bim $dir/ukb/ukb_snp_chr${chr}_v2.bim --fam $sam_typ --remove negative.sam.id --make-bed --out chr$chr
    plink2 --bgen $dir/ukb/ukb22438_c${chr}_b0_v2.bgen ref-first --sample $sam_hap --remove negative.sam.id --oxford-single-chr $chr --make-pgen --out chr$chr.hap
    plink2 --bgen $dir/ukb/ukb22828_c${chr}_b0_v3.bgen ref-first --sample $sam_imp --remove negative.sam.id --make-pgen --out chr$chr
    " > $outdir/chr$chr.cmd
	cd $outdir # short, medium, large
    bsub -q smp -J ukb.chr$chr -o chr$chr.LOG -e chr$chr.ERR < chr$chr.cmd
done
exit
seq 2 22 | xargs -n1 -I % echo chr%.bed chr%.bim chr%.fam > merge-list.txt
plink --bfile chr1 --merge-list merge-list.txt --make-bed --out ukb
plink --bfile ukb --keep /restricted/projectnb/ukbiobank/jiehuang/files/ukb.keep --make-bed --out ukb.EUR


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Extract genotype of VIP snps
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#!/bin/bash -l
dir=/data/sph-huangj
label=vip
snps=/work/sph-huangj/data/ukb/phe/common/ukb.$label.snp
outdir=/work/sph-huangj/tmp/$label
if [ -d $outdir ]; then rm -r $outdir; fi; mkdir -p $outdir
echo -e "#!/bin/bash -l
for chr in {1..22} X; do
	plink2 --pfile $dir/ukb/imp/chr\$chr --extract $snps --make-pgen --out chr\$chr
done
ls -1 chr*.pgen | awk '{print \$1}' | sort -k 1,1 -V | sed 's/\.pgen//' > merge-list.txt
plink2 --pmerge-list merge-list.txt -delete-pmerge-result --update-name $snps 1 2 --make-pgen --out $label
	plink2 --pfile $label --maj-ref --export A --out $label
	cut -f 2,7- $label.raw | sed 's/\t/ /g' | gzip -f > $label.raw.gz
" > $outdir/$label.cmd
cd $outdir
bsub -q smp -J $label.gen -o $label.LOG -e $label.ERR < $label.cmd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract APOE haplotype
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## e1: C-T;  e2: T-T;  e3: T-C;  e4: C-C
echo -e "rs429358\nrs7412" > apoe.snps
plink2 --pfile $dir/data/ukb/imp/hap/chr19 --extract apoe.snps --export vcf id-paste=iid bgz --out apoe
	tabix apoe.vcf.gz
	bcftools query -l apoe.vcf.gz | awk 'BEGIN{print "IID"}; {print $1}' > sample.ids
	bcftools query -f '%CHROM:%POS:%REF:%ALT[ %GT]\n' -r 19:45411941 apoe.vcf.gz | tr ' ' '\n' > snp1.txt
	bcftools query -f '%CHROM:%POS:%REF:%ALT[ %GT]\n' -r 19:45412079 apoe.vcf.gz | tr ' ' '\n' > snp2.txt
paste sample.ids snp1.txt snp2.txt | sed 's/|/ /g' > apoe.hap.tmp
	awk 'NR>1 {print $1, $2$4, $3$5}' apoe.hap.tmp > apoe.hap.tmp2
	sed -e 's/ 00/ e3/g' -e 's/ 10/ e4/g' -e 's/ 11/ e1/g' -e 's/ 01/ e2/g' apoe.hap.tmp2 > apoe.hap.tmp3
	awk 'BEGIN{print "IID apoe"}{print $1,$2$3}' apoe.hap.tmp3 | sed 's/e2e1/e1e2/;s/e3e1/e1e3/;s/e3e2/e2e3/;s/e4e1/e1e4/;s/e4e2/e2e4/;s/e4e3/e3e4/' > apoe.hap
	awk '{print $2,$3,$4}' apoe.hap | sort | uniq -c


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Download WES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir=/restricted/projectnb/ukbiobank/jiehuang
pwd=`pwd`
awk '{print $1,"23164_0_0"}' $dir/data/ukb/wes/fe.fam | split -d -a 3 -l 100 - list
cnt=0
for dat in `ls list*`; do
	let "cnt=$cnt+1"
	let "gp=($cnt-1)/9"
	let "gp2=$gp +1"
	if [[ $gp2 == 1 ]]; then
		qsub_str="-N g$gp2.$cnt"
	else
		qsub_str="-N g$gp2.$cnt -hold_jid g$gp.*"
	fi
	raw=${dat/list/raw}
	outdir=$dir/data/ukb/wes/BAM/$raw; mkdir -p $outdir
	mv $pwd/$dat $outdir
	echo "#!/bin/bash -l
	echo -n > $dat.res.txt
	module load samtools
	sed 's/23164/23163/' $dat > $dat.2
	ukbfetch -a$dir/files/ukb.key -b$dat
	ukbfetch -a$dir/files/ukb.key -b$dat.2
	for d in \`awk '{printf \" \"\$1}' $dat\`; do
		mv \${d}_23163_0_0.cram \$d.cram
		mv \${d}_23164_0_0.cram.crai \$d.cram.crai
	#	samtools view -L $dir/files/apoe.b38.bed -O BAM -o \$d.bam \$d.cram; samtools index \$d.bam; rm \$d.cram*
	#	samtools view -h \$d.bam | awk -v b=44908684 -v e=44908822 '{if (\$1~/^@/) print \$0; else if (\$1 in readname && pos[\$1]<=b && (\$4+76-1)>=e) print readname[\$1]\"\n\"\$0; readname[\$1]=\$0; pos[\$1]=\$4}' > \$d.tmp.sam
	#	samtools sort \$d.tmp.sam -O BAM -o \$d.pairs.bam; samtools index \$d.pairs.bam; samtools view \$d.pairs.bam > \$d.pairs.sam
	#	cat \$d.pairs.sam | awk '\$6==\"76M\" {if (\$1 in readname) print readname[\$1]\"\t|\t\"\$0; readname[\$1]=\$0}' | awk -v b=44908684 -v e=44908822 '{pos1=b-\$4+1; pos2=e-\$23+1; split(\$10,seq1,\"\"); split(\$29,seq2,\"\"); print seq1[pos1] \"-\" seq2[pos2]}' | sort | uniq -c | awk -v dat=$dat -v d=\$d '{print dat, d, \$2, \$1}' >> $dat.res.txt
	done
	" > $outdir/$dat.cmd
	cd $outdir
	qsub -P ukbiobank $qsub_str -o $dat.LOG -e $dat.ERR < $dat.cmd
done
### re-run failed jobs
ls -l raw*/*ERR | awk '$5 !=0 {print $NF}'| sed 's/raw[0-9]*\///; s/\.ERR//' > failed.list
dir=/restricted/projectnb/ukbiobank/jiehuang
cnt=0
for dat in `cat failed.list`; do
	let "cnt=$cnt+1"
	let "gp=($cnt-1)/9"
	let "gp2=$gp +1"
	if [[ $gp2 == 1 ]]; then
		qsub_str="-N g$gp2.$cnt"
	else
		qsub_str="-N g$gp2.$cnt -hold_jid g$gp.*"
	fi
	raw=${dat/list/raw}
	outdir=$dir/data/ukb/wes/BAM/$raw;
	echo process $cnt, $gp, $gp2, $dat.cmd
	cd $outdir
	qsub -P ukbiobank $qsub_str -o $dat.LOG -e $dat.ERR < $dat.cmd
done
