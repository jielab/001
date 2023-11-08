#!/bin/bash -l

## sort -k 1 -V SFARI_cnv.b37.bed | bgzip -c > SFARI_cnv.b37.bed.gz; tabix -fpbed SFARI_CNV.b37.bed.gz
## bcftools view -h abc.vcf.gz | awk '$1 !~/^##contig/' > new.header.txt; bcftools reheader -h new.header.txt
## bcftools annotate freeze.7a.chr$chr.pass_and_fail.sites.vcf.gz --include 'INFO/AC >1' --remove ^INFO/AVGDP,INFO/AC,INFO/AF
## vep --vcf adds "--symbol" by default, not working with "--most_severe", --pick --pick_order tsl,appris,rank

## ClinVAR ##
cat raw.txt | sed 's/;not provided//; s/;not specified//; s/\t\t/\tNA\t/g; s/\t\t/\tNA\t/g' | awk 'BEGIN{FS="\t"; OFS="\t"}{$1="rs"$1; split($NF,a,";"); for (var in a) {$NF=a[var]; print $0 }}'
bcftools query -f '%CLNSIG %CLNREVSTAT\n' clinvar.b38.vcf.gz | sort | uniq -c
bcftools filter --include 'CLNSIG ~"athogenic" || CSQ[*] ~ "missense_variant.*deleterious"' # CLNREVSTAT=="reviewed_by_expert_panel" || CLNREVSTAT=="practice_guideline"
seq 1 22 | xargs -n1 -I % echo % chr% > chr_name_conv.txt; echo -e "X chrX\nY chrY\nMT chrMT" >> chr_name_conv.txt
bcftools annotate clinvar.vcf.gz --rename-chrs chr_name_conv.txt -Oz -o clinvar.b38.vcf.gz


join_file=/mnt/d/scripts/library/join_file.py
dir=/mnt/d/data/sequencing
adir=/mnt/d/files/annotation
project=N3305
ver=38
samples=`awk '{printf " "$1}' $dir/$project/sample.ped` # proband, father, mother
samples2=`echo $samples | sed 's/ /,/g'`


## merge files ##
for sample in $samples; do
	for type in snp indel; do
        bcftools annotate $dir/$project/$type/Vcf/$sample.GATK.$type.vcf.gz --include 'FORMAT/DP[*] >10' --remove ^INFO/DP -Oz -o $sample.$type.vcf.gz
        tabix $sample.$type.vcf.gz
    done
    bcftools concat $sample.snp.vcf.gz $sample.indel.vcf.gz -a -Oz -o $sample.concat.vcf.gz 
    tabix $sample.concat.vcf.gz
done
bcftools merge *concat.vcf.gz -m both -Oz -o all.tmp1.vcf.gz
tabix all.tmp1.vcf.gz


## merge AF and Annotations ##
bcftools annotate all.tmp1.vcf.gz -a $adir/Topmed.vcf.gz -c Topmed_AF:=AF --mark-sites In_Topmed -Oz -o all.tmp2.vcf.gz; tabix all.tmp2.vcf.gz
bcftools query -l all.tmp2.vcf.gz; bcftools query -f '%CHROM\n' all.tmp2.vcf.gz | sort | uniq -c
bcftools filter all.tmp2.vcf.gz --exclude '(Topmed_AF >=0.01 && Topmed_AF <=0.99)' -Oz -o all.tmp3.vcf.gz
echo '##INFO=<ID=SFARI_gene,Number=1,Type=String,Description="disease region">' > SFARI.gene.h
echo '##INFO=<ID=SFARI_cnv,Number=1,Type=String,Description="disease region">' > SFARI.cnv.h
bcftools annotate all.tmp3.vcf.gz -a $adir/SFARI_genes.b$ver.bed.gz -c CHROM,FROM,TO,SFARI_gene -h SFARI.gene.h --mark-sites In_SFARI_genes -Oz -o all.tmp4.vcf.gz; tabix all.tmp4.vcf.gz
bcftools annotate all.tmp4.vcf.gz -a $adir/SFARI_cnv.b$ver.bed.gz -c CHROM,FROM,TO,SFARI_cnv -h SFARI.cnv.h --mark-sites In_SFARI_cnvs -Oz -o all.tmp5.vcf.gz; tabix all.tmp5.vcf.gz
bcftools annotate all.tmp5.vcf.gz -a $adir/clinvar.b$ver.vcf.gz -c CLNDISDB,CLNDN,CLNREVSTAT,CLNSIG,MC,RS --mark-sites In_CLINVAR -Ov -o all.tmp6.vcf


## add VEP prediction ##
vep -i all.tmp6.vcf --vcf --cache --dir /mnt/d/data/vep --offline --assembly GRCh$ver -o all.tmp7.vcf --force_overwrite --sift b --canonical --symbol --pick --af_1kg --freq_pop 1KG_EUR --af_gnomad --pubmed --failed 1
bgzip -c all.tmp7.vcf > all.tmp7.vcf.gz
bcftools +split-vep all.tmp7.vcf.gz -l
bcftools +split-vep all.tmp7.vcf.gz -c Consequence,IMPACT,SYMBOL,CANONICAL,SIFT -s worst | bcftools view -Oz -o all.final.vcf.gz; tabix all.final.vcf.gz
bcftools query all.final.vcf.gz -f '%IMPACT %Consequence\n' | sort | uniq -c
# !! view VEP summary !!


## For each proband, find top SNV ##
bcftools filter all.final.vcf.gz --include '(In_CLINVAR=1 && CLNSIG !~ "onflict" && CLNSIG ~ "athogenic") || (IMPACT="HIGH")' -Oz -o all.top.vcf.gz # && (In_SFARI_genes=1 || In_SFARI_cnvs=1)
bcftools view all.top.vcf.gz --samples $samples2 --no-update | \
bcftools filter --include 'GT[0] !~"\."' | \
bcftools query -H -f '%ID\t%CHROM:%POS:%REF:%ALT\t%DP\t%Topmed_AF\t%SFARI_gene\t%SFARI_cnv\t%CLNDISDB\t%CLNDN\t%CLNREVSTAT\t%CLNSIG\t%MC\t%Consequence\t%IMPACT\t%SYMBOL\t%CANONICAL\t%SIFT[\t%SAMPLE=%GT]\t%CSQ\n' > final.txt


## For each proband, find top CNV ##
#vep --plugin StructuralVariantOverlap,file=gnomad_v2_sv.sites.vcf.gz
for sample in $samples; do
#	cp $dir/$project/cnv/cnvnator/$sample.cnv.bed ./
#	awk 'BEGIN{OFS="\t"}{if ($1 !~/^@/ && $1 ~/^chr/) print $1,$2,$3,"GATK"}' $dir/$project/cnv/modelFinal/$sample.modelFinal.seg > $sample.cnv.bed
#	bcftools query $dir/$project/gcnv/$sample/genotyped-segments-*-$sample.vcf.gz -f '%CHROM %POS %INFO/END [%CN]\n' | awk '{if ($NF<=1) $NF="deletion"; else if ($NF>2) $NF="duplication"; if ($NF !=2) print $0}' | sed 's/ /\t/g' > $sample.cnv.bed
	awk -F "\t" 'NR>1 && $24 !="breakpoint" {print "chr"$1,$2,$3,$21}' $dir/$project/CNV/Annotation/$sample.final.bam_CNVs.final.hg38_multianno.xls | sed 's/ /\t/g' > $sample.cnv.bed
	bedtools intersect -a $sample.cnv.bed -b $adir/SFARI_cnv.b$ver.bed -wo > $sample.cnv.tmp1
	awk '{len1=$3-$2; len2=$7-$6; lap=$9; r1=len1/lap; r2=len2/lap; print FILENAME"\t"$0"\t|", len1, len2, r1,r2}' $sample.cnv.tmp1 | sed 's/.cnv.tmp1//' > $sample.cnv.tmp2
done
cat *.cnv.tmp2 | awk '{print $9,$1}' | sort | uniq | awk '{array[$1]=array[$1]","$2} END{for(key in array) print key,array[key]}' > vip.cnv.tmp
python $join_file -i "vip.cnv.tmp,SPACE,0 $adir/SFARI_cnv.b$ver.bed,TAB,3" -o vip.cnv
sort -k 7nr vip.cnv > vip.cnv.txt


### check inheritance patterns (de novo vs. heritable) ##
module load genmod/3.7.2
module load gcc/7.2.0
module load denovogear/2018-05-1_6723027
dng dnm auto --ped ../vcf/sample.ped --bcf step1.final.bcf > step2.dng.out
awk 'BEGIN{print "#CHROM POS REF ALT DNG"}{if ($5 >=1 && $5<=22) print $5,$7,$9,$11,"Y"}' step2.dng.out | \
	sort -k 1,1n -k 2,2n | sed 's/ /\t/g' | bgzip > step2.dng.txt.gz
tabix -f -s 1 -b 2 -e 2 -S 1 step2.dng.txt.gz
echo '##INFO=<ID=DNG,Number=1,Type=String,Description="denovogear analysis results">' > dng.h
bcftools annotate -a step2.dng.txt.gz -c 'CHROM,POS,REF,ALT,DNG' -h dng.h -Oz -o step2a.vcf.gz step1.final.bcf
# $dir/software/from_SCC/triodenovo.0.05/bin/triodenovo --ped ../vcf/sample.ped --in_vcf step2a.vcf.gz --out_vcf step2b.tmp.vcf
genmod annotate step2a.vcf.gz -r ?? --annotate_regions | genmod models - --family_file ../vcf/sample.ped | bcftools view -Ob -o step2b.bcf
tabix step2b.bcf
bcftools filter --include 'DNG ="Y" || GeneticModels ~ "_dn"' -Ob -o step2.final.bcf step2b.bcf
bcftools query -f '%CSQ\n' step2.final.bcf | awk -F "|" '{print $2}' | sort | uniq -c 

