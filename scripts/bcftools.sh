# 技术讨论 https://github.com/freeseek/score/issues/26

#!/bin/bash

dir0=/work/sph-huangj
trait=height

for dat in AFR EUR HIS SAS EAS; do
	bcftools +munge -Ou -f $dir0/data/gatk-bundle/hg19/ucsc.hg19.fasta -C $dir0/files/colheaders.tsv -s $dat $dir0/data/gwas/main/clean/$trait.$dat.gz | \
	bcftools +liftover -Ou -- -s $dir0/data/gatk-bundle/hg19/ucsc.hg19.fasta -f $dir0/data/gatk-bundle/hg38/Homo_sapiens_assembly38.fasta -c $dir0/files/hg19ToHg38.over.chain.gz | \
	bcftools sort -Ob -o $trait.$dat.bcf --write-index
#	bcftools view $trait.$dat.bcf -O v -o $trait.$dat.vcf
done

bcftools merge $trait.{AFR,EUR,HIS,SAS,EAS}.bcf -m none -o $trait.bcf -Ob --write-index
bcftools +pgs -Ov -o $trait.pgs.vcf --log $trait.log --beta-cov 5e-08 \
	--samples AFR,EUR,HIS,SAS,EAS $trait.bcf $dir0/data/ldref/1kg_ldgm.{AFR,EUR,HIS,SAS,EAS}.bcf

