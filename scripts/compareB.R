## 如果进行比较的两个文件有一个文件太大（几百万行），会导致R读进去很慢，这个时候先用下面的 Linux 代码把大的那个文件变小
## trait1=XXX; trait2=YYY; zcat $trait2.gwas.txt | fgrep -wf $trait1.top.snps > $trait1.$trait2.tmp

setwd('C:/Users/jiehu/Desktop')
source('D:/scripts/library/compareB.f.R')

pdf('compareB.pdf')
par(mfrow=c(3,2), mai=c(1,1,0.5,0.5))
compareB(
  f1='D:/data/gwas/pheweb/heart_f2020.gwas.txt.gz', f1_name='hf2020', f1_snp='SNP', f1_ea='EA', f1_nea='NEA', f1_eaf=NA, f1_beta='BETA', f1_se='SE', f1_p='P',
  f2='D:/data/gwas/pheweb/heart_failure.gwas.txt.gz', f2_name='hf2018', f2_snp='SNP', f2_ea='EA', f2_nea='NEA', f2_eaf=NA, f2_beta='BETA', f2_se='SE', f2_p='P'
)
dev.off()

