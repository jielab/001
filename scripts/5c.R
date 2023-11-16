setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, dplyr, tidyverse, reshape2, ggplot2, TwoSampleMR, MendelianRandomizaiton)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1: Correlation 主要指在基因水平上
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read.table('D:/analysis/ldsc/all.rg.res', header=T)
rg <- dat %>% select(p1, p2, rg) %>% acast(p1 ~ p2, value.var='rg'); rg[is.na(rg)] =0;  rg=round(rg,1)
pval <- dat %>% select(p1, p2, p) %>% acast(pval, p1 ~ p2, value.var='p')
plt <- ggcorrplot(rg, lab=T, p.mat=pval, sig.level=5e-4, insig ='blank') 
	plt + theme(axis.title=element_text(size=15, face='bold'), axis.text=element_text(size=12, face='bold'))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2: Causation 主要是通过MR分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 可以用来比较两个GWAS的beta
IV <- read.table("D:/data/gwas/pheno/bb_CRP.top.snps", header=T)
dat_X <- read.table("D:/data/gwas/pheno/bb_CRP.gz", header=T) %>% merge(IV, by="SNP") %>% 
	format_data(type ="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
dat_Y <- extract_outcome_data(dat_X$SNP, "ieu-b-35") # CRP
dat <- harmonise_data(dat_X, dat_Y, action=2)
plot(dat$beta.exposure, dat$beta.outcome)
res <- mr(dat); p1 <- mr_scatter_plot(res, dat); p1[[1]]
bsmr_dat <- TwoSampleMR::dat_to_MRInput(dat); bsmr_dat <- bsmr_dat[[1]]
MendelianRandomization::mr_plot(bsmr_dat)
# TwoSampleMR 标准流程
out_dat <- read_outcome_data(filename="D:/data/gwas/penguin/y.t2d.sub.txt", sep="\t", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
exp_dat0 <- read.table("D:/Downloads/ng2022.txt", sep="\t", header=T)
exp_dat <- format_data(exp_dat0[1,], type ="exposure", phenotype_col="PHENOTYPE", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat)
mr(dat)
run_mr_presso(dat)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C3: Colocalization 从全局到局部local
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 “详细解读！利用ieu数据库进行GWAS-GWAS共定位分析”
pacman::p_load(gwasglue, dplyr, gassocplot, coloc)
top <- ieugwasr::tophits('ieu-a-300') %>% arrange(p) #读取服务器上的数据
	chrpos <- paste0(top$chr[1], ":", top$position[1] - 90000, "-", top$position[1] + 90000)
	out <- ieugwasr_to_coloc(id1='ieu-a-300', id2='ieu-a-7', chrompos=chrpos)
	res <- coloc::coloc.abf(out[[1]], out[[2]])
chrpos <- "19:11112306-11292306" # 下载VCF文件到本地
	out <- gwasvcf_to_coloc("ieu-a-300.vcf.gz", "ieu-a-7.vcf.gz", chrpos)
	res <- coloc::coloc.abf(vout[[1]], vout[[2]])
temp <- coloc_to_gassocplot(out)
gassocplot::stack_assoc_plot(temp$markers, temp$z, temp$corr, traits=temp$traits)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C4: Coevolution mRNA及蛋白质互作
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C5: Community 共同体
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~