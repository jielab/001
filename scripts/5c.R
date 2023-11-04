setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, dplyr, tidyverse, reshape2, ggplot2, corrplot, ggcorrplot, hyprcoloc, mediation)

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
# TwoSampleMR 标准流程
out_dat <- read_outcome_data(filename="D:/data/gwas/penguin/y.t2d.sub.txt", sep="\t", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
exp_dat0 <- read.table("D:/Downloads/ng2022.txt", sep="\t", header=T)
exp_dat <- format_data(exp_dat0[1,], type ="exposure", phenotype_col="PHENOTYPE", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", eaf_col="EAF", beta_col="BETA", se_col="SE", pval_col="P")
dat <- harmonise_data(exposure_dat=exp_dat, outcome_dat=out_dat)
mr(dat)
run_mr_presso(dat)
  exposure_dat <- mv_extract_exposures(c("ieu-a-299", "ieu-a-300", "ieu-a-302")) #提取多个暴露数据
  outcome_dat <- extract_outcome_data(exposure_dat$SNP, "ieu-a-7") #提取结局数据
  mvdat <- mv_harmonise_data(exposure_dat, outcome_dat)  #合并暴露数据与结局数据
  res <- mv_multiple(mvdat)  #进行多变量MR分析
# 多个分析的结果合并、整理
dat <- read.table('D:/analysis/mr/pheno.res', sep='\t', header=F, as.is=T) %>% subset(V6 %in% c("Wald ratio", "Inverse variance weighted"))
	names(dat) <- c('exp_name', 'exp_cnt', 'out_name', 'out_cnt', 'dat_cnt', 'method', 'beta', 'se', 'p')
	dat$p=signif(dat$p,2)
	# dat %>% subset(exp_name %like% "Bifido" & out_name =="t2d") 
	dat.beta <- dat %>% dplyr::select(exp_name, out_name, beta) %>% acast(exp_name ~ out_name, value.var='beta'); dat.beta[is.na(dat.beta)] =0; dat.beta=round(dat.beta,1)
	dat.p <- dat %>% dplyr::select(exp_name, out_name, p) %>% acast(exp_name ~ out_name, value.var='p'); dat.p[is.na(dat.p)] =1
	write.table(dat.p, file="mr.p.txt", sep='\t', row.names=T, col.names=T, append=F, quote=F)
plt <- ggcorrplot(dat.beta, lab=T, p.mat=dat.p, sig.level=.25e-4, insig ='blank') 
	plt + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
#corrplot(beta, is.corr=F, method='shade', bg='black', col=colorRampPalette(c('white','green','gold'))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig='pch', pch.cex=2, tl.srt=45, outline=T)

	
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