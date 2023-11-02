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
# C4: Community 共同体
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#fit.totaleffect <- lm(varY ~ varX, dat); summary(fit.totaleffect) # coef=0.12984, about 35% *35%.
fit.varM <- lm(varM ~ varX, dat); summary(fit.varM) # beta=0.30429
fit.full <- lm(varY ~ varX+varM, dat); summary(fit.full) # varM beta=0.37194; varX beta=0.01667, but no more significant, i.e., "complete mediation"
res <- mediate(fit.varM, fit.full, treat='varX', mediator='varM', boot=T); summary(res) #ACME=0.1132 (i.e., 0.30429 * 0.37194); ADE=0.01667
fit.snp <- lm(varX ~ snp, dat); summary(fit.snp) 
# ukb data
dat <- dat0 %>% filter(ethnicity_gen==1) %>% dplyr::select(grep("age|sex|bmi|covid|^bb_|^bc|age_mn",names(dat0),value=T)) %>% 
	mutate (
	age_mn_3p = ifelse(age_mn %in% 12:14, "normal", ifelse(age_mn>14, "late", "early")),
	age_mn_3p = factor(age_mn_3p, levels=c("normal", "early", "late")),
	outcome = ifelse(!is.na(icdDate_covid), 1, ifelse(covid_inf==1,0, NA))
	) %>% rename(CYS=bb_CYS)
varX="age_mn_3p"; varY="outcome"
form <- formula(paste(varY, "~", varX, "+age+sex+bmi+CYS"))
form <- formula(paste(varY, "~", varX, "+age+sex+bmi+CYS+", paste(grep("^bb_", names(dat), value=T), sep="", collapse="+")))
summary(glm(form, data=dat))
for (varM in grep("CYS|^bb_|^bc_", names(dat), value=T)) {
	print(paste("M变量:", varM))
	dat1 <- subset(dat, select=c(varX, varY, varM, "age", "sex", "bmi")) %>% na.omit()
	dat1[[varM]] <- inormal(dat1[[varM]]) # normal transformation
	names(dat1) <- c("varX", "varY", "varM", "age","sex", "bmi")
	dat1$varX <- as.factor(dat1$varX)
	fit.med = lm(varM ~ varX + age+sex, data=dat1)
	fit.out = glm(varY ~ varX + varM +age+sex, data=dat1, family=binomial("probit"))
	med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", sims=100, boot=T)
	print(summary(med.out))
}
fit.med = lm(varM ~ varX, data=dat1)
fit.out = glm(varY ~ varX + varM, data=dat1, family=binomial("probit"))
med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", sims=100, boot=T)
