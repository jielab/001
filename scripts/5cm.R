pacman::p_load(data.table, dplyr, tidyverse, reshape2, ggplot2, corrplot, ggcorrplot, hyprcoloc, mediation)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C1: genetic correlation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read.table('D:/analysis/ldsc/all.rg.res', header=T)
rg <- dat %>% select(p1, p2, rg) %>% acast(p1 ~ p2, value.var='rg'); rg[is.na(rg)] =0;  rg=round(rg,1)
pval <- dat %>% select(p1, p2, p) %>% acast(pval, p1 ~ p2, value.var='p')
plt <- ggcorrplot(rg, lab=T, p.mat=pval, sig.level=5e-4, insig ='blank') 
	plt + theme(axis.title=element_text(size=15, face='bold'), axis.text=element_text(size=12, face='bold'))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C2: Mendelian Randomization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 个人数据与“无人”数据
set.seed(12345)
dat <- iris %>% rename(varX=Sepal.Length) 
dat <- dat %>% mutate(
	ran = runif(nrow(dat), min=min(varX), max=max(varX)),
	varY = varX*0.35 + ran*0.65,
	snp = ifelse(varX > quantile(varX, probs=0.90), 2, ifelse(varX < quantile(varX, probs=0.10), 0, 1))
	)
dat$varX.pred = predict.lm( lm( varX ~ snp, data=dat)) 
summary(lm(varY ~ varX.pred, data=dat))
beta.snp2varX <- summary(lm( varX ~ snp, data=dat))$coef[2,1]
beta.snp2varY <- summary(lm( varY ~ snp, data=dat))$coef[2,1]
se.snp2varX <- summary(lm( varX ~ snp, data=dat))$coef[2,2]
se.snp2varY <- summary(lm( varY ~ snp, data=dat))$coef[2,2]
beta.wald <- beta.snp2varY / beta.snp2varX; beta.wald
se.wald <- se.snp2varY / beta.snp2varX
p.wald <- 2*pnorm(abs(beta/se), lower.tail=F)
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
# 多多个分析的结果合并、整理
dat <- read.table('D:/analysis/psy.mr.res', sep='\t', header=F, as.is=T) %>% subset(V6 %in% c("Wald ratio", "Inverse variance weighted"))
	names(dat) <- c('exp_name', 'exp_cnt', 'out_name', 'out_cnt', 'dat_cnt', 'method', 'beta', 'se', 'p')
#	dat %>% subset(exp_name %like% "Bifido" & out_name =="t2d") 
dat.beta <- dat %>% dplyr::select(exp_name, out_name, beta) %>% acast(exp_name ~ out_name, value.var='beta'); dat.beta[is.na(dat.beta)] =0; dat.beta=round(dat.beta,1)
dat.p <- dat %>% dplyr::select(exp_name, out_name, p) %>% acast(exp_name ~ out_name, value.var='p'); dat.p[is.na(dat.p)] =1
plt <- ggcorrplot(dat.beta, lab=T, p.mat=dat.p, sig.level=.25e-4, insig ='blank') 
	plt + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
#corrplot(beta, is.corr=F, method='shade', bg='black', col=colorRampPalette(c('white','green','gold'))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig='pch', pch.cex=2, tl.srt=45, outline=T)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# C3: Colocalization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- read.table("D:/analysis/coloc/coloc.txt.gz", header=T, as.is=T) %>% rename(locus=ChrPos.m); names(dat)
sink("output.txt")
for (loc in unique(dat$locus)) {
	print(loc)
	dat1 <- dat %>% subset(locus==loc) %>% na.omit()
	betas <- subset(dat1, select=grepl("BETA", names(dat))); betas <- as.matrix(betas)
	ses <- subset(dat1, select=grepl("^SE", names(dat))); ses <- as.matrix(ses)
	traits <- gsub("BETA.", "", grep("BETA", names(dat1), value=T))
	rsid <- as.matrix(dat1[,1])
	res <- hyprcoloc(betas, ses, trait.names=traits, trait.subset=c("x", "y.dementia", "y.depress"), snp.id=rsid)
	print(res)
}
sink()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(data.table, tidyverse, mediation)
set.seed(12345)
dat <- iris %>% rename(varX=Sepal.Length) 
dat <- dat %>% 
	mutate(
	ran = runif(nrow(dat), min=min(varX), max=max(varX)),
	varY = varX*0.35 + ran*0.65,
	snp = ifelse(varX > quantile(varX, probs=0.90), 2, ifelse(varX < quantile(varX, probs=0.10), 0, 1))
	)
dat$varX.pred = predict.lm( lm( varX ~ snp, data=dat))
summary(lm(varY ~ varX.pred, data=dat))
beta.snp2varX <- summary(lm( varX ~ snp, data=dat))$coef[2,1]; beta.snp2varX
beta.snp2varY <- summary(lm( varY ~ snp, data=dat))$coef[2,1]; beta.snp2varY
se.snp2varX <- summary(lm( varX ~ snp, data=dat))$coef[2,2]; se.snp2varX
se.snp2varY <- summary(lm( varY ~ snp, data=dat))$coef[2,2]; se.snp2varY
beta.wald <- beta.snp2varY / beta.snp2varX; beta.wald
se.wald <- se.snp2varY / beta.snp2varX; se.wald
2*pnorm(abs(beta.wald / se.wald), lower.tail=FALSE); 2*pt(-abs(beta.wald / se.wald),df=150-1)
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
