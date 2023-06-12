setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, lubridate, tidyverse, dplyr, ggplot2, CMplot, TwoSampleMR, MendelianRandomization, survival, survminer, mediation)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read and check data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
	naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
	prop.table(table(dat0$abo, dat0$fut2.rs601338_A),1)
	table(dat0$sp1)
	hist(dat0$lac.rs4988235_A)
	hist(dat0$age_mn)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 万有引力之生存分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) %>% 
	mutate (
		sp1.M = ifelse(sp1=="MM", 2, ifelse(grepl("M", sp1), 1, 0)),
		sp1.S = ifelse(grepl("SS", sp1), 2,  ifelse(grepl("S", sp1), 1, 0)),
		sp1.Z = ifelse(grepl("ZZ", sp1), 2,  ifelse(grepl("Z", sp1), 1, 0))
	)
varXs <- c("o_type", "fut2.rs601338_A", "lac.rs4988235_A", "sp1.Z")
for (varY in c("icdDate_copd", "icdDate_covid", "fod_t1dm", "fod_t2dm", "fod_dm")) {
	print(paste("Y变量:", varY))
	dat1 <- dat %>% 
	mutate(
		outcome_date = dat[[varY]],
		outcome_yes = ifelse( is.na(dat[[varY]]), 0,1),
		follow_end_day = data.table::fifelse(!is.na(outcome_date), outcome_date, fifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 )
	print(table(dat1$outcome_yes))
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$outcome_yes)
	for (varX in varXs) {
		print(paste("X变量:", varX))
		dat1$varX <- dat1[[varX]]
		fit.suv <- survival::survfit(surv.obj ~ varX + sex, data=dat1)
			#print(ggsurvplot(fit.suv, ylim=c(0.95,1), data=dat1)); dev.off() 
		fit.cox <- coxph(surv.obj ~ varX + age+sex, data=dat1)
			#fit.cox <- coxph(surv.obj ~ varX + varX * fut2.rs601338_A + age+sex, data=dat1)
			print(summary(fit.cox))
			png(file=paste(varX,varY,"frt.png",sep="."), w=600, h=800)
			print(ggforest(fit.cox, data=dat1)); dev.off()
	}
}
for (varY in grep("^bb_|^bc_", names(dat), value=T)) {
	print(noquote(varY))
	dat[[varY]] = inormal(dat[[varY]]) # normal transformation
	dat$varY = dat[[varY]] # assign temporary variable name "varY"
	for (varX in varXs) {
		dat$varX = dat[[varX]]
		fit.lm <- lm(varY ~ varX +age+sex, data=dat)
		print(noquote(paste("  --", varX, signif(coef(summary(fit.lm))[2,4],3))))
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mediation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) %>% 
	mutate (
	age_mn_3p = ifelse(age_mn %in% 12:14, "normal", ifelse(age_mn>14, "late", "early")),
	age_mn_3p = factor(age_mn_3p, levels=c("normal", "early", "late")),
	outcome = ifelse(!is.na(icdDate_covid), 1, ifelse(covid_inf==1,0, NA))
	)
varX="age_mn_3p"; varY="outcome"
for (varM in grep("^bb_|^bc_", names(dat), value=T)) {
	print(paste("M变量:", varM))
	dat1 <- subset(dat, select=c(varX, varY, varM, "age", "sex", "bmi")) %>% na.omit()
	dat1[[varM]] <- inormal(dat1[[varM]]) # normal transformation
	names(dat1) <- c("varX", "varY", "varM", "age","sex", "bmi")
	dat1$varX <- as.factor(dat1$varX)
	fit.med = lm(varM ~ varX, data=dat1)
	fit.out = glm(varY ~ varX + varM, data=dat1, family=binomial("probit"))
	med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", sims=100, boot=T)
	print(summary(med.out))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mendelian Randomization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exp_dat0 <- read.table("D:/files/posCtrls/bmi.EUR77.top.txt", header=T, as.is=T)
	exp_iv <- subset(exp_dat0, select="SNP")
	exp_dat <- format_data(exp_dat0, type ="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P")
exp_dat0 <- read.table("D:/Downloads/GCST90019494_buildGRCh37.tsv.gz", header=T, as.is=T) 
	exp_iv <- read.table("D:/data/gwas/bb/bb_ALP.top.snps", header=T, as.is=T)
	exp_dat <- merge(exp_dat0, exp_iv, by.x="variant_id", by.y="SNP")
	exp_dat <- format_data(exp_dat, type ="exposure", snp_col="variant_id", effect_allele_col="effect_allele", other_allele_col="other_allele", beta_col="beta", se_col="standard_error", pval_col="p_value") 
out_dat0 <- read.table("D:/data/gwas/bb/cad.2017.gwas.gz", header=T, as.is=T)
	out_dat <- merge(out_dat0, exp_iv, by.x="SNP", by.y="SNP")
	out_dat <- format_data(out_dat, type ="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
	out_dat <- format_data(out_dat, type ="outcome", snp_col="rsid", effect_allele_col="Allele1", other_allele_col="Allele2", beta_col="Effect", se_col="StdErr", pval_col="P.value") 
	out_dat <- format_data(out_dat, type ="outcome", snp_col="rsid", effect_allele_col="ALT", other_allele_col="REF", beta_col="all_inv_var_meta_beta", se_col="all_inv_var_meta_sebeta", pval_col="all_inv_var_meta_p") 
tsmr_dat <- harmonise_data(exp_dat, out_dat)
	mr(tsmr_dat)
bsmr_dat0 <- dat_to_MRInput(tsmr_dat)
	bsmr_dat <- bsmr_dat0[[1]]
	mr_ivw(bsmr_dat)
	plt <- mr_plot(bsmr_dat, interactive=F)
		plt + theme(axis.title=element_text(size=15, face="bold"), axis.text=element_text(size=12, face="bold"))
	plt <- mr_forest(bsmr_dat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) 
		plt + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))
	plt <- mr_funnel(bsmr_dat)
		plt + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=14), axis.text=element_text(size=12, face="bold"))	
	ggsave(plt, file=paste0(exp,".exp.png"), w=10, h=10)
#Sum Many MR results: 
#	pdftk */*.cB.pdf cat output 001.cB.pdf; cat */*.mr.res.txt | awk 'NF ==10' > 001.mr.txt
pacman::p_load(corrplot, reshape2)
dat <- read.table('001.mr.txt', header=F, as.is=T)
	names(dat) <- c('ex_name', 'ou_name', 'cnt', 'p.ivw', 'p.median', 'p.maxlik', 'p.mbe', 'p.conmix', 'pleio.egger', 'p.egger')
dat1 <- dat %>% group_by(ex_name, ou_name) %>% summarize(Pmin=min(p)) %>% mutate(p= -log10(Pmin))
	dat1 <- dat %>% select(ex_name, ou_name, p.ivw) %>% mutate (p=-log10(p.ivw))
	pmat <- acast(dat1, ou_name ~ ex_name, value.var='p')
	pmat[is.na(pmat)] =0;  pmat=round(pmat,1); pmat[pmat >7] =7; 
pdf('mr.pval.pdf', w=20, h=20)
	corrplot(pmat, is.corr=F, method='shade', bg='black', col=colorRampPalette(c("white","green","gold"))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig="pch", pch.cex=2, tl.srt=45, outline=T)
	dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Circular Manhattan Plot 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- read.table('D:/data/gwas/bb/bb_ALP.p03', header=T) %>% dplyr::select(SNP, CHR, POS, P) %>% filter(!is.na(SNP)) %>% mutate(P=ifelse(P<1e-20,1e-20,P)) %>% rename(ALP=P)
	dat2 <- read.table('D:/data/gwas/bb/covid_severe.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.2=CHR, POS.2=POS, COVID_severe=all_inv_var_meta_p)
	dat3 <- read.table('D:/data/gwas/bb/covid_inpatient.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.3=CHR, POS.3=POS, COVID_inpatient=all_inv_var_meta_p)
	dat4 <- read.table('D:/data/gwas/bb/covid_infected.p03', header=T) %>% 
	dplyr::select(rsid, CHR, POS, all_inv_var_meta_p) %>% filter(!is.na(rsid)) %>% mutate(all_inv_var_meta_p=ifelse(all_inv_var_meta_p<1e-20,1e-20,all_inv_var_meta_p)) %>% rename(SNP=rsid, CHR.4=CHR, POS.4=POS, COVID_infected=all_inv_var_meta_p)
dat0 <- Reduce(function(x,y) merge(x,y,by='SNP',all=T, no.dups=T), list(dat1, dat2, dat3, dat4)) # !! Must Remove NA SNPs before merging
dat <- dat0 %>%
	mutate (chrom=apply( dat0[,grepl('CHR',names(dat0))], 1, FUN=min, na.rm=T ), position=apply( dat0[,grepl('POS',names(dat0))], 1, FUN=min, na.rm=T )) %>% 
	dplyr::select(SNP, chrom, position, COVID_severe, COVID_inpatient, COVID_infected, ALP)
	write.table(dat, "plot.dat.txt", quote=F, row.names=F, append=F)
SNPs <- subset(dat, (chrom==1 & position >=21835857 & position <=21904905) | (chrom==9 & position >=136130562 & position<=136150630))$SNP
CMplot(dat, type="p", plot.type="c", r=0.4, col=matrix(c("gold","orange", "gray","lightblue", "olivedrab3", "orange", "grey30","grey60"), nrow=4, byrow=T),
	highlight=SNPs, highlight.col="red",highlight.cex=1, highlight.pch=10,
	threshold=5e-8, cir.chr.h=1, amplify=F, threshold.col="red", signal.line=1, signal.col="black",
	bin.size=1e6, outward=T, file='jpg', memo='', dpi=600, file.output=T, verbose=T, width=12, height=12
)