setwd("C:/Users/jiehu/Desktop")
pacman::p_load(data.table, lubridate, tidyverse, dplyr, ggplot2, CMplot, TwoSampleMR, MendelianRandomization, survival, survminer, mediation)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
table(dat0$sp1); table(dat0$blood_group); table(dat0$abo1); hist(dat0$bb_ALP)
dat <- dat0 %>% filter(ethnicity_gen==1) %>% 
	mutate (
	o_no = ifelse(blood_group !="OO", 1,0),
	sp1_dose = 
	inpatient_yes = ifelse(!is.na(icd_covid_date), 1, ifelse(covid_inf==1,0, NA))
	)
round(prop.table(table(dat$abo1, dat$covid_inf), 1), 5) # sp1, age_cat
tab1 <- dat %>% group_by(abo1) %>% summarise_at(vars(grep("bb_", names(dat), value=T)), list(mean=mean), na.rm=T) %>% as.data.frame() %>% na.omit()
	tab1.long <- melt(tab1, id.vars="abo1")
	ggplot(tab1.long, aes(x=variable, y=value, fill=factor(abo1))) + geom_bar(stat="identity", position="dodge") + theme(axis.text.x = element_text(angle = 90))
naniar::gg_miss_var(subset(dat, select=grep("rheu|oest", grep("sex|bb_", names(dat), value=T), invert=T, value=T)), facet=sex)
for (varY in grep("^bb_|^bc_", names(dat), value=T)) {
	# dat[[var]] = inormal(dat[[var]])
	dat$varY = dat[[varY]]
	print(noquote(varY))
	for (varX in grep("^sp1\\.", names(dat), value=T)) {
		dat$varX <- dat[[varX]]
		fit.lm <- lm(varY ~ varX +age+sex, data=dat)
		print(noquote(paste("  --", varX, signif( coef(summary(fit.lm))[2,4], 3))))
	}
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mediation 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
varX="abo1"; varY="inpatient_yes"; varM="bb_ALP"
dat1 <- subset(dat, select=c(varX, varY, varM, "age", "sex", "bmi")) %>% na.omit()
names(dat1) <- c("varX", "varY", "varM", "age","sex", "bmi")
fit.med = lm(varM ~ varX, data=dat1); summary(fit.med)	
fit.out = glm(varY ~ varX + varM, data=dat1, family=binomial("probit")); summary(fit.out) 
med.out = mediation::mediate(fit.med, fit.out, treat="varX", mediator="varM", control.value="O", treat.value="A", sims=500, boot=T)
	print(summary(med.out)) # plot(med.out)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mendelian Randomization
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
exp="ALP" # ("ALB", "ALP", "ALT", "APOA", "APOB", "AST", "AST2ALT", "BILD", "BUN", "CA", "CHOL", "CRE", "CRP", "CYS", "EGFR", "GGT", "GLU", "HBA1C", "HDL", "IGF1", "LDLD", "LPA", "NAP", "PHOS", "SHBG", "TBIL", "TES", "TP", "TRIG", "UA", "UCR", "URMA", "URNA", "VITD")
	exp_file=paste0("D:/data/gwas/bb/bb_", exp, ".top.txt")
	exp_phe="Trait"; exp_snp="SNP"; exp_ea="EA"; exp_nea="NEA"; exp_beta="BETA"; exp_se="SE"; exp_p="P"
	exp_dat <- read_exposure_data(filename=exp_file, phenotype_col=exp_phe, snp_col=exp_snp, effect_allele_col=exp_ea, other_allele_col=exp_nea, beta_col=exp_beta, se_col=exp_se, pval_col=exp_p) # %>% filter (pval.exposure <1e-3)
out="covid_severe"
	out_dat <- read_outcome_data(filename=out_file, sep="\t", snps=exp_dat$SNP, snp_col=out_snp, effect_allele_col=out_ea, other_allele_col=out_nea, beta_col=out_beta, se_col=out_se, pval_col=out_p) 
	out_file="D:/data/gwas/bb/covid_severe.small" # can NOT be in .gz format
	out_snp="rsid"; out_ea="ALT"; out_nea="REF"; out_beta="all_inv_var_meta_beta"; out_se="all_inv_var_meta_sebeta"; out_p="all_inv_var_meta_p"
dat <- harmonise_data(exp_dat, out_dat)
	res <- mr(dat, method_list=c('mr_ivw', 'mr_simple_median', 'mr_weighted_median', 'mr_two_sample_ml', 'mr_penalised_weighted_median', 'mr_egger_regression')); print(res)
	plt <- mr_scatter_plot(res, dat)
		#plt + theme(legend.text=element_text(size=14, face="bold"), axis.title=element_text(size=10, face="bold"), axis.text=element_text(size=15, face="bold"))
mrdat <- mr_input(dat$beta.exposure, dat$se.exposure, dat$beta.outcome, dat$se.outcome); # mr_funnel(mrdat) 
	plt <- mr_plot(mrdat, interactive=F)
		plt + theme(axis.title=element_text(size=15, face="bold"), axis.text=element_text(size=12, face="bold"))
	plt <- mr_forest(mrdat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) 
		plt + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))
	plt <- mr_funnel(mrdat)
		plt + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=14), axis.text=element_text(size=12, face="bold"))	
	ggsave(plt, file=paste0(exp,".exp.png"), w=10, h=10)
## reverse-direction
exp_snp="rsid"; exp_ea="ALT"; exp_nea="REF"; exp_beta="all_inv_var_meta_beta"; exp_se="all_inv_var_meta_sebeta"; exp_p="all_inv_var_meta_p"
out_phe="Trait"; out_snp="SNP"; out_ea="EA"; out_nea="NEA"; out_beta="BETA"; out_se="SE"; out_p="P"
exp_file="D:/data/gwas/covid_sus.gwas"
exp_dat <- read_exposure_data(filename=exp_file, sep="\t", snp_col=exp_snp, effect_allele_col=exp_ea, other_allele_col=exp_nea, beta_col=exp_beta, se_col=exp_se, pval_col=exp_p) # %>% filter (pval.exposure <1e-3)
exp_dat <- extract_instruments(exp_dat, p1=5e-08, clump=T, r2=0.1, kb=1000, force_server=F) 
out_file=paste0("D:/data/gwas/", out, ".top.txt")
out_dat <- read_exposure_data(filename=out_file, phenotype_col=out_phe, snp_col=out_snp, effect_allele_col=out_ea, other_allele_col=out_nea, beta_col=out_beta, se_col=out_se, pval_col=out_p) # %>% filter (pval.exposure <1e-3)
dat <- harmonise_data(exp_dat, out_dat)
res <- mr(dat, method_list=c('mr_ivw', 'mr_simple_median', 'mr_weighted_median', 'mr_two_sample_ml', 'mr_penalised_weighted_median', 'mr_egger_regression')) 
pp <- mr_scatter_plot(res, dat); ggsave(pp[[1]], file=paste0(out,".out.png"), w=7, h=7)
## Summarize Many MR results
# pdftk */*.cB.pdf cat output 001.cB.pdf
# cat */*.mr.res.txt | awk 'NF ==10' > 001.mr.txt
pacman::p_load(data.table, dplyr, tidyverse, magrittr, corrplot, reshape2)
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
# Manhattan Plot 
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