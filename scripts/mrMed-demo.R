# devtools::install_github("scllin/mrMed") 
# https://rdrr.io/github/scllin/toypack/src/R/mrMed.R
pacman::p_load(data.table, stringi, tidyverse, OneSampleMR, TwoSampleMR, MendelianRandomization, mrMed, RadialMR, psych, bda)
pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))


dir.Y = 'D:/data/gwas/main/clean'; Y = 'covid_C2' # 🙍
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.EUR.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
dir.X = dir.Y; X = 'ABO' # 🍷
	dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
	dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
# Total effect
dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
	fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2Y 
	beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- rp(fit.X2Y$pval)
dir.M = 'D:/data/gwas/bb/clean'; M = 'bb_ALP' # 🐎
	dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T)
	names(dat.M.raw) <- stri_replace_all_regex(toupper(names(dat.M.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- "SNP" 
	} else {
		dat.M.sig <- dat.M.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.M.iv <- dat.M.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
	}
# Step1: X --> M # 
	dat.M.4X <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
	dat.XM <- harmonise_data(dat.X, dat.M.4X, action=1)
	fit.X2M <- mr(dat.XM, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2M 
	beta.X2M <- fit.X2M$b; se.X2M <- fit.X2M$se; p.X2M <- rp(fit.X2M$pval)
# Harmonize dat.X + dat.M + dat.Y
	dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by="SNP")
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
# Step 2: M (+X) --> Y 
	dat.M2Y <- MendelianRandomization::mr_mvinput(bx=cbind(dat$beta.mediator, dat$beta.exposure), bxse=cbind(dat$se.mediator, dat$se.exposure), by=dat$beta.outcome, byse=dat$se.outcome)
	fit.M2Y <- MendelianRandomization::mr_mvivw(dat.M2Y, model="default", correl=FALSE, distribution="normal", alpha=0.05); fit.M2Y
	beta.M2Y <- fit.M2Y$Estimate[["Bx1"]]; se.M2Y <- fit.M2Y$StdError[["Bx1"]] ; p.M2Y <- fit.M2Y$Pvalue[["Bx1"]]
	beta.X2Y.dr <- fit.M2Y$Estimate[["Bx2"]]; se.X2Y.dr <- fit.M2Y$StdError[["Bx2"]] ; p.X2Y.dr <- fit.M2Y$Pvalue[["Bx2"]]
# mediation self-made method # 🏮
	beta.me <- beta.X2M * beta.M2Y # 🐕 product
	se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2); p.me <- 2*pnorm(-abs(beta.me/se.me))
	beta.dr <- beta.X2Y - beta.me # difference
	se.dr <- sqrt(se.X2Y^2 + se.me^2); p.dr <- 2*pnorm(-abs(beta.dr/se.dr))
	beta.prop <- beta.me / beta.X2Y
	se.prop <- sqrt(se.me^2/beta.X2Y^2 + beta.me^2*se.X2Y^2/beta.X2Y^4); p.prop <- 2*pnorm(-abs(beta.prop/se.prop))
# mrMed method # 🏮
	dat1 <- dat
	names(dat1) <- stri_replace_all_regex(names(dat1), pattern=c("exposure", "mediator", "outcome"), replacement=c("X", "M", "Y"), vectorize_all=FALSE)
	dat1 <- dat1 %>% mutate (
		Gx = ifelse(pval.X >5e-08 | pval.M >5e-08, 0, 1), # 仅仅用于 *_0 的 function 里面
		Gx_plum = Gx, # (Gx==1 & [XY]keep=TRUE), 用于 X -> M或Y
		Gm = ifelse(pval.M <5e-08, 1, 0), # (... & [XM]keep=TRUE), 没用 
		Gm_plum = Gm, # (... & [MY]keep=TRUE), 用于 M -> Y
		G_mvmr = ifelse(Gx_plum==0 & Gm_plum==0, 0, 1)
	)
	res <- mrMed(dat_mrMed=dat1, method_list="Prod_IVW") # 还有 Diff_IVW, Prod_Median
# key results # ✒
	X_str=paste0(X, "(", nrow(dat.X.iv), ")"); M_str=paste0(M, "(", nrow(dat.M.iv), ")");
	print(paste("RES:", X_str, M_str, Y, "|[X2Y]", nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		"[X2M]", nrow(dat.XM), rb(beta.X2M), rp(p.X2M), rb(beta.dr), rp(p.dr),
		"[M2Y]", nrow(dat), rb(beta.M2Y), rp(p.M2Y), 
		"[me]", rb(beta.me), rp(p.me), rb(beta.X2Y.dr), rp(p.X2Y.dr), rb(beta.prop), rp(p.prop), # ACME, ADE, Prop(rho) 
		"[mrMed]", rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p) # TE, ACME, ADE, Prop(rho) 
	))
