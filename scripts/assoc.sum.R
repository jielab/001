pacman::p_load(tidyr, TwoSampleMR, MendelianRandomization, survival, survminer)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于MR的mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.Y <- dir.X = 'D:/data/gwas/main'
dir.M <- "D:/data/gwas/bb"
Xs <- "?"
Ys <- "?" # grep('^y.', names(XYs), value=T)
Ms <- "?" # c('bb_CRE', 'bb_ALB', 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_LPA', 'bb_SHBG', 'bb_TES', 'bb_VITD')
sink("?.log")

for (Y in Ys) {
	writeLines(paste('\n\n\t-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)

	for (X in Xs) { # X
		dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T) 
		if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
		} else {
			dat.X.iv <- dat.X.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
		}
		
		for (M in Ms) { # M
			writeLines(paste('\t\tRun:', Y, X, M))
			dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T) 
			if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
				dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- "SNP" 
			} else {
				dat.M.iv <- dat.M.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
			}
			
			# 合并 dat.X 和 dat.M 
			dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
			dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
			dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
			dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
			dat.M.4X <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
			# Step1: X --> M
			dat.XM <- harmonise_data(dat.X, dat.M.4X, action=1)
			dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
			res.X2M <- mr(dat.XM, method_list=c("mr_wald_ratio", "mr_ivw")); res.X2M 
			beta.X2M <- res.X2M %>% pull(b); se.X2M <- res.X2M %>% pull(se); p.X2M <- signif(res.X2M %>% pull(pval),2)
	
			# 合并 dat.X 和 dat.M 和 dat.Y
			dat.Y <- merge(dat.Y.raw, dat.XnM.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
			dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y, action=1)
			dat <- merge(dat.XM.mv, dat.XY.mv, by="SNP")
			bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
			if (bad_row !=0) {print(paste("ERROR:", X, Y, M, "dat.XM与dat.XY中的 effect_allele.exposure不同")); next}
			names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
			if (nrow(dat)==0) next
			# Total effect
			dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
			res.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); res.X2Y 
			beta.X2Y <- res.X2Y %>% pull(b); se.X2Y <- res.X2Y %>% pull(se); p.X2Y <- signif(res.X2Y %>% pull(pval),2)	
			
			# Step 2: M (+X) --> Y 
			dat.mvmr <- MendelianRandomization::mr_mvinput(bx=cbind(dat$beta.mediator, dat$beta.exposure), bxse=cbind(dat$se.mediator, dat$se.exposure), by=dat$beta.outcome, byse=dat$se.outcome)
			res.mvmr <- MendelianRandomization::mr_mvivw(dat.mvmr, model="default", correl=FALSE, distribution="normal", alpha=0.05); res.mvmr
			beta.M2Y <- as.numeric(res.mvmr$Estimate[1]); se.M2Y <- as.numeric(res.mvmr$StdError[1]); p.M2Y <- as.numeric(res.mvmr$Pvalue[1])
			
			# Mediation
			beta.me <- beta.X2M * beta.M2Y # 🐕 乘法 
			se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2) # 🐕 POE法，也称Delta法 
			p.me <- 2*pnorm(-abs(beta.me/se.me))
			print(paste("RES:", X,M,Y, round(beta.X2Y,3),round(se.X2Y,3),p.X2Y, round(beta.X2M,3),round(se.X2M,3),p.X2M, round(beta.M2Y,3),round(se.M2Y,3),p.M2Y, round(beta.me,3),round(se.me,3),p.me))

		}	
	}
}
sink()