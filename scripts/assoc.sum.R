pacman::p_load(data.table, stringi, dplyr, tidyr, TwoSampleMR, MendelianRandomization, survival, survminer)

pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS','EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于MR的mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir.Y = dir.X = '/data/sph-huangj/data/gwas/main'
dir.M = '/data/sph-huangj/data/gwas/?/clean'
Xs = '?'
Ys = '?' 
Ms = gsub('.gz', '', list.files(path=dir.M, pattern='.gz$')); Ms
sink("?.log")

for (Y in Ys) { # 🙍
	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

	for (X in Xs) { # 🍷
		dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
		names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
		if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
		} else {
			dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+06))
			dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
		}
		dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		# Total effect
		dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
		res.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); res.X2Y 
		beta.X2Y <- res.X2Y$b; se.X2Y <- res.X2Y$se; p.X2Y <- signif(res.X2Y$pval,2)
		if(p.X2Y >0.05) {print(paste("RES: XY不显著|", X,Y)); next} # 🛑

		for (M in Ms) { # 🐎
			writeLines(paste('\n\tRun:', Y, X, M))
			dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T)
			names(dat.M.raw) <- stri_replace_all_regex(toupper(names(dat.M.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
			if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
				dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- "SNP" 
			} else {
				dat.M.sig <- dat.M.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+06))
				dat.M.iv <- dat.M.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
			}
			# 合并 dat.X 和 dat.M 
			dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
			dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
			dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
			dat.M.4X <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
			# Step1: X --> M
			dat.XM <- harmonise_data(dat.X, dat.M.4X, action=1)
			dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
			res.X2M <- mr(dat.XM, method_list=c("mr_wald_ratio", "mr_ivw")); res.X2M 
			beta.X2M <- res.X2M$b; se.X2M <- res.X2M$se; p.X2M <- signif(res.X2M$pval,2)
	
			# 合并 dat.X 和 dat.M 和 dat.Y
			dat.Y <- merge(dat.Y.raw, dat.XnM.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
			dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y, action=1)
			dat <- merge(dat.XM.mv, dat.XY.mv, by="SNP")
			bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
			if (bad_row !=0) {print(paste("ERROR:", X, Y, M, "dat.XM与dat.XY中的 effect_allele.exposure不同")); next}
			names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
			if (nrow(dat)==0) next
			
			# Step 2: M (+X) --> Y 
			dat.mvmr <- MendelianRandomization::mr_mvinput(bx=cbind(dat$beta.mediator, dat$beta.exposure), bxse=cbind(dat$se.mediator, dat$se.exposure), by=dat$beta.outcome, byse=dat$se.outcome)
			res.mvmr <- MendelianRandomization::mr_mvivw(dat.mvmr, model="default", correl=FALSE, distribution="normal", alpha=0.05); res.mvmr
			beta.M2Y <- res.mvmr$Estimate[["Bx1"]]; se.M2Y <- res.mvmr$StdError[["Bx1"]] ; p.M2Y <- res.mvmr$Pvalue[["Bx1"]]
			
			# Mediation
			beta.me <- beta.X2M * beta.M2Y # 🐕 乘法
			prop.me <- beta.me / beta.X2Y
			se.me <- sqrt(beta.X2M^2 * se.X2M^2 + beta.M2Y^2 * se.M2Y^2) # 🐕 POE法，也称Delta法 
			p.me <- 2*pnorm(-abs(beta.me/se.me))
			print(paste("RES: XY|", X,M,Y, nrow(dat.XY),round(beta.X2Y,3),round(se.X2Y,3),p.X2Y, "|", nrow(dat.XM),round(beta.X2M,3),round(se.X2M,3),p.X2M, "|", nrow(dat),round(beta.M2Y,3),round(se.M2Y,3),signif(p.M2Y,2), "|", round(prop.me,3),round(beta.me,3),round(se.me,3),signif(p.me,2)))

		}	
	}
}
sink()
