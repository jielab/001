setwd("C:/Users/jiehu/Desktop")
pacman::p_load(readxl, dplyr, tidyverse, TwoSampleMR, MendelianRandomization, MVMR)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 最简单原始的方法
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
iv.snp <- read.table('D:/data/gwas/main/walk_pace.top.snp', header=T)
dat_X <- read.table('D:/data/gwas/main/walk_pace.gz', header=T) %>% merge(iv.snp, by="SNP") %>% 
	format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
dat_Y <- read.table('D:/data/gwas/main/y.vte.gz', header=T) %>%
	merge(iv.snp, by="SNP") %>% 
	format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') 
dat <- harmonise_data(dat_X, dat_Y, action=1) 
res <- mr(dat); res
	res_single <- mr_singlesnp(dat)
	mr(dat, method_list = c("mr_ivw", "mr_egger_regression", "mr_egger_regression_bootstrap"))
mr_pleiotropy_test(dat); mr_heterogeneity(dat)
mr_scatter_plot(res, dat)[[1]]
	mr_forest_plot(res_single)[[1]]
	mr_funnel_plot(res_single)
bsmr_dat <- dat_to_MRInput(dat); bsmr_dat <- bsmr_dat[[1]]
	mr_plot(bsmr_dat, interactive=F)
	mr_forest(bsmr_dat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于MR的mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
label = 'pheno'
ieu_X = ''; ieu_M = ''; ieu_Y = ''
dir_X = dir_M = dir_Y = 'D:/data/gwas/pheno'
XYs = as.data.frame(read_excel(paste0('D:/analysis/mr/', label,'.xlsx'))) %>% filter(!grepl("bb_", trait)); rownames(XYs) <- XYs$trait; XYs$trait <- NULL # EXCEL文件第一个格子应该是trait
X_list = rownames(XYs)
M_list = c('bb_CRE', 'bb_ALB', 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_LPA', 'bb_SHBG', 'bb_TES', 'bb_VITD')
Y_list = grep('^y.', names(XYs), value=T)
if (ieu_X !='') {Xs=ieu_X; ieu_X_use=TRUE} else {Xs=X_list; ieu_X_use=FALSE}
if (ieu_M !='') {Ms=ieu_M; ieu_M_use=TRUE} else {Ms=M_list; ieu_M_use=FALSE} 
if (ieu_Y !='') {Ys=ieu_Y; ieu_Y_use=TRUE} else {Ys=Y_list; ieu_Y_use=FALSE} 

sink("mediation.txt")
for (M in Ms) { # M
	writeLines(paste('\n\nRun:', M))
	if (ieu_M_use) {
		dat_M.snp <- extract_instruments(outcomes=M, clump=F) %>% select(SNP)
	} else {
		dat_M.raw <- read.table(paste0(dir_M, '/', M, '.gz'), header=T) 
		dat_M.snp <- dat_M.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
		dat_M.iv <- read.table(paste0(dir_M, '/', M, '.top.snp'), header=T); names(dat_M.iv) <- "SNP"  
	}
	
	for (X in Xs) { # X
		if (ieu_X_use) {
			dat_X.snp <- extract_instruments(outcomes=X, clump=F) %>% select(SNP)
		} else {
			dat_X.raw <- read.table(paste0(dir_X, '/', X, '.gz'), header=T)
			dat_X.snp <- dat_X.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
			dat_X.iv <- read.table(paste0(dir_X, '/', X, '.top.snp'), header=T); names(dat_X.iv) <- "SNP" 
		}
		dat_XnM.snp <- rbind(dat_X.iv, dat_M.iv) %>% unique() # 最理想的是把 dat_X.snp和dat_M.snp合并然后 %>% clump_data() %>% select(SNP)
		if (ieu_X_use) {
			dat_X <- extract_outcome_data(dat_XnM.snp, X) %>% convert_outcome_to_exposure()
		} else {
			dat_X <- dat_X.raw %>% merge(dat_XnM.snp) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		}
		if (ieu_M_use) {
			dat_M <- extract_outcome_data(dat_XnM.snp, M)
		} else {
			dat_M <- dat_M.raw %>% merge(dat_XnM.snp) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
		}
		dat_XM <- harmonise_data(dat_X, dat_M, action=1)
		# Step1: X --> M
		res_X2M <- mr(dat_XM, method_list=c("mr_wald_ratio", "mr_ivw")); res_X2M 
			beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- signif(res_X2M %>% pull(pval),2)
		names(dat_XM) <- gsub("outcome", "mediator", names(dat_XM))

		for (Y in Ys) {
			ivw_p <- as.numeric(XYs[X,Y]); if (is.na(ivw_p) | ivw_p > 1e-4) next
			writeLines(paste('\n\n\t-->Run:', X, M, Y))
			if (ieu_Y_use) {
				dat_Y <- extract_outcome_data(dat_XnM.snp, Y)
			} else {
				dat_Y.raw <- read.table(paste0(dir_Y, '/', Y, '.gz'), header=T)
				dat_Y <- merge(dat_Y.raw, dat_XnM.snp) %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)
			}
			dat_XY <- harmonise_data(dat_X, dat_Y, action=1) 
			# Total effect
			res_X2Y <- mr(dat_XY, method_list=c("mr_wald_ratio", "mr_ivw")); res_X2Y 
				beta_X2Y <- res_X2Y %>% pull(b); se_X2Y <- res_X2Y %>% pull(se); p_X2Y <- signif(res_X2Y %>% pull(pval),2)
			
			# Merge X, M, Y 
			dat <- merge(dat_XM, dat_XY, by="SNP") %>% drop_na(beta.exposure.x, beta.mediator, beta.outcome)
			bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
			if (bad_row !=0) {print(paste(X, Y, M, "dat_XM与dat_XY中的 effect_allele.exposure不同")); next}
			names(dat) <- gsub("\\.x$", "", names(dat))
			dat <- subset(dat, select=!grepl("\\.y", names(dat)))
			if (nrow(dat)==0) next
			# Step 2: M (+X) --> Y, 使用 MVMR from WSpiller Package
			dat1 <- format_mvmr(RSID=dat$SNP, BXGs=subset(dat, select=c(beta.mediator, beta.exposure)), BYG=dat$beta.outcome, seBXGs=subset(dat, select=c(se.mediator, se.exposure)), seBYG=dat$se.outcome)
				res_mvmr <- ivw_mvmr(r_input=dat1); res_mvmr
				beta_M2Y.adjX=res_mvmr[1,1]; se_M2Y.adjX=res_mvmr[1,2]; p_M2Y.adjX=signif(res_mvmr[1,4],2)
			
			# Mediation 乘法计算
			beta_2step <- round(beta_X2M * beta_M2Y.adjX, 4)
				CIs = RMediation::medci(beta_X2M, beta_M2Y.adjX, se_X2M, se_M2Y.adjX, type='dop'); se_2step = CIs$SE; p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F),2)		
			print(paste(X,M,Y, round(beta_X2Y,3),round(se_X2Y,3),p_X2Y, round(beta_X2M,3),round(se_X2M,3),p_X2M, round(beta_M2Y.adjX,3),round(se_M2Y.adjX,3),p_M2Y.adjX, round(beta_2step,3),round(se_2step,3),p_2step))
		}	
	}
}
sink()