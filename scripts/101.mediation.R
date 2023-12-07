setwd("C:/Users/jiehu/Desktop")
pacman::p_load(readxl, dplyr, tidyverse, TwoSampleMR, MVMR)
label = 'height'
dir_X = 'D:/data/gwas/pheno'
dir_M = 'D:/data/gwas/pheno'
dir_Y = 'D:/data/gwas/pheno' 
X_ieus = ''  # BMI
M_ieus = '' # CRP
Y_ieus = '' # CAD 
# XYs = as.data.frame(read_excel(paste0('D:/analysis/mr/', label,'.xlsx'))); rownames(XYs) <- XYs$variable; XYs$variable <- NULL
X_list = 'x.height' # rownames(XYs)
M_list = 'bb_CRE' # c('bb_ALB', 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_LPA', 'bb_SHBG', 'bb_TES', 'bb_VITD')
Y_list = 'y.vte' # grep('^y.', names(XYs), value=T)

if (X_ieus !='') {Xs=X_ieus; X_use_ieu=TRUE} else {Xs=X_list; X_use_ieu=FALSE}
if (M_ieus !='') {Ms=M_ieus; M_use_ieu=TRUE} else {Ms=M_list; M_use_ieu=FALSE} 
if (Y_ieus !='') {Ys=Y_ieus; Y_use_ieu=TRUE} else {Ys=Y_list; Y_use_ieu=FALSE} 

for (M in Ms) { # M
	writeLines(paste('\n\nRun:', M))
	# The following codes create dat_M1, the non-clumped version of dat_M
	if (M_use_ieu) {
		dat_M1 <- extract_instruments(outcomes=M, clump=F) %>% dplyr::select(-c(samplesize.exposure, data_source.exposure))
	} else {
		dat_M0 <- read.table(paste0(dir_M, '/', M, '.gz'), header=T) %>% mutate(exposure=M)
		IV_M <- read.table(paste0(dir_M, '/', M, '.top.snp'), header=T); names(IV_M) <- "SNP"  
		dat_M1 <- dat_M0 %>% merge(IV_M) %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=M)
	}

	for (X in Xs) { # X
		# The following codes create dat_X1, the non-clumped version of dat_X
		if (X_use_ieu) {
			dat_X1 <- extract_instruments(outcomes=X, clump=F) %>% dplyr::select(-c(samplesize.exposure, data_source.exposure))
			dat_X1.clumped <- dat_X1 %>% clump_data()
		} else {
			dat_X0 <- read.table(paste0(dir_X, '/', X, '.gz'), header=T)
			IV_X <- read.table(paste0(dir_X, '/', X, '.top.snp'), header=F); names(IV_X) <- "SNP" 
			dat_X1 <- dat_X0 %>% merge(IV_X) %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
			dat_X1.clumped <- dat_X1
		}
		if (M_use_ieu) {
			dat_M4x <- extract_outcome_data(dat_X1.clumped$SNP, M)
		} else {
			dat_M4x <- dat_M0 %>% merge(subset(dat_X1.clumped, select='SNP')) %>% format_data(type='outcome', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
		}
		dat_X2M <- harmonise_data(dat_X1.clumped, dat_M4x)
		# The following codes merge X and M, then clump, then create: dat_X, dat_M, dat_X8M
		dat_X8M1 <- rbind(dat_X1, dat_M1) %>% clump_data() # clump takes time! may use unique() as a shortcut.
		if (X_use_ieu) {
			dat_X <- extract_outcome_data(dat_X8M1$SNP, X) %>% convert_outcome_to_exposure()
		} else {
			dat_X <- dat_X0 %>% merge(subset(dat_X8M1, select='SNP')) %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		}
		if (M_use_ieu) {
			dat_M <- extract_outcome_data(dat_X8M1$SNP, M) %>% convert_outcome_to_exposure()
		} else {
			dat_M <- dat_M0 %>% merge(subset(dat_X8M1, select='SNP')) %>% format_data(type='exposure', phenotype_col='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure = M)
		}
		dat_X8M <- rbind(dat_X, dat_M)

		for (Y in Ys) { # Y
			# ivw_p <- XYs[X,Y]; if (!is.na(ivw_p) & ivw_p > 1e-6) next
			writeLines(paste('\n\nRun:', X, M, Y))
			if (Y_use_ieu) {
				dat_Y <- extract_outcome_data(dat_X8M$SNP, Y)
			} else {
				dat_Y <- read.table(paste0(dir_Y, '/', Y, '.gz'), header=T) %>% merge(subset(dat_X8M, select='SNP')) %>% format_data(type='outcome', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)
			}
			if (nrow(dat_Y)==0) next

			dat <- harmonise_data(dat_X8M, dat_Y)
			# Total effect
			res_X2Y <- dat %>% filter(id.exposure==X) %>% mr() %>% filter(method=='Inverse variance weighted') 
				beta_X2Y <- res_X2Y %>% pull(b); se_X2Y <- res_X2Y %>% pull(se); p_X2Y <- signif(res_X2Y %>% pull(pval),2)
			# Two-step MR 
			res_X2M <- mr(dat_X2M) %>% filter(method=='Wald ratio' | method=='Inverse variance weighted') # 1st step 三角形的上坡
				beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- res_X2M %>% pull(pval)
			res_M2Y <- dat %>% filter(id.exposure==M) %>% mr() %>% filter(method=='Wald ratio' | method=='Inverse variance weighted') # 2nd step 三角形的下坡
				beta_M2Y <- res_M2Y %>% pull(b); se_M2Y <- res_M2Y %>% pull(se); p_M2Y <- res_M2Y %>% pull(pval)
			beta_2step <- round(beta_X2M * beta_M2Y, 4); beta_2step # product method
				CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y, type='MC'); se_2step = CIs$SE; p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F), 3); p_2step			
			print(paste(X, M, Y, p_X2Y, p_2step))
			
			# MVMR from TwoSampleMR Package
			dat <- mv_harmonise_data(dat_X8M, dat_Y)
			mv_multiple(dat) 
			# MVMR from WSpiller Package
			dat <- format_mvmr(RSID=rownames(dat$exposure_beta), BXGs=dat$exposure_beta, BYG=dat$outcome_beta, seBXGs=dat$exposure_se, seBYG=dat$outcome_se)
				ivw_mvmr(r_input=dat)
				strength_mvmr(r_input=dat) 
				pleiotropy_mvmr(r_input=dat) 
		}
	
	}
}
