setwd("C:/Users/jiehu/Desktop")
pacman::p_load(readxl, dplyr, tidyverse, TwoSampleMR, MVMR)

label = 'pheno'
dir_X = dir_M = dir_Y = 'D:/data/gwas/pheno'
XYs = as.data.frame(read_excel(paste0('D:/analysis/mr/', label,'.xlsx'))); rownames(XYs) <- XYs$trait; XYs$trait <- NULL # EXCEL文件第一个格子应该是trait
Xs = 'x.height' # rownames(XYs)
Ms = 'bb_CRE' # c('bb_ALB', 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_LPA', 'bb_SHBG', 'bb_TES', 'bb_VITD')
Ys = 'y.vte' # grep('^y.', names(XYs), value=T)

for (M in Ms) { # M
	writeLines(paste('\n\nRun:', M))
	dat_M.raw <- read.table(paste0(dir_M, '/', M, '.gz'), header=T) 
	dat_M.snp <- dat_M.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
	dat_M.iv <- read.table(paste0(dir_M, '/', M, '.top.snp'), header=T); names(dat_M.iv) <- "SNP"  

	for (X in Xs) { # X
		dat_X.raw <- read.table(paste0(dir_X, '/', X, '.gz'), header=T)
		dat_X.snp <- dat_X.raw %>% filter(P<5e-8) %>% dplyr::select(SNP)
		dat_X.iv <- read.table(paste0(dir_X, '/', X, '.top.snp'), header=T); names(dat_X.iv) <- "SNP" 
		dat_XnM.snp <- rbind(dat_X.iv, dat_M.iv) %>% unique() # 最理想的是把 dat_X.snp和dat_M.snp合并然后 %>% clump_data() %>% select(SNP)
		dat_X <- dat_X.raw %>% merge(dat_XnM.snp) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		dat_M <- dat_M.raw %>% merge(dat_XnM.snp) %>% format_data(type='outcome', phenotype_col='exposure', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
		dat_XM <- harmonise_data(dat_X, dat_M)
		names(dat_XM) <- gsub("outcome", "mediator", names(dat_XM))

		for (Y in Ys) { # Y
			ivw_p <- XYs[X,Y]; if (!is.na(ivw_p) & ivw_p > 1e-4) next
			writeLines(paste('\n\nRun:', X, M, Y))
			dat_Y.raw <- read.table(paste0(dir_Y, '/', Y, '.gz'), header=T)
			dat_Y <- merge(dat_Y.raw, dat_XnM.snp) %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)
			dat_XY <- harmonise_data(dat_X, dat_Y) 
			
			# merged data
			dat0 <- merge(dat_XM, dat_XY, by="SNP") %>% drop_na(beta.exposure.x, beta.mediator, beta.outcome); # subset(dat_XMY_tmp, effect_allele.exposure.x != effect_allele.exposure.y)
			names(dat0) <- gsub("\\.x$", "", names(dat0))
			if (nrow(dat0)==0) next

			# Total effect
			res_X2Y <- mr(dat0) %>% filter(method=='Wald ratio' | method=='Inverse variance weighted'); res_X2Y 
				beta_X2Y <- res_X2Y %>% pull(b); se_X2Y <- res_X2Y %>% pull(se); p_X2Y <- signif(res_X2Y %>% pull(pval),2)
			# step1: M ~ X
			dat <- dat0; names(dat) <- gsub(".mediator$", ".outcome", gsub(".outcome$", ".outCOPY", names(dat0)))
				res_X2M <- mr(dat) %>% filter(method=='Wald ratio' | method=='Inverse variance weighted'); res_X2M 
				beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- signif(res_X2M %>% pull(pval),2)
			# step 2: Y ~ M + X, 使用 MVMR from WSpiller Package
			dat <- format_mvmr(RSID=dat0$SNP, BXGs=subset(dat0, select=c(beta.mediator, beta.exposure)), BYG=dat0$beta.outcome, seBXGs=subset(dat0, select=c(se.mediator, se.exposure)), seBYG=dat0$se.outcome)
				res_mvmr <- ivw_mvmr(r_input=dat); res_mvmr
				beta_M2Y.adjX=res_mvmr[1,1]; se_M2Y.adjX=res_mvmr[1,2]; p_M2Y.adjX=signif(res_mvmr[1,4],2)
			# mediation 乘法计算
			beta_2step <- round(beta_X2M * beta_M2Y.adjX, 4)
				CIs = RMediation::medci(beta_X2M, beta_M2Y.adjX, se_X2M, se_M2Y.adjX, type='MC'); se_2step = CIs$SE; p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F),2)		
			print(paste(X,M,Y, round(beta_X2Y,3),round(se_X2Y,3),p_X2Y, round(beta_X2M,3),round(se_X2M,3),p_X2M, round(beta_M2Y.adjX,3),round(se_M2Y.adjX,3),p_M2Y.adjX, round(beta_2step,3),round(se_2step,3),p_2step))
		}	
	}
}
