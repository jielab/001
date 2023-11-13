setwd("D:/")
pacman::p_load(readxl, tidyverse, TwoSampleMR, RMediation)

analysis <- "mb"
dir_mr <- paste0("D:/analysis/mr/", analysis) # MR分析结果所在地
dir_M <- paste0("D:/data/gwas/", "pheno") # mediator的GWAS和top.snp所在地
dir_Y <- paste0("D:/data/gwas/", "pheno") #
xys <- as.data.frame(read_excel(paste0(dir_mr,".xlsx"))); rownames(xys) <- xys$variable; xys$variable <- NULL
for (m in c('bb_ALB')){ #, 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_AST', 'bb_AST2ALT', 'bb_BILD', 'bb_BUN', 'bb_CA', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_GGT', 'bb_GLU', 'bb_HBA1C', 'bb_HDL', 'bb_IGF1', 'bb_LDLD', 'bb_LPA', 'bb_NAP', 'bb_PHOS', 'bb_SHBG', 'bb_TBIL', 'bb_TES', 'bb_TP', 'bb_TRIG', 'bb_UA', 'bb_UCR', 'bb_URK', 'bb_URMA', 'bb_URNA', 'bb_VITD')) {
	writeLines(paste('\n\nRun:', m))
	dat_M0 <- read.table(paste0(dir_M, "/", m, ".gz"), header=T)
	snp_M <- read.table(paste0(dir_M, "/", m, ".top.snps"), header=T) # dir_M 或者 dir_mr
	dat_M4y <- merge(dat_M0, snp_M) %>% format_data(type="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P")
	for (x in rownames(xys)) {
	for (y in grep("^y.", names(xys), value=T)) { # invert=T
		ivw_p <- xys[x,y]
		if (ivw_p <= 1e-6) {
			writeLines(paste('\n\nRun:', x, m, y))
			dat_X <- read.table(paste0(dir_mr, "/", x, "/", x, ".txt"), header=T) %>% format_data(type ="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_Y0 <- read.table(paste0(dir_Y, "/", y, ".gz"), header=T) 
			dat_M4x <- dat_M0 %>% merge(subset(dat_X, select="SNP"))
			dat_Y4x <- dat_Y0 %>% merge(subset(dat_X, select="SNP")) 
			dat_Y4m <- dat_Y0 %>% merge(subset(dat_M4y, select="SNP")) 
			if (nrow(dat_M4x) * nrow(dat_Y4x) * nrow(dat_Y4m) ==0) {
				write.table(paste(x, m, y, nrow(dat_M4x), nrow(dat_Y4x), nrow(dat_Y4x), "NA ... ...", sep='|'), paste0(analysis,".out"), append=T, quote=F, row.names=F, col.names=F)
				next
			}
			dat_M4x <- format_data(dat_M4x, type ="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_Y4x <- format_data(dat_Y4x, type="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P")
			dat_Y4m <- format_data(dat_Y4m, type="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_XM <- harmonise_data(dat_X, dat_M4x, action=2) 
			dat_MY <- harmonise_data(dat_M4y, dat_Y4m, action=2) 
			dat_XY <- harmonise_data(dat_X, dat_Y4x) 
			res_X2M <- TwoSampleMR::mr(dat_XM) %>% filter(method=="Wald ratio" | method=="Inverse variance weighted") # 1st step 三角形的上坡
				beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- res_X2M %>% pull(pval)
			res_M2Y <- TwoSampleMR::mr(dat_MY) %>% filter(method=="Wald ratio" | method=="Inverse variance weighted") # 1st step 三角形的下坡
				beta_M2Y <- res_M2Y %>% pull(b); se_M2Y <- res_M2Y %>% pull(se); p_M2Y <- res_M2Y %>% pull(pval)
			beta_2step <- round(beta_X2M * beta_M2Y, 4); beta_2step
				CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y, type='MC'); se_2step = CIs$SE 
				p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F), 3); p_2step
			res_str <- paste(x, m, y, nrow(dat_XM), nrow(dat_MY), nrow(dat_XY), beta_X2M, beta_M2Y, beta_2step, se_2step, p_2step, sep='|')
			write.table(res_str, paste0(analysis,".out"), append=T, quote=F, row.names=F, col.names=F)
		}
	}
	}
}
