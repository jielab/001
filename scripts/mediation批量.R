setwd("D:/")
pacman::p_load(readxl, tidyverse, TwoSampleMR)

dir_gwas <- "D:/data/gwas/pheno/"
dir_mr <- "D:/analysis/mr/pheno/"
xys <- as.data.frame(read_excel("D:/analysis/mr/pheno.xlsx")); rownames(xys) <- xys$variable; xys$variable <- NULL
for (m in c('bb_ALB', 'bb_ALP', 'bb_ALT', 'bb_APOA', 'bb_APOB', 'bb_AST', 'bb_AST2ALT', 'bb_BILD', 'bb_BUN', 'bb_CA', 'bb_CHOL', 'bb_CRE', 'bb_CRP', 'bb_CYS', 'bb_EGFR', 'bb_GGT', 'bb_GLU', 'bb_HBA1C', 'bb_HDL', 'bb_IGF1', 'bb_LDLD', 'bb_LPA', 'bb_NAP', 'bb_PHOS', 'bb_SHBG', 'bb_TBIL', 'bb_TES', 'bb_TP', 'bb_TRIG', 'bb_UA', 'bb_UCR', 'bb_URK', 'bb_URMA', 'bb_URNA', 'bb_VITD')) {
	print(paste('Run:', m))
	dat_M0 <- read.table(paste0(dir_gwas, m, ".gz"), header=T) # raw
	snp_M <- read.table(paste0(dir_gwas, m, ".top.snps"), header=T)
	dat_M4y <- merge(dat_M0, snp_M) %>% format_data(type="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P")
	for (x in rownames(xys)) {
	for (y in names(xys)) {
		ivw_p <- xys[x,y]
		if (ivw_p <= 5e-8) {
			dat_X <- read.table(paste0(dir_mr, x, "/", x, ".txt"), header=T) %>% format_data(type ="exposure", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_M4x <- dat_M0 %>% merge(subset(dat_X, select="SNP")) %>% format_data(type ="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_XM <- harmonise_data(dat_X, dat_M4x, action=2) 
			dat_Y0 <- read.table(paste0(dir_gwas, y, ".gz"), header=T) 
			dat_Y4x <- dat_Y0 %>% merge(subset(dat_X, select="SNP")) %>% format_data(type="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P")
			dat_Y4m <- dat_Y0 %>% merge(subset(dat_M4y, select="SNP")) %>% format_data(type="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
			dat_MY <- harmonise_data(dat_M4y, dat_Y4m, action=2) 
			#dat_XY <- harmonise_data(dat_X, dat_Y4x) 
			res_X2M <- TwoSampleMR::mr(dat_XM) %>% filter(method=="Inverse variance weighted") # 1st step 三角形的上坡
				beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- res_X2M %>% pull(pval)
			res_M2Y <- TwoSampleMR::mr(dat_MY) %>% filter(method=="Inverse variance weighted") # 1st step 三角形的下坡
				beta_M2Y <- res_M2Y %>% pull(b); se_M2Y <- res_M2Y %>% pull(se); p_M2Y <- res_M2Y %>% pull(pval)
			beta_2step <- round(beta_X2M * beta_M2Y, 4); beta_2step
			se_2step = round((se_X2Y^2 + se_M2Y^2)^0.5, 4)
			p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F), 3); p_2step
			res_str <- paste(x, m, y, nrow(dat_XM), nrow(dat_MY), beta_2step, se_2step, p_2step, sep='|')
			write.table(res_str, paste0('jie.out'), append=T, quote=F, row.names=F, col.names=F)
		}
	}
	}
}