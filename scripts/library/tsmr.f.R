tsmr_fn <-function( iv_X,
	pheno_X, ieu_X, file_X, 
	pheno_Y, ieu_Y, file_Y, 	
	X_snp="SNP", X_ea="EA", X_nea="NEA", X_beta="BETA", X_se="SE", X_p="P",
	Y_snp="SNP", Y_ea="EA", Y_nea="NEA", Y_beta="BETA", Y_se="SE", Y_p="P"
	) {			
	pacman::p_load(TwoSampleMR, MendelianRandomization, MRInstruments, dplyr, ggcorrplot)
	if (file.exists(file_X)) {
		X_dat0 <- read.table(file_X, header=T, as.is=T)
		if (file.exists(iv_file_X)) {
			X_iv <- read.table(iv_file_X, header=T, as.is=T)
		} else {
			X_iv <- subset(X_dat0, select=X_snp)
		}
		X_dat <- merge(X_dat0, X_iv, by.x=X_snp, by.y=names(X_iv)) 
		X_dat <- format_data(X_dat0, type ="exposure", phenotype_col=pheno_X, snp_col=X_snp, effect_allele_col=X_ea, other_allele_col=X_nea, beta_col=X_beta, se_col=X_se, pval_col=X_p) # %>% filter (pval.exposure <1e-3)
	} else if (ieu_X != 'NA') {
		X_dat <- extract_instruments(ieu_X, p1=5e-08, clump=T, r2=0.1, kb=1000) 
		X_iv <- subset(X_dat, select=X_snp)
	} else {
		next
	}
	if (file.exists(file_Y)) {
		Y_dat0 <- read.table(file_Y, header=T, as.is=T)
		Y_dat <- merge(Y_dat0, X_iv, by.x=Y_snp, by.y=names(X_iv))
		if (nrow(Y_dat) >0) {
			Y_dat <- format_data(Y_dat, type="outcome", phenotype_col=pheno_Y, snp_col=Y_snp, effect_allele_col=Y_ea, other_allele_col=Y_nea, beta_col=Y_beta, se_col=Y_se, pval_col=Y_p) # %>% filter (pval.exposure <1e-3)
		}
	} else if (ieu_Y != 'NA') {
		Y_dat <- extract_outcome_data(snps=X_iv[[X_snp]], outcomes=ieu_Y)	
	} else {	
		next
	}
	if (nrow(Y_dat) >0) {
		tsmr_dat <- harmonise_data(X_dat, Y_dat, action=1) # accept palindromic SNPs
		if (nrow(tsmr_dat) >0) {
		res_mr <- mr(tsmr_dat) 
		res_mr.str <- paste(pheno_X, nrow(X_dat), pheno_Y, nrow(Y_dat), res_mr$nsnp, res_mr$method, res_mr$b, res_mr$se, res_mr$pval, sep='|')
		write.table(res_mr.str, paste0(pheno_X,'.',pheno_Y,'.tsmr.out'), append=F, quote=F, row.names=F, col.names=F)
	#	#下面代码用 MendelianRandomizaiton 包画图
	#	bsmr_dat <- dat_to_MRInput(tsmr_dat); bsmr_dat <- bsmr_dat[[1]]
	#	mr_ivw(bsmr_dat)
	#	png(paste0(pheno_X,'.',pheno_Y,'.beta.png'), w=800, h=800)
	#		plt <- mr_plot(bsmr_dat, interactive=F) + theme(axis.title=element_text(size=15, face="bold"), axis.text=element_text(size=12, face="bold"))
	#		print(plt); dev.off()
	#	png(paste0(pheno_X,'.',pheno_Y,'.forest.png'), w=800, h=800)
	#		plt <- mr_forest(bsmr_dat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))
	#		print(plt); dev.off()
	#	png(paste0(pheno_X,'.',pheno_Y,'.funnel.png'), w=800, h=800)	
	#		plt <- mr_funnel(bsmr_dat) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=14), axis.text=element_text(size=12, face="bold"))	
	#		print(plt); dev.off()
		}
	}
}
