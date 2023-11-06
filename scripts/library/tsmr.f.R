tsmr_fn <-function( exp_iv_file,
	exp_name, exp_pheno, exp_ieu, exp_file, 
	out_name, out_pheno, out_ieu, out_file, 
	exp_snp="SNP", exp_ea="EA", exp_nea="NEA", exp_beta="BETA", exp_se="SE", exp_p="P",
	out_snp="SNP", out_ea="EA", out_nea="NEA", out_beta="BETA", out_se="SE", out_p="P"
	) {			
	pacman::p_load(TwoSampleMR, MendelianRandomization, MRInstruments, dplyr, ggcorrplot)
	if (file.exists(exp_file)) {
		exp_dat0 <- read.table(exp_file, header=T, as.is=T)
		if (file.exists(exp_iv_file)) {
			exp_iv <- read.table(exp_iv_file, header=T, as.is=T)
		} else {
			exp_iv <- subset(exp_dat0, select=exp_snp)
		}
		exp_dat <- merge(exp_dat0, exp_iv, by.x=exp_snp, by.y=names(exp_iv)) 
		exp_dat <- format_data(exp_dat0, type ="exposure", phenotype_col=exp_pheno, snp_col=exp_snp, effect_allele_col=exp_ea, other_allele_col=exp_nea, beta_col=exp_beta, se_col=exp_se, pval_col=exp_p) # %>% filter (pval.exposure <1e-3)
	} else if (exp_ieu != 'NA') {
		exp_dat <- extract_instruments(exp_ieu, p1=5e-08, clump=T, r2=0.1, kb=1000) 
		exp_iv <- subset(exp_dat, select=exp_snp)
	} else {
		next
	}
	if (file.exists(out_file)) {
		out_dat0 <- read.table(out_file, header=T, as.is=T)
		out_dat <- merge(out_dat0, exp_iv, by.x=out_snp, by.y=names(exp_iv))
		if (nrow(out_dat) >0) {
			out_dat <- format_data(out_dat, type="outcome", phenotype_col=out_pheno, snp_col=out_snp, effect_allele_col=out_ea, other_allele_col=out_nea, beta_col=out_beta, se_col=out_se, pval_col=out_p) # %>% filter (pval.exposure <1e-3)
		}
	} else if (out_ieu != 'NA') {
		out_dat <- extract_outcome_data(snps=exp_iv[[exp_snp]], outcomes=out_ieu)	
	} else {	
		next
	}
	if (nrow(out_dat) >0) {
		tsmr_dat <- harmonise_data(exp_dat, out_dat, action=1) # accept palindromic SNPs
		if (nrow(tsmr_dat) >0) {
		res_mr <- mr(tsmr_dat) 
		res_mr.str <- paste(exp_name, nrow(exp_dat), out_name, nrow(out_dat), res_mr$nsnp, res_mr$method, res_mr$b, res_mr$se, res_mr$pval, sep='|')
		write.table(res_mr.str, paste0(exp_name,'.',out_name,'.tsmr.out'), append=F, quote=F, row.names=F, col.names=F)
	#	#下面代码用 MendelianRandomizaiton 包，画图
	#	bsmr_dat <- dat_to_MRInput(tsmr_dat); bsmr_dat <- bsmr_dat[[1]]
	#	mr_ivw(bsmr_dat)
	#	png(paste0(exp_name,'.',out_name,'.beta.png'), w=800, h=800)
	#		plt <- mr_plot(bsmr_dat, interactive=F) + theme(axis.title=element_text(size=15, face="bold"), axis.text=element_text(size=12, face="bold"))
	#		print(plt); dev.off()
	#	png(paste0(exp_name,'.',out_name,'.forest.png'), w=800, h=800)
	#		plt <- mr_forest(bsmr_dat, snp_estimates=F, methods=c("ivw", "median", "wmedian", "egger", "maxlik", "conmix")) + scale_colour_manual(values = c("IVW estimate"="red")) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=10, face="bold"), axis.text.x=element_text(size=10, face="bold"), axis.text.y=element_text(size=15, face="bold"))
	#		print(plt); dev.off()
	#	png(paste0(exp_name,'.',out_name,'.funnel.png'), w=800, h=800)	
	#		plt <- mr_funnel(bsmr_dat) + theme(axis.title.y=element_text(size=20, face="bold"), axis.title.x=element_text(size=14), axis.text=element_text(size=12, face="bold"))	
	#		print(plt); dev.off()
		if (res_mr$pval < 1e-6) {
			for (m in c('bb_ALB.gz', 'bb_ALB.top.snps', 'bb_ALP.gz', 'bb_ALP.new.snps', 'bb_ALP.new.txt', 'bb_ALP.p03', 'bb_ALP.top.snps', 'bb_ALT.gz', 'bb_ALT.top.snps', 'bb_APOA.gz', 'bb_APOA.top.snps', 'bb_APOB.gz', 'bb_APOB.top.snps', 'bb_AST.gz', 'bb_AST.top.snps', 'bb_AST2ALT.gz', 'bb_AST2ALT.top.snps', 'bb_BILD.gz', 'bb_BILD.top.snps', 'bb_BUN.gz', 'bb_BUN.top.snps', 'bb_CA.gz', 'bb_CA.top.snps', 'bb_CHOL.gz', 'bb_CHOL.top.snps', 'bb_CRE.gz', 'bb_CRE.top.snps', 'bb_CRP.4.snps', 'bb_CRP.gz', 'bb_CRP.top.snps', 'bb_CYS.gz', 'bb_CYS.top.snps', 'bb_EGFR.gz', 'bb_EGFR.top.snps', 'bb_GGT.gz', 'bb_GGT.top.snps', 'bb_GLU.gz', 'bb_GLU.top.snps', 'bb_HBA1C.gz', 'bb_HBA1C.top.snps', 'bb_HDL.gz', 'bb_HDL.top.snps', 'bb_IGF1.gz', 'bb_IGF1.top.snps', 'bb_LDLD.gz', 'bb_LDLD.top.snps', 'bb_LPA.gz', 'bb_LPA.top.snps', 'bb_NAP.gz', 'bb_NAP.top.snps', 'bb_PHOS.gz', 'bb_PHOS.top.snps', 'bb_SHBG.gz', 'bb_SHBG.top.snps', 'bb_TBIL.gz', 'bb_TBIL.top.snps', 'bb_TES.gz', 'bb_TES.top.snps', 'bb_TP.gz', 'bb_TP.top.snps', 'bb_TRIG.gz', 'bb_TRIG.top.snps', 'bb_UA.gz', 'bb_UA.top.snps', 'bb_UCR.gz', 'bb_UCR.top.snps', 'bb_URK.gz', 'bb_URK.top.snps', 'bb_URMA.gz', 'bb_URMA.top.snps', 'bb_URNA.gz', 'bb_URNA.top.snps', 'bb_VITD.gz', 'bb_VITD.top.snps')) {
				## 跑 mediation 分析
			}
		}
	}
	}
}
