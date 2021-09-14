mrfn <-function( indir, ex_name, ex_ieu, ex_gcat, ou_name, ou_ieu, ou_gcat, top_snps_file ) {
	
	pacman::p_load(TwoSampleMR, MRInstruments, dplyr, tidyverse, magrittr)
		
	ex_gcat_file <- paste0(indir,'/',ex_gcat,'.gz')
	ou_gcat_file <- paste0(indir,'/',ou_gcat,'.gz')

	if (ex_gcat != 'NA' & file.exists(top_snps_file)) {
		col_names_ex <- as.character(read.table(ex_gcat_file, header=F, nrow=1))
		SNP_ex <- grep('^SNP$|^variant_id$', col_names_ex, ignore.case=T, value=T)
		CHR_ex <- grep('^CHR$|^chromosome$', col_names_ex, ignore.case=T, value=T)
		POS_ex <- grep('^POS$|^bp$|^base_pair_location$', col_names_ex, ignore.case=T, value=T)
		EA_ex <- grep('^EA$|^effect_allele$|^allele1$|^a1$', col_names_ex, ignore.case=T, value=T)
		NEA_ex <- grep('^NEA$|^other_allele$|^allele0$|^allele2$|^a2$', col_names_ex, ignore.case=T, value=T)
	#	EAF_ex <- grep('^EAF$|^freq$|^effect_allele_frequency$|^a1freq$', col_names_ex, ignore.case=T, value=T)
		BETA_ex <- grep('^BETA$|^effect$', col_names_ex, ignore.case=T, value=T)
		SE_ex <- grep('^SE$|^standard_error$', col_names_ex, ignore.case=T, value=T)
		P_ex <- grep('^P$|^p_value$|^pvalP_LINREG$', col_names_ex, ignore.case=T, value=T)
		ex_dat <- read_exposure_data(filename=ex_gcat_file, sep='\t', snp_col=SNP_ex, effect_allele_col=EA_ex, other_allele_col=NEA_ex, beta_col=BETA_ex, se_col=SE_ex, pval_col=P_ex) # %>% filter (pval.exposure <1e-3)
		top_snps <- read.table(top_snps_file, header=F); names(top_snps) <- 'SNP'
		ex_dat <- merge(ex_dat, top_snps, by='SNP') 
	} else if (ex_ieu != 'NA') {
		ex_dat <- extract_instruments(ex_ieu, p1=5e-08, clump=T, r2=0.1, kb=1000) 
	} else {	
		next
	}
	if (ou_gcat != 'NA') { # for outcome, use gwas-catalog data first
		col_names_ou <- as.character(read.table(ou_gcat_file, header=F, nrow=1))
		SNP_ou <- grep('^SNP$|^variant_id$', col_names_ou, ignore.case=T, value=T)
		CHR_ou <- grep('^CHR$|^chromosome$', col_names_ou, ignore.case=T, value=T)
		POS_ou <- grep('^POS$|^bp$|^base_pair_location$', col_names_ou, ignore.case=T, value=T)
		EA_ou <- grep('^EA$|^effect_allele$|^allele1$|^a1$', col_names_ou, ignore.case=T, value=T)
		NEA_ou <- grep('^NEA$|^other_allele$|^allele0$|^allele2$|^a2$', col_names_ou, ignore.case=T, value=T)
	#	EAF_ou <- grep('^EAF$|^freq$|^effect_allele_frequency$|^a1freq$', col_names_ou, ignore.case=T, value=T)
		BETA_ou <- grep('^BETA$|^b$|^effect$', col_names_ou, ignore.case=T, value=T)
		SE_ou <- grep('^SE$|^standard_error$', col_names_ou, ignore.case=T, value=T)
		P_ou <- grep('^P$|^p_value$|^pvalP_LINREG$', col_names_ou, ignore.case=T, value=T)
		ou_dat <- read_outcome_data(snps=ex_dat$SNP, filename=ou_gcat_file, sep='\t', snp_col=SNP_ou, effect_allele_col=EA_ou, other_allele_col=NEA_ou, beta_col=BETA_ou, se_col=SE_ou, pval_col=P_ou) 
	} else if (ou_ieu != 'NA') {
		ou_dat <- extract_outcome_data(snps=ex_dat$SNP, outcomes=ou_ieu)	
	} else {	
		next
	}
	dat <- harmonise_data(ex_dat, ou_dat)
	write.table(dat, paste0(ex_name,'.',ou_name,'.dat.txt'), append=F, quote=F, row.names=F, col.names=T)
	res <- mr(dat, method_list=c('mr_ivw', 'mr_simple_median', 'mr_weighted_median', 'mr_two_sample_ml', 'mr_penalised_weighted_median', 'mr_egger_regression')) 
	res_str <- paste(ex_name, nrow(ex_dat), ou_name, nrow(ou_dat), res$method, res$pval, sep='|')
	write.table(res_str, paste0(ex_name,'.',ou_name,'.res.txt'), append=F, quote=F, row.names=F, col.names=F)
}


	# pl_res <- mr_scatter_plot(res, dat) # 
	# mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))
	# res <- mr_pleiotropy_test(dat); pl_mr <- mr_scatter_plot(res, dat); ggsave(p1[[1]], file='filename.png', width=7, height=7)
	# res_single <- mr_singlesnp(dat, all_method=c("mr_ivw", "mr_two_sample_ml")); pl_forest <- mr_forest_plot(res_single); ggsave(p_forest[[1]], file="forest.png", width=7, height=7)
	# res_loo <- mr_leaveoneout(dat); pl_loo <- mr_leaveoneout_plot(res_loo); ggsave(p_loo[[1]], file="loo.png", width=7, height=7)
	# p_funnel <- mr_funnel_plot(res_single); ggsave(p_forest[[1]], file="forest.png", width=7, height=7)
	
	next
	
	#### compareB
	source('$dir/scripts/library/compareB.f.R')
	pdf('$ex_name.$ou_name.compareB.pdf')
	par(mfrow=c(3,2), mai=c(0.5,1,0.5,0.5))
	compareB(
		f1='$ex_name.top.txt', f1_name='$ex_name', f1_snp='SNP', f1_ea='A1', f1_nea='A2', f1_eaf='freq', f1_beta='b', f1_se='se', f1_p='p',
		f2='$ou_name.for.$ex_name.top.txt', f2_name='$ou_name', f2_snp='SNP', f2_ea='A1', f2_nea='A2', f2_eaf='freq', f2_beta='b', f2_se='se', f2_p='p'
	)
	dev.off()
	
	### MendelianRandomization 
	library('MendelianRandomization')
	# read regular file
	dat = read.table('D:/data/gwas/posCtrls/cad.w.t2d.txt', header=T)
	dat$t2d.BETA[dat$BETA<0] = 0-dat$t2d.BETA[dat$BETA<0]; dat$BETA[dat$BETA<0] = 0-dat$BETA[dat$BETA<0] 
	XGb = dat$BETA; XGse = dat$SE; YGb = dat$t2d.BETA; YGse = dat$t2d.SE
	plot(XGb, YGb); abline(h=0); abline(v=0, col='blue')
	# read compareB merged file
	dat = read.table('asthma.hf.merged.txt', header=T); head(dat)
	XGb = dat$BETA.aligned_1; XGse = dat$SE_1; YGb = dat$BETA.aligned_2; YGse = dat$SE_2
	# read GCTA output
	source('D:/scripts/library/gsmr_plot.r')
	gsmr.data = read_gsmr_data('D:/analysis/gsmr/crp.copd.eff_plot.gz'); str(gsmr.data)
	dat = gsmr_snp_effect(gsmr.data, 'crp', 'copd'); str(dat)
	dat$bzy[dat$bzx<0] = 0-dat$bzy[dat$bzx<0]; dat$bzx[dat$bzx<0] = 0-dat$bzx[dat$bzx<0] 
	XGb = dat$bzx; XGse = dat$bzx_se; YGb = dat$bzy; YGse = dat$bzy_se
	# MR and plot
	png('$ex_name.$ou_name.mr.png', w=800, res=128)
	dat = read.table('$ex_name.$ou_name.merged.txt', header=T)
	XGb = dat$BETA.aligned_1; XGse = dat$SE_1; YGb = dat$BETA.aligned_2; YGse = dat$SE_2
	mr_forest(mr_input(XGb, XGse, YGb, YGse), snp_estimates=F, methods=c('ivw', 'median', 'wmedian', 'egger', 'maxlik', 'mbe', 'conmix'))
	dev.off()
}