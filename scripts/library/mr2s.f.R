mr2sfn <-function( top_snps_file, 
	ex_name, ex_ieu, ex_gwas, ex_snp, ex_ea, ex_nea, ex_eaf, ex_beta, ex_se, ex_p,
	ou_name, ou_ieu, ou_gwas, ou_snp, ou_ea, ou_nea, ou_eaf, ou_beta, ou_se, ou_p 
	) {
			
	pacman::p_load(TwoSampleMR, MRInstruments, dplyr, tidyverse)

	if (file.exists(ex_gwas_file) & file.exists(top_snps_file)) {
		ex_col_names <- as.character(read.table(ex_gwas_file, header=F, nrow=1))
		for (var in c('ex_snp', 'ex_ea', 'ex_nea', 'ex_eaf', 'ex_beta', 'ex_se', 'ex_p')) {eval(parse(text= paste0(var, '=grep(', var, ', ex_col_names, ignore.case=T, value=T)')))}
		ex_dat <- read_exposure_data(filename=ex_gwas_file, sep='\t', snp_col=ex_snp, effect_allele_col=ex_ea, other_allele_col=ex_nea, beta_col=ex_beta, se_col=ex_se, pval_col=ex_p) # %>% filter (pval.exposure <1e-3)
		top_snps <- read.table(top_snps_file, header=F); names(top_snps) <- 'SNP'
		ex_dat <- merge(ex_dat, top_snps, by='SNP') 
	} else if (ex_ieu != 'NA') {
		ex_dat <- extract_instruments(ex_ieu, p1=5e-08, clump=T, r2=0.1, kb=1000) 
	} else {	
		next
	}
	if (file.exists(ou_gwas_file)) { # for outcome, use gwas-catalog data first
		ou_col_names <- as.character(read.table(ou_gwas_file, header=F, nrow=1))
		for (var in c('ou_snp', 'ou_ea', 'ou_nea', 'ou_eaf', 'ou_beta', 'ou_se', 'ou_p')) {eval(parse(text= paste0(var, '=grep(', var, ', ou_col_names, ignore.case=T, value=T)')))}
		ou_dat <- read_outcome_data(snps=ex_dat$SNP, filename=ou_gwas_file, sep='\t', snp_col=SNP_ou, effect_allele_col=EA_ou, other_allele_col=NEA_ou, beta_col=BETA_ou, se_col=SE_ou, pval_col=P_ou) 
	} else if (ou_ieu != 'NA') {
		ou_dat <- extract_outcome_data(snps=ex_dat$SNP, outcomes=ou_ieu)	
	} else {	
		next
	}
	dat <- harmonise_data(ex_dat, ou_dat)
	write.table(dat, paste0(ex_name,'.',ou_name,'.mr2s.dat.txt'), sep='\t', append=F, quote=F, row.names=F, col.names=T)
	mr.res <- mr(dat, method_list=c('mr_ivw', 'mr_simple_median', 'mr_weighted_median', 'mr_two_sample_ml', 'mr_penalised_weighted_median', 'mr_egger_regression')) 
	mr.res.str <- paste(ex_name, nrow(ex_dat), ou_name, nrow(ou_dat), mr.res$method, mr.res$pval, sep='|')
	write.table(mr.res.str, paste0(ex_name,'.',ou_name,'.mr2s.res.txt'), append=F, quote=F, row.names=F, col.names=F)

}
