pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')

IV.filter <- function(dat, CHR_col, POS_col, AF_col, BETA_col, N_col) {
	filter1 <- dat[[CHR_col]]==6 & dat[[POS_col]] >=28510120 & dat[[POS_col]] <=33480577 
	dat$R2 <- 2 * dat[[AF_col]] * (1-dat[[AF_col]]) * (dat[[BETA_col]]^2)
	dat$F_stat <- dat$R2 * (dat[[N_col]]-2) / (1-dat$R2)
	filter2 <- dat$F_sta < 10
	return(dat[!(filter1 | filter2), ])
}


run_mr2s <- function( # 🏮 常规的 TwoSampleMR
  analysis, p_t, label, X, dir.X, dat.X.raw, IV.file, IV_filter, Y, dir.Y, dat.Y.raw) {
  
	pval_t <- as.numeric(paste0('5e-', p_t))
	log_file <- paste0(label, '.mr.log')
	if (IV.file != 'NO') {
		dat.X.iv <- read.table(IV.file, header=T); names(dat.X.iv) <- 'SNP'
	} else if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- 'SNP' 
	} else if (file.exists(paste0(dir.X, '/', X, '.NEW.', p_t, '.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.', p_t, '.top.snp'), header=T)
	} else {
		dat.X.iv <- dat.X.raw %>% filter(P<=pval_t) %>% mutate(mb=ceiling(POS/1e+05))
				dat.X.iv <- dat.X.iv %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select('SNP')
		write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.', p_t, '.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
	}

	if (nrow(dat.X.iv)==0) {write(paste('SKIP:', analysis, X,Y, 'X has no IV !!'), file=log_file, append=TRUE); return(NULL)}	
	dat.X <- dat.X.raw %>% merge(dat.X.iv); if (!'N' %in% names(dat.X)) {dat.X$N=100000}
		if (IV_filter==1) { dat.X <- IV.filter(dat.X, 'CHR', 'POS', 'EAF', 'BETA', 'N') }
		if (nrow(dat.X)==0) {write(paste('SKIP:', analysis, X,Y, 'X is empty after IV filter!!'), file=log_file, append=TRUE); return(NULL)}
		dat.X <- dat.X %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.Y.4x <- subset(dat.Y.raw, SNP %in% dat.X$SNP); if (!'N' %in% names(dat.Y.4x)) {dat.Y.4x$N=100000} 
		if(nrow(dat.Y.4x)==0) {write(paste('SKIP:', analysis, X,Y, 'IV SNPs do NOT exist in Y!!'), file=log_file, append=TRUE); return(NULL)}	
		dat.Y.4x <- dat.Y.4x %>% format_data(type='outcome', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.XY <- harmonise_data(dat.X, dat.Y.4x, action=1) 
	write.table(dat.XY, paste0(label, '.dat'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)


	if (nrow(dat.XY) >1) {
		fit.X2Y <- mr(dat.XY, method_list=c('mr_ivw', 'mr_two_sample_ml', 'mr_egger_regression', 'mr_weighted_median', 'mr_ivw_radial'))
		fit.pres <- try( {run_mr_presso(dat.XY, NbDistribution=1000, SignifThreshold=0.05)}, silent=TRUE )
		fit.pres.p <- try( fit.pres[[1]]$`Main MR results`$`P-value`, silent=TRUE )
	} else if (nrow(dat.XY) ==1) {
		fit.X2Y <- mr(dat.XY, method_list=c('mr_wald_ratio'))
	} else {
		return(NULL)
	}
	beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- fit.X2Y$pval
	p.hetero <- mr_heterogeneity(dat.XY)$Q_pval
	p.pleio <- mr_pleiotropy_test(dat.XY)$pval
	method <- str_replace_all(fit.X2Y$method, c("Inverse variance weighted"="ivw", "Maximum likelihood"="ml", "MR Egger"="egger", "Weighted median"="median", "IVW radial"="radial"))
	write(cat(analysis, p_t, X,Y, nrow(dat.XY), '|METHOD|', method, '|BETA|', rb(beta.X2Y), '|SE|', rb(se.X2Y), '|P|', rp(p.X2Y), '|PRESSO|', fit.pres.p, '|HETERO|', rp(p.hetero), '|PLEIO|', rp(p.pleio), '|P.Y.min|', rp(min(dat.XY$pval.outcome))), file=log_file, append=TRUE)
	
	return(min(p.X2Y))
}


run_cisMr <- function( # 🏮 技术讨论见 https://github.com/ZhaotongL/cisMRcML/issues/6
  analysis, p_t, label, chr, flank, pos0, pos1, X, dat.X.raw, dir.X.cojo, Y, dat.Y.raw, run_coloc) {

	pval_t <- as.numeric(paste0('5e-', p_t))
	log_file <- paste0(label, '.cisMr.log')
	pos0=as.numeric(pos0); pos1=as.numeric(pos1)
	bfile <- paste('--bfile', paste0(dir.X.cojo, '/', X)) 
	
	dat.X <- dat.X.raw %>% filter(CHR==chr, POS>=(pos0-flank), POS<=(pos1+flank)); if (!'N' %in% names(dat.X)) {dat.X$N=100000}
	dat.Y <- dat.Y.raw %>% filter(SNP %in% dat.X$SNP) %>% mutate(N=100000); if (!'N' %in% names(dat.Y)) {dat.X$N=100000}
	dat.X <- dat.X %>% format_data(type='exposure', snp_col='SNP', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.Y <- dat.Y %>% format_data(type='outcome',  snp_col='SNP', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.X[is.na(dat.X)] <- 0; dat.Y[is.na(dat.Y)] <- 0 # 上面的命令会将0变成NA，所以调回来🏮
	dat.XY <- harmonise_data(dat.X, dat.Y, action=1)
		if (nrow(dat.XY)==0) {write(paste('SKIP:', analysis, X, Y, 'merged data is empty'), file=log_file, append=TRUE); return(NULL)}
		if (min(dat.XY$pval.exposure) > pval_t) { write(paste('SKIP:', analysis, X, Y, 'X .4gcta file has no significant SNPs'), file=log_file, append=TRUE); return(NULL) }
		if (all(is.na(dat.XY$eaf.exposure)) & all(is.na(dat.XY$eaf.outcome))) {write(paste('ERROR:', analysis, X, Y, 'neither data has EAF'), file=log_file, append=TRUE); return(NULL)} 
		dat.XY <- dat.XY %>% mutate(eaf.exposure = ifelse(is.na(eaf.exposure), eaf.outcome, eaf.exposure), eaf.outcome = ifelse(is.na(eaf.outcome), eaf.exposure, eaf.outcome))
		dat.X.coloc <- dat.XY %>% {list(type='quant', snp=.$SNP, position=.$pos.exposure, MAF=.$eaf.exposure, N=.$samplesize.exposure, beta=.$beta.exposure, varbeta =.$se.exposure^2, pvalues=.$pval.exposure)}
		dat.Y.coloc <- dat.XY %>% {list(type='cc', snp=.$SNP, position=.$pos.outcome, MAF=.$eaf.outcome, N=.$samplesize.outcome, beta=.$beta.outcome, varbeta =.$se.outcome^2, pvalues=.$pval.outcome)}
		dat.X <- dat.XY %>% dplyr::select(SNP, effect_allele.exposure, other_allele.exposure, eaf.exposure, beta.exposure, se.exposure, pval.exposure, samplesize.exposure) 
		dat.Y <- dat.XY %>% dplyr::select(SNP, effect_allele.outcome, other_allele.outcome, eaf.outcome, beta.outcome, se.outcome, pval.outcome, samplesize.outcome) 
		names(dat.X) <- c('SNP', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'N'); fwrite(dat.X, paste0(label,'.X.4gcta'), sep='\t', na=0, quote=FALSE)
		names(dat.Y) <- c('SNP', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'N'); fwrite(dat.Y, paste0(label,'.Y.4gcta'), sep='\t', na=0, quote=FALSE)
	
	X.prefix <- paste0(X,'.',p_t)				
	X.cojo.cmd <- paste0('gcta ', bfile, ' --cojo-file ', label,'.X.4gcta --cojo-slct --cojo-p ', pval_t, ' --out ', X.prefix); ecode = system(X.cojo.cmd, intern=FALSE)
		if(ecode!=0) { write(paste('ERROR:', analysis, X, Y, '--cojo-slct on X does not run'), file=log_file, append=TRUE); return(NULL) }
		X.cojo.fn <- paste0(X.prefix,'.jma.cojo')
		X.cojo <- fread(X.cojo.fn); if(length(X.cojo$SNP)<3) {write(paste('SKIP:', analysis, X, Y, 'jma.cojo file has less than 3 records'), file=log_file, append=TRUE); return(NULL)}
	
	Y.prefix <- paste0(Y,'.',X,'.',p_t)
	Y.cojo.cmd <- paste0('gcta ', bfile, ' --cojo-file ', label,'.Y.4gcta --cojo-slct --cojo-p ', pval_t, ' --out ', Y.prefix); system(Y.cojo.cmd, intern=FALSE)
		Y.cojo.fn <- paste0(Y.prefix,'.jma.cojo')
		if (file.exists(Y.cojo.fn) && file.size(Y.cojo.fn) >0) {
			Y.cojo = fread(Y.cojo.fn); IV = union(X.cojo$SNP, Y.cojo$SNP); IV.fn = paste0(label,'.',p_t,'.iv')
			write.table(IV, IV.fn, append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
			XY.cojo.cmd = paste0('gcta ', bfile, ' --cojo-file ', label, '.X.4gcta --extract ', IV.fn, ' --cojo-slct --cojo-p ', pval_t, ' --out ', X.prefix); ecode=system(XY.cojo.cmd, intern=FALSE)
			if (ecode==1) { write(paste('ERROR:', analysis, X, Y, '--extract --cojo-joint does not run'), file=log_file, append=TRUE); return(NULL) }
		} else {
			IV=X.cojo$SNP
		}
			
	LD_mat.fn = paste0(X.prefix, '.ldr.cojo')
	LD_mat_cmd = paste("sed -i 's/\t$//'", LD_mat.fn); system(LD_mat_cmd)
	LD_mat = fread(LD_mat.fn)
	LD_mat = LD_mat[,2:ncol(LD_mat)]; LD_mat = as.matrix(LD_mat); rownames(LD_mat) = colnames(LD_mat)
	
	dat.X <- dat.X %>% filter(SNP %in% colnames(LD_mat))
		dat.X = dat.X[match(colnames(LD_mat), dat.X$SNP),] # 必须的🏮
		dat.X$cor = dat.X$BETA / sqrt(dat.X$BETA^2 + (dat.X$N-2) * dat.X$SE^2)
		dat.X$se_cor = dat.X$SE * dat.X$cor/dat.X$BETA
		dat.X$bJ = solve(LD_mat) %*% dat.X$cor
	
	dat.Y = dat.Y %>% filter(SNP %in% colnames(LD_mat))
		dat.Y = dat.Y[match(colnames(LD_mat), dat.Y$SNP),] # 必须的🏮
		dat.Y$cor = dat.Y$BETA / sqrt(dat.Y$BETA^2 + (dat.Y$N-2) * dat.Y$SE^2)
		dat.Y$se_cor = dat.Y$SE * dat.Y$cor/dat.Y$BETA
		dat.Y$bJ = solve(LD_mat) %*% dat.Y$cor
	
	dat.mr = list(b_exp=dat.X$bJ, b_out=dat.Y$bJ, N1=median(dat.X$N), N2=median(dat.Y$N), LD_mat=LD_mat, exp_df=dat.X, out_df=dat.Y, exp_IV=X.cojo$SNP, out_IV=setdiff(IV, X.cojo$SNP))
	
	Sig_exp1 = solve(dat.mr$LD_mat) %*% (dat.mr$exp_df$se_cor %o% dat.mr$exp_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat); Sig_exp_inv=solve(Sig_exp1)
	Sig_out1 = solve(dat.mr$LD_mat) %*% (dat.mr$out_df$se_cor %o% dat.mr$out_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat); Sig_out_inv=solve(Sig_out1)
	res = cismr_cML_DP(
		b_exp=dat.mr$exp_df$bJ, b_out=dat.mr$out_df$bJ, Sig_exp_inv=Sig_exp_inv, Sig_out_inv=Sig_out_inv, maxit=200, n=dat.mr$N1, 
		random_start=5, min_theta_range=-0.1, max_theta_range=0.1, num_pert=100, random_start_pert=5, random_seed=12345
	)
	write(paste(analysis, p_t, X, Y, length(IV), rb(res$BIC_DP_theta), rb(res$BIC_DP_se), rp(res$BIC_DP_p)), file=log_file, append=TRUE)

	if (run_coloc != 1) {return(NULL)}
	log_file <- paste0(label, '.coloc.log')
	if (min(dat.XY$pval.outcome) >1e-2) { write(paste('SKIP:', analysis, X, Y, 'Y minimum P', min(dat.XY$pval.outcome), 'is not significant'), file=log_file, append=TRUE); return(NULL)}
	fit.coloc <- coloc::coloc.abf(dat.X.coloc, dat.Y.coloc)
	res <- paste(signif(fit.coloc$summary, 3),collapse=' ') # nsnps H0 H1 H2 H3 H4 
	write(paste(X, Y, res), file=log_file, append=TRUE)
	if (fit.coloc$summary[6] >0.5) { pdf(paste0(label,'.coloc.pdf')); sensitivity(fit.coloc,'H4>0.5'); dev.off() }
}


run_mrMed <- function( # 🏮 devtools::install_github('scllin/mrMed') 
  label, X, dat.X.raw, file.M, Y, dat.Y.raw) { 

	log_file <- paste0(label, '.mrMed.log')
	dir.M <- dirname(file.M)
	M <- sub('\\.gz$', '', basename(file.M))
	
	if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- 'SNP' 
	} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
	} else {
		dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select('SNP')
		write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
	}
	
	dat.M.raw <- read.table(file.M, header=T)
	names(dat.M.raw) <- stri_replace_all_regex(toupper(names(dat.M.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- 'SNP' 
	} else if (file.exists(paste0(dir.M, '/', M, '.NEW.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.NEW.top.snp'), header=T)
	} else {
		dat.M.sig <- dat.M.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.M.iv <- dat.M.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select('SNP')
		write.table(dat.M.iv, paste0(dir.M, '/', M, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
	}
	
	dat.XnM.iv <- rbind(IV, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) 
	if (nrow(dat.M.mv)==0) {
		write(paste('SKIP:', X, Y, M, 'no SNP in dat.M.mv'), file=log_file, append=TRUE); return()
	} else {
		dat.M.mv <- dat.M.mv %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	}
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub('outcome', 'mediator', names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by='SNP')
	if (nrow(dat)==0) { write(paste('SKIP:', X, Y, M, 'X-M-Y harmonized data is empty'), file=log_file, append=TRUE); return()}
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	if (bad_row !=0) { write(paste(X, Y, M, 'ERR: X-M-Y have inconsistent alleles'), file=log_file, append=TRUE); return()}
	names(dat) <- gsub('\\.x$', '', names(dat)); dat <- subset(dat, select=!grepl('\\.y', names(dat)))
	
	dat1 <- dat # mrMed 方法 
	names(dat1) <- stri_replace_all_regex(names(dat1), pattern=c('exposure', 'mediator', 'outcome'), replacement=c('X', 'M', 'Y'), vectorize_all=FALSE)
	dat1 <- dat1 %>% mutate (
		Gx = ifelse(dat$SNP %in% IV$SNP, 1, 0), Gx_plum = Gx, # 可加上 “& pval.M >5e-08”, 用于 X -> M或Y
		Gm = 1, Gm_plum = Gm, # 可考虑 ifelse(dat$SNP %in% dat.M.iv$SNP, 1, 0), 用于 M -> Y 
		G_mvmr = ifelse(Gx_plum==0 & Gm_plum==0, 0, 1)
	)
	
	if ( max(dat1$Gx) * max(dat1$Gm)==0 ) { write(paste(X, Y, M, 'mrMed data has 0 rows for Gx or Gm'), file=log_file, append=TRUE); return() }
	res <- try( {mrMed(dat_mrMed=dat1, method_list='Prod_IVW_0')}, silent=TRUE )	
	if (inherits(res, 'try-error')) { write(paste('ERROR:', X, Y, M, 'mrMed gives ERROR message'), file=log_file, append=TRUE); return() }
	X_str=paste0(X, '(', nrow(IV), ')'); M_str=paste0(M, '(', nrow(dat.M.iv), ')');
	write(paste(X_str, M_str, Y, nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		paste(rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p))
		), file=log_file, append=TRUE
	)

}
