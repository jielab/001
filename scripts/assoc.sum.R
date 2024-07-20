pacman::p_load(foreach, doParallel, data.table, tidyverse, stringi, TwoSampleMR, MendelianRandomization, mrMed)

pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
pval <- function(b,se) (2*pnorm(-abs(b/se)))
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

dir='/work/sph-huangj'
dir.X = paste0(dir,'/data/gwas/?/clean')
dir.Y = paste0(dir,'/data/gwas/?/clean')
dir.M = paste0(dir,'/data/gwas/?/clean')
Xs = '?'
Ys = '?' 
Ms = gsub('.gz', '', list.files(path=dir.M, pattern='.gz$')); Ms
log_file="?.log"
n_cores=40


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 先定义function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mediation <- function(X, dat.X.raw, dat.X.iv, M, Y, dat.Y.raw) {
	# Harmonize dat.X + dat.M + dat.Y
	dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T)
	dat.M.4x <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
	names(dat.M.raw) <- stri_replace_all_regex(toupper(names(dat.M.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	if (file.exists(paste0(dir.M, '/', M, '.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.top.snp'), header=T); names(dat.M.iv) <- "SNP" 
	} else if (file.exists(paste0(dir.M, '/', M, '.NEW.top.snp'))) {
		dat.M.iv <- read.table(paste0(dir.M, '/', M, '.NEW.top.snp'), header=T)
	} else {
		dat.M.sig <- dat.M.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.M.iv <- dat.M.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
		write.table(dat.M.iv, paste0(dir.M, '/', M, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
	}
	dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub("outcome", "mediator", names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by="SNP")
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	if (bad_row !=0) { write(paste("ERR:", X, Y, M, "X-M-Y have inconsistent alleles !!"), file=log_file, append=TRUE); break }
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
	if (nrow(dat)==0) { write(paste("SKIP:", X, Y, M, "X-M-Y harmonized data is empty !!"), file=log_file, append=TRUE); break }
	# mrMed 方法 🏮🏮
	dat1 <- dat
	names(dat1) <- stri_replace_all_regex(names(dat1), pattern=c("exposure", "mediator", "outcome"), replacement=c("X", "M", "Y"), vectorize_all=FALSE)
	dat1 <- dat1 %>% mutate (
		Gx = ifelse(dat$SNP %in% dat.X.iv$SNP, 1, 0), Gx_plum = Gx, # 可加上 “& pval.M >5e-08”, 用于 X -> M或Y
		Gm = 1, Gm_plum = Gm, # 可考虑 ifelse(dat$SNP %in% dat.M.iv$SNP, 1, 0), 用于 M -> Y 
		G_mvmr = ifelse(Gx_plum==0 & Gm_plum==0, 0, 1)
	)
	if ( max(dat1$Gx) * max(dat1$Gm) ==0 ) {
		mrMed_str <- paste("mrMed has no SNP !!")
	} else {
		res <- mrMed(dat_mrMed=dat1, method_list="Prod_IVW_0")
	}
	# key results # 🔦
	X_str=paste0(X, "(", nrow(dat.X.iv), ")"); M_str=paste0(M, "(", nrow(dat.M.iv), ")");
	write(paste("RES:", X_str, M_str, Y, "|X2Y", nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		"|mrMed", paste(rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p)) # TE, ACME, ADE, Prop(rho)
		), file=log_file, append=TRUE
	)
}
		
			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 基于MR的mediation
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

	for (X in Xs) { # 🍷
		dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
		names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
		if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
		} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
		} else {
			dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
			dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
			write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
		}
		if(nrow(dat.X.iv) ==0) {write(paste("SKIP:", X,"NA",Y, "X no IV !!"), file=log_file, append=TRUE); next}
		
		dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		

		# Harmonize dat.X + dat.Y, 并计算 Total effect
		dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
		fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2Y 
		beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- fit.X2Y$pval
		if(p.X2Y >0.05) { write(paste("SKIP:", X,"NA",Y, nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), "X2Y not significant !!"), file=log_file, append=TRUE); next } # 🛑

#		numCores <- detectCores()-1; cl <- makeCluster(numCores); registerDoParallel(cl)
#		foreach::foreach(M = Ms) %dopar% { # 🐎
		for (M in Ms) {
			mediation(X, dat.X.raw, dat.X.iv, M, Y, dat.Y.raw)
		}
#		stopCluster(cl)
	}
}
