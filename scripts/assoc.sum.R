# devtools::install_github("scllin/mrMed") 
pacman::p_load(data.table, tidyverse, stringi, TwoSampleMR, MendelianRandomization, mrMed, cisMRcML) # foreach, doParallel
pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
pval <- function(b,se) (2*pnorm(-abs(b/se)))
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

dir='/work/sph-huangj'
dir.X = paste0(dir,'/data/gwas/?/clean')
dir.Y = paste0(dir,'/data/gwas/?/clean')
dir.M = paste0(dir,'/data/gwas/?/clean')
cis=0
Xs = '?' # 可以仅是一个变量
Ys = '?' # 可以仅是一个变量
label='?'; log_mr=paste0(label,'.mr.log'); log_mrMed=paste0(label,'.mrMed.log')
if (grepl("ALL", Xs)) {Xs = gsub('.gz', '', list.files(path=dir.X, pattern='.gz$'))}
if (grepl("ALL", Ys)) {Ys = gsub('.gz', '', list.files(path=dir.Y, pattern='.gz$'))}
if (!grepl("NA", dir.M)) {Ms = gsub('.gz', '', list.files(path=dir.M, pattern='.gz$'))}
n_cores=40


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 中介function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mediation <- function(X, dat.X.raw, dat.X.iv, M, Y, dat.Y.raw, log_mrMed) {
	# Harmonize dat.X + dat.M + dat.Y
	dat.M.raw <- read.table(paste0(dir.M, '/', M, '.gz'), header=T)
	dat.M.4x <- dat.M.raw %>% merge(dat.X.iv) %>% format_data(type='outcome', snp_col='SNP',  effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome=M)
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
	dat.XnM.iv <- rbind(dat.X.iv, dat.M.iv) %>% unique() # 最理想的是把 dat.X和dat.M中所有的显著性snp合并然后 %>% clump_data() %>% select(SNP)
	dat.X.mv <- dat.X.raw %>% merge(dat.XnM.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
	dat.M.mv <- dat.M.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub('outcome', 'mediator', names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by='SNP')
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	if (bad_row !=0) { write(paste("ERR:", X, Y, M, "X-M-Y have inconsistent alleles !!"), file=log_mrMed, append=TRUE); break }
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
	if (nrow(dat)==0) { write(paste("SKIP:", X, Y, M, "X-M-Y harmonized data is empty !!"), file=log_mrMed, append=TRUE); break }
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
	X_str=paste0(X, "(", nrow(dat.X.iv), ")"); M_str=paste0(M, "(", nrow(dat.M.iv), ")");
	write(paste("RES:", X_str, M_str, Y, "|X2Y", nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		"|mrMed", paste(rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p)) # TE, ACME, ADE, Prop(rho)
		), file=log_mrMed, append=TRUE
	)
}

			
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 进行分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

	for (X in Xs) { # 🍷
		dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
		names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
		if (cis==1) { dat.X.raw <- subset(dat.X.raw, CHR==chr & POS>=(pos_begin-flank) & POS<=(pos_end+flank)) }
		
		# 基于某个locus所有SNP进行cisMR-cML分析 🏮🏮
		if (cisMR==1) {
			mr_dat$exp_df$cor = mr_dat$exp_df$b / sqrt(mr_dat$exp_df$b^2 + (mr_dat$exp_df$N - 2) * mr_dat$exp_df$se^2)
			mr_dat$exp_df$bJ = solve(mr_dat$LD_mat) %*% mr_dat$exp_df$cor
			mr_dat$out_df$cor = mr_dat$out_df$b / sqrt(mr_dat$out_df$b^2 + (mr_dat$out_df$N - 2) * mr_dat$out_df$se^2)
			mr_dat$out_df$bJ = solve(mr_dat$LD_mat) %*% mr_dat$out_df$cor
			b_exp_cond=mr_dat$exp_df$bJ
			b_out_cond=mr_dat$out_df$bJ 
			Sig_exp1 = solve(mr_dat$LD_mat) %*% (mr_dat$exp_df$se_cor %o% mr_dat$exp_df$se_cor * mr_dat$LD_mat) %*% solve(mr_dat$LD_mat)
			Sig_out1 = solve(mr_dat$LD_mat) %*% (mr_dat$out_df$se_cor %o% mr_dat$out_df$se_cor * mr_dat$LD_mat) %*% solve(mr_dat$LD_mat)
			Sig_exp_inv=solve(Sig_exp1)
			Sig_out_inv=solve(Sig_out1)
			res = cismr_cML_DP(
				b_exp=b_exp_cond,b_out=b_out_cond, Sig_exp_inv=Sig_exp_inv,Sig_out_inv=Sig_out_inv,maxit=200, n = mr_dat$N1,random_start = 5,
                min_theta_range=-0.1,max_theta_range=0.1, num_pert=100,random_start_pert=5,random_seed = 12345
			)			  
			c(res$BIC_DP_theta, res$BIC_DP_se, res$BIC_DP_p)
		}

		# 基于显著SNP工具变量进行常规MR分析 🏮🏮
		if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
		} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
		} else {
			dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
			dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
			write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
		}
		if(nrow(dat.X.iv) ==0) {write(paste("SKIP:", X,"NA",Y, "X no IV !!"), file=log_mr, append=TRUE); next}	
		dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
		dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
		write.table(dat.XY, paste0(label, '.dat'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
		fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2Y 
			beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- fit.X2Y$pval
			p.hetero <- mr_heterogeneity(dat.XY)$Q_pval
            p.pleio <- mr_pleiotropy_test(dat.XY)$pval
			write("X Y M BETA SE P.ivw p.hetero.egger p.hetero.weighted p.pleio", file=log_mr, append=FALSE)
            write(paste(X,Y, nrow(dat.XY), rb(beta.X2Y), rb(se.X2Y), rp(p.X2Y), rp(p.hetero[1]), rp(p.hetero[2]), rp(p.pleio)), file=log_mr, append=TRUE)

		## 中介分析
		next # 🛑 if(p.X2Y >0.05)
		# numCores <- detectCores()-1; cl <- makeCluster(numCores); registerDoParallel(cl)
		# foreach::foreach(M = Ms) %dopar% { # 🐎
		for (M in Ms) { mediation(X, dat.X.raw, dat.X.iv, M, Y, dat.Y.raw, log_mrMed) }
		# stopCluster(cl)
	}
}

