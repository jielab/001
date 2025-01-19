pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')

rb <- function(x) (if(is.null(x)) return(NA) else round(x,3))
rp <- function(x) (if(is.null(x)) return(NA) else signif(x,2))
inormal <- function(x) qnorm((rank(x, na.last="keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
remove_outlier <- function(x) (ifelse((x > (mean(x,na.rm=TRUE) + 3*sd(x,na.rm=TRUE)) | x < (mean(x,na.rm=TRUE) - 3*sd(x,na.rm=TRUE))), NA, x))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))
expo <- function(x) 1.1^x


valcano <- function( # 🌋
  res, prot, BETA, P, sig.level) {
	res$BETA.std <- scale(res$BETA)
	res$color <- with(res, ifelse(P < sig.level & BETA.std > 0, "positive", ifelse(P < sig.level & BETA.std < 0, "negative", "NS")))
	res <- res[order(res$P), ]
	top5 <- head(res, 5); top5$label=toupper(gsub("prot_", "", top5[[prot]]))
	ggplot(res, aes(x=BETA.std, y=-log10(P), color=color)) + geom_point(size=2) +
		scale_color_manual(values=c("positive"="purple", "negative"="green", "NS"="gray")) +
		geom_text(data=top5, aes(label=label), hjust=0, vjust=-1.5, size=3.5, color="black", fontface="bold") +
	#	geom_segment(data=top5, aes(x=BETA.std + 0.05, y=-log10(P) + 1, xend=BETA.std, yend=-log10(P)), color="black") +
		geom_hline(yintercept=-log10(sig.level), linetype="dashed", color="red", linewidth=1.2) +
		geom_vline(xintercept=0, linetype="dotted", linewidth=1.2) + theme_minimal() +
		labs(x="Standardized beta", y="-log10(P)", title=Y) +
		theme(axis.text=element_text(size=12, face="bold"), axis.title=element_text(size=14, face="bold"), axis.line=element_line(linewidth=1.2), legend.position="none", plot.title=element_text(size=16, face="bold", hjust=0.5))
}


allele.qc <- function(a1, a2, ref1, ref2) {
    flip_single <- function(allele) {
	if (allele=="A") return("T")
	if (allele=="T") return("A")
	if (allele=="G") return("C")
	if (allele=="C") return("G")
	return(allele)
    }
	flip1 <- sapply(ref1, flip_single)
	flip2 <- sapply(ref2, flip_single)
	snp <- list()
	snp$keep <- !(a1=="A" & a2=="T") & !(a1=="T" & a2=="A") & !(a1=="C" & a2=="G") & !(a1=="G" & a2=="C")
	snp$keep[!grepl("^[ATGC]+$", a1)] <- F
	snp$keep[!grepl("^[ATGC]+$", a2)] <- F
	snp$exact_match <- (a1==ref1 & a2==ref2)
	snp$sign_flip <- (a1==ref2 & a2==ref1) | (a1==flip2 & a2==flip1)
	snp$strand_flip <- (a1==flip1 & a2==flip2) | (a1==flip2 & a2==flip1)
	snp$keep[!(snp$exact_match | snp$sign_flip | snp$strand_flip)] <- F
	return(snp)
}


IV.filter <- function(dat, CHR_col, POS_col, AF_col, BETA_col, N_col) {
	filter1 <- dat[[CHR_col]]==6 & dat[[POS_col]] >=28510120 & dat[[POS_col]] <=33480577
	dat$R2 <- 2 * dat[[AF_col]] * (1-dat[[AF_col]]) * (dat[[BETA_col]]^2)
	dat$F_stat <- dat$R2 * (dat[[N_col]]-2) / (1-dat$R2)
	filter2 <- dat$F_sta < 10
	return(dat[!(filter1 | filter2), ])
}


run_mr2s <- function( # 🏮 常规的 TwoSampleMR
  analysis, label, X, dir.X, dat.X.raw, IV.file, IV_filter, Y, dir.Y, dat.Y.raw) { 
  
	log_file <- paste0(label, ".mr.log")
	if (IV.file != "NO") {
		dat.X.iv <- read.table(IV.file, header=T); names(dat.X.iv) <- "SNP"
	} else if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
	} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
	} else {
		dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
		write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
	}

	if (nrow(dat.X.iv)==0) {write(paste(analysis, X,Y, "SKIP: data X has no IV !!"), file=log_file, append=TRUE); return(NULL)}	
	dat.X <- dat.X.raw %>% merge(dat.X.iv) 
		if (IV_filter==1) { dat.X <- IV.filter(dat.X, "CHR", "POS", "EAF", "BETA", "N") }
		if (nrow(dat.X)==0) {write(paste(analysis, X,Y, "SKIP: data X is empty after IV filter!!"), file=log_file, append=TRUE); return(NULL)}
		dat.X <- dat.X %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.Y.4x <- subset(dat.Y.raw, SNP %in% dat.X$SNP) 
		if(nrow(dat.Y.4x)==0) {write(paste(analysis, X,Y, "SKIP: IV SNPs don't exist in Y!!"), file=log_file, append=TRUE); return(NULL)}	
		dat.Y.4x <- dat.Y.4x %>% format_data(type='outcome', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P')
	dat.XY <- harmonise_data(dat.X, dat.Y.4x, action=1) 
	write.table(dat.XY, paste0(label, '.dat'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)

	fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw"))
	beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- fit.X2Y$pval
	p.hetero <- mr_heterogeneity(dat.XY)$Q_pval
	p.pleio <- mr_pleiotropy_test(dat.XY)$pval
	write(paste(analysis, X,Y, nrow(dat.XY), rb(beta.X2Y), rb(se.X2Y), rp(p.X2Y), rp(p.hetero[1]), rp(p.hetero[2]), rp(p.pleio), rp(min(dat.XY$pval.outcome))), file=log_file, append=TRUE)
}


run_cisMr <- function( # 🏮 参考 https://github.com/ZhaotongL/cisMR-paper/blob/main/process_ARIC.R
  label, chr, flank, pos0, pos1, X, dat.X.raw, dir.X.cojo, Y, dat.Y.raw) {

	p_t <- 5e-06
	log_file <- paste0(label, ".cisMr.log")
	
	cols <- c('SNP', 'POS', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'P', 'N')
	cols.y <- c(cols[1], paste0(cols[-1], ".y"))
	bfile <- paste("--bfile", paste0(dir.X.cojo, "/", X)) 
	
	dat.X <- dat.X.raw %>% filter(CHR==chr, POS>=(pos0-flank), POS<=(pos1+flank)) %>% select(all_of(cols))
	dat.Y <- dat.Y.raw %>% filter(SNP %in% dat.X$SNP) %>% mutate(N=100000) %>% select(all_of(cols))

	names(dat.Y) <- cols.y
	dat.XY <- merge(dat.X, dat.Y, by='SNP') %>% filter(!is.na(BETA), !is.na(BETA.y), !is.nan(BETA), !is.nan(BETA.y), BETA !=0, BETA.y !=0)
		if (nrow(dat.XY)==0) {write(paste('SKIP:', X, Y, 'merged data is empty'), file=log_file, append=TRUE); return(NULL)}
		remove_flip = allele.qc(dat.XY$EA, dat.XY$NEA, dat.XY$EA.y, dat.XY$NEA.y)
		dat.XY$BETA[which(remove_flip$sign_flip)] = -dat.XY$BETA[which(remove_flip$sign_flip)]
		dat.XY$EAF[which(remove_flip$sign_flip)] = 1-dat.XY$EAF[which(remove_flip$sign_flip)]
		dat.XY <- dat.XY[remove_flip$keep,]
		if (min(dat.XY$P) > 5e-06) { write(paste('SKIP:', X, Y, 'X .4gcta file has no significant SNPs'), file=log_file, append=TRUE); return(NULL) }
		dat.X <- dat.XY %>% select(all_of(cols)); dat.X$POS <- NULL # 后面的分析，不能有POS
		dat.Y <- dat.XY %>% select(all_of(cols.y)); names(dat.Y) <- cols; dat.Y$POS <- NULL
		fwrite(dat.X, paste0(X,'.4gcta'), sep='\t') 
		fwrite(dat.Y, paste0(Y,'.4gcta'), sep='\t') 
					
	X.cojo.cmd <- paste("gcta", bfile, "--cojo-file", paste0(X,'.4gcta'), "--cojo-slct --cojo-p ", p_t, " --out", X)
		X.cojo.fn <- paste0(X, '.jma.cojo')	
		ecode = system(X.cojo.cmd, intern=FALSE); if(ecode!=0) { write(paste('ERROR:', X, Y, '--cojo-slct on X does not run'), file=log_file, append=TRUE); return(NULL) }
		X.cojo <- fread(X.cojo.fn); if(length(X.cojo$SNP)<3) {write(paste('SKIP:', X, 'jma.cojo files less than 3 records'), file=log_file, append=TRUE); return(NULL)}
	
	Y.cojo.cmd <- paste("gcta", bfile, "--cojo-file", paste0(Y,'.4gcta'), "--cojo-slct --cojo-p ", p_t, " --out", Y); system(Y.cojo.cmd, intern=FALSE)
		Y.cojo.fn <- paste0(Y, '.jma.cojo')
		if (file.exists(Y.cojo.fn) && file.info(Y.cojo.fn)$size > 0) {Y.cojo=fread(Y.cojo.fn); IV=union(X.cojo$SNP, Y.cojo$SNP)} else {IV=X.cojo$SNP}
		write.table(IV, paste0(X, '.iv'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)
	
	XY.cojo.cmd = paste("gcta", bfile, "--cojo-file", paste0(X,'.4gcta'), "--extract", paste0(X, ".iv --cojo-joint --out"), X) # 🏮
		ecode = system(XY.cojo.cmd, intern=FALSE); if(ecode==1) { write(paste('ERROR:', X, Y, '--extract --cojo-joint on IV does not run'), file=log_file, append=TRUE); return(NULL) }
	
	LD_mat = fread(paste0(X, '.ldr.cojo'))
	LD_mat = LD_mat[,2:(ncol(LD_mat)-1)]; LD_mat = as.matrix(LD_mat); rownames(LD_mat) = colnames(LD_mat)
	
	dat.X <- dat.X %>% filter(SNP %in% colnames(LD_mat))
		dat.X$cor = dat.X$BETA / sqrt(dat.X$BETA^2 + (dat.X$N-2) * dat.X$SE^2)
		dat.X$se_cor = dat.X$SE * dat.X$cor/dat.X$BETA
		dat.X$bJ = solve(LD_mat) %*% dat.X$cor
	
	dat.Y = dat.Y %>% filter(SNP %in% colnames(LD_mat))
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
	write(paste(X, Y, nrow(dat.X), rb(res$BIC_DP_theta), rb(res$BIC_DP_se), rp(res$BIC_DP_p)), file=log_file, append=TRUE)

	if (res$BIC_DP_p > 1e-05) {return(NULL)}
	log_file <- paste0(label, ".coloc.log")
	write("X Y n-SNP H0 H1 H2 H3 H4", file=log_file, append=FALSE)		
	if (min(dat.XY$P.y) >1e-5) { write(paste(X, Y, min(dat.XY.$P.y), "SKIP: the minimum P-value in data Y is not significant"), file=log_file, append=TRUE); return(NULL)}
	dat.X.coloc <- dat.XY %>% {list(type="quant", snp=.$SNP, position=.$POS, MAF=.$EAF, N=.$N, beta=.$BETA, varbeta =.$SE^2, pvalues=.$P)}
	dat.Y.coloc <- dat.XY %>% {list(type="cc", snp=.$SNP, position=.$POS.y, MAF=.$EAF.y, N=.$N.y, beta=.$BETA.y, varbeta =.$SE.y^2, pvalues=.$P.y)}
	fit.coloc <- coloc::coloc.abf(dat.X.coloc, dat.Y.coloc)
	res <- paste(signif(fit.coloc$summary, 3),collapse=" ") # nsnps H0 H1 H2 H3 H4 
	write(paste(X, Y, nrow(dat.XY), res), file=log_file, append=TRUE)
	if (fit.coloc$summary[6] >0.5) { pdf(paste0(label,".coloc.pdf")); sensitivity(fit.coloc,"H4>0.5"); dev.off() }
}


run_mrMed <- function( # 🏮 devtools::install_github("scllin/mrMed") 
  label, X, dat.X.raw, file.M, Y, dat.Y.raw) { 

	log_file <- paste0(label, ".mrMed.log")
	dir.M <- dirname(file.M)
	M <- sub('\\.gz$', '', basename(file.M))
	
	if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
	} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
		dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
	} else {
		dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
		dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
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
		write(paste(X, Y, M, "SKIP: there is no SNP in dat.M.mv"), file=log_file, append=TRUE); return()
	} else {
		dat.M.mv <- dat.M.mv %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =M)
	}
	dat.Y.mv <- dat.Y.raw %>% merge(dat.XnM.iv) %>% format_data(type='outcome',  snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
	dat.XY.mv <- harmonise_data(dat.X.mv, dat.Y.mv, action=1)
	dat.XM.mv <- harmonise_data(dat.X.mv, dat.M.mv, action=1); names(dat.XM.mv) <- gsub('outcome', 'mediator', names(dat.XM.mv))
	dat <- merge(dat.XM.mv, dat.XY.mv, by='SNP')
	if (nrow(dat)==0) { write(paste(X, Y, M, "SKIP: X-M-Y harmonized data is empty"), file=log_file, append=TRUE); return()}
	bad_row <- subset(dat, effect_allele.exposure.x != effect_allele.exposure.y) %>% nrow()
	if (bad_row !=0) { write(paste(X, Y, M, "ERR: X-M-Y have inconsistent alleles"), file=log_file, append=TRUE); return()}
	names(dat) <- gsub("\\.x$", "", names(dat)); dat <- subset(dat, select=!grepl("\\.y", names(dat)))
	
	dat1 <- dat # mrMed 方法 
	names(dat1) <- stri_replace_all_regex(names(dat1), pattern=c("exposure", "mediator", "outcome"), replacement=c("X", "M", "Y"), vectorize_all=FALSE)
	dat1 <- dat1 %>% mutate (
		Gx = ifelse(dat$SNP %in% IV$SNP, 1, 0), Gx_plum = Gx, # 可加上 “& pval.M >5e-08”, 用于 X -> M或Y
		Gm = 1, Gm_plum = Gm, # 可考虑 ifelse(dat$SNP %in% dat.M.iv$SNP, 1, 0), 用于 M -> Y 
		G_mvmr = ifelse(Gx_plum==0 & Gm_plum==0, 0, 1)
	)
	
	if ( max(dat1$Gx) * max(dat1$Gm)==0 ) { write(paste(X, Y, M, "mrMed data has 0 rows for Gx or Gm"), file=log_file, append=TRUE); return() }
	res <- try( {mrMed(dat_mrMed=dat1, method_list="Prod_IVW_0")}, silent=TRUE )	
	if (inherits(res, "try-error")) { write(paste(X, Y, M, "ERR: mrMed gives ERROR message"), file=log_file, append=TRUE); return() }
	X_str=paste0(X, "(", nrow(IV), ")"); M_str=paste0(M, "(", nrow(dat.M.iv), ")");
	write(paste(X_str, M_str, Y, nrow(dat.XY), rb(beta.X2Y), rp(p.X2Y), 
		paste(rb(res$TE$b), rp(res$TE$p), rb(res$IE$b), rp(res$IE$p), rb(res$DE$b), rp(res$DE$p), rb(res$rho$b), rp(res$rho$p))
		), file=log_file, append=TRUE
	)

}
