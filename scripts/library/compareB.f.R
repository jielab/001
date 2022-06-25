compareB <-function(
	f1, f1_name, f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_se, f1_p,
	f2, f2_name, f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_se, f2_p
	){

	dat1 <- read.table(f1, header=T, as.is=T)
	dat2 <- read.table(f2, header=T, as.is=T) 
	f1_col_names <- as.character(names(dat1)); for (var in c('f1_snp', 'f1_ea', 'f1_nea', 'f1_eaf', 'f1_beta', 'f1_se', 'f1_p')) { eval(parse(text= paste0(var, '=grep(', var, ', f1_col_names, ignore.case=T, value=T)')[1]))}
	f2_col_names <- as.character(names(dat2)); for (var in c('f2_snp', 'f2_ea', 'f2_nea', 'f2_eaf', 'f2_beta', 'f2_se', 'f2_p')) { eval(parse(text= paste0(var, '=grep(', var, ', f2_col_names, ignore.case=T, value=T)')[1]))}
	#if (is.na(f1_nea)) {f1_nea='NEA'; dat1$NEA=NA}
	#if (is.na(f2_nea)) {f2_nea='NEA'; dat2$NEA=NA}
	if (length(f1_eaf)==0) {f1_eaf='EAF'; dat1$EAF=0}
	if (length(f2_eaf)==0) {f2_eaf='EAF'; dat2$EAF=0}
	if (is.null(dat1[[f1_beta]]) & !is.null(dat1$OR)) dat1[[f1_beta]]=log(dat1$OR)
	if (is.null(dat2[[f2_beta]]) & !is.null(dat2$OR)) dat2[[f2_beta]]=log(dat2$OR)
	if (is.na(f1_se)) {f1_se='SE'; dat1$SE=round(abs(dat1[[f1_beta]]/qnorm(dat1[[f1_p]]/2)),4); dat1$SEcal=NA} else {dat1$SEcal=round(abs(dat1[[f1_beta]]/qnorm(dat1[[f1_p]]/2)),4)}
	if (is.na(f2_se)) {f2_se='SE'; dat2$SE=round(abs(dat2[[f2_beta]]/qnorm(dat1[[f2_p]])),4); dat2$SEcal=NA} else {dat2$SEcal=round(abs(dat2[[f2_beta]]/qnorm(dat2[[f2_p]]/2)),4)}
	if (is.null(dat1[f1_p])) dat1[[f1_p]]=NA; dat1$logP <- -log10(dat1[[f1_p]])
	if (is.null(dat2[f2_p])) dat2[[f2_p]]=NA; dat2$logP <- -log10(dat2[[f2_p]])
	
	dat1 <- subset(dat1, select=c(f1_snp, f1_ea, f1_nea, f1_eaf, f1_beta, f1_se, 'SEcal', f1_p, 'logP'))
	dat2 <- subset(dat2, select=c(f2_snp, f2_ea, f2_nea, f2_eaf, f2_beta, f2_se, 'SEcal', f2_p, 'logP'))	
	names(dat1) <- paste0( c('SNP', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'SEcal', 'P', 'logP'), '_1')
	names(dat2) <- paste0( c('SNP', 'EA', 'NEA', 'EAF', 'BETA', 'SE', 'SEcal', 'P', 'logP'), '_2')
	dat <- merge(dat1, dat2, by.x='SNP_1', by.y='SNP_2')
	for (var in c('EA_1', 'NEA_1', 'EA_2', 'NEA_2')) { dat[[var]] <- toupper(dat[[var]]) }
	if(nrow(dat) > 1000) { dat <- subset(dat, P_1>5) }
	#dat <- na.omit(dat)
	#dat <- subset(dat, (EA_1==EA_2 & NEA_1==NEA_2) | (EA_1==NEA_2 & NEA_1==EA_2))
	#dat[which(dat$EA_1 != dat$EA_2), 'EAF_2']  <- 1- dat[which(dat$EA_1 != dat$EA_2), 'EAF_2']
	dat$EAF_2 = ifelse(dat$EA_1 == dat$EA_2, dat$EAF_2, 1- dat$EAF_2)
	dat$BETA_2 = ifelse(dat$EA_1 == dat$EA_2, dat$BETA_2, 0- dat$BETA_2)
	dat$BETA.aligned_1 = ifelse(dat$BETA_1>0, dat$BETA_1, 0- dat$BETA_1)
	dat$BETA.aligned_2 = ifelse(dat$BETA_1>0, dat$BETA_2, 0- dat$BETA_2)	
	write.table(dat, paste(f1_name,f2_name,'merged.txt',sep='.'), sep='\t', na='NA', append=F, quote=F, col.names=T, row.names=F)	
	
	pdf(paste0(f1_name,'.',f2_name,'.cB.pdf'))
	par(mfrow=c(3,2), mai=c(0.5,1,0.5,0.5))
	for (var in c('EAF', 'logP', 'BETA', 'BETA.aligned')) { 
		datt <- subset(dat, select=c(paste0(var,'_1'), paste0(var,'_2'))) # Var 1 on X-axis
		datt <- na.omit(datt)
		if (nrow(datt)>0) { # bty='n'
			plot(datt, main=paste0(var, ' (N=', nrow(dat), ')'), xlab=f1_name, ylab=f2_name, pch=20, cex=0.8, cex.axis=1.1, cex.lab=1.2, col=colors()[75], font=2, font.lab=2)
		}
	}
	plot(dat1$SE_1,dat1$SEcal_1, main=paste0(f1_name, ' (N=', nrow(dat1), ')'), xlab='SE original', ylab='SE calculated', pch=20, cex=0.8, cex.axis=1.1, cex.lab=1.2, col=colors()[75], font=2, font.lab=2)
	plot(dat2$SE_2,dat2$SEcal_2, main=paste0(f2_name, ' (N=', nrow(dat2), ')'), xlab='SE original', ylab='SE calculated', pch=20, cex=0.8, cex.axis=1.1, cex.lab=1.2, col=colors()[75], font=2, font.lab=2)
	dev.off()
	
	# MR
	library(MendelianRandomization); library(ggplot2)
	dat <- subset(dat, select=c('BETA.aligned_1', 'SE_1', 'BETA.aligned_2', 'SE_2'))
	dat <- na.omit(dat)
	mr.in <- mr_input(dat[,1], dat[,2], dat[,3], dat[,4])
#	mr_presso(data=dat, BetaOutcome='BETA_2', BetaExposure='BETA_1', SdOutcome='SE_2', SdExposure='SE_1', OUTLIERtest=T, DISTORTIONtest=T, NbDistribution=1000, SignifThreshold=0.05)
	mr.res <- paste(mr_ivw(mr.in)$Pvalue, mr_median(mr.in)$Pvalue, mr_maxlik(mr.in)$Pvalue, mr_mbe(mr.in)$Pvalue, mr_conmix(mr.in)$Pvalue, mr_egger(mr.in)$Pleio.pval, mr_egger(mr.in)$Causal.pval)
	mr.res.str <- paste(f1_name, f2_name, nrow(dat), mr.res) 
	write.table(mr.res.str, paste0(f1_name,'.',f2_name,'.mr.res.txt'), append=F, quote=F, row.names=F, col.names=F)
	jpeg(paste0(f1_name,'.',f2_name,'.mr.jpg'), w=1600, h=1000, res=256)
	mr_forest(mr.in, snp_estimates=F, methods=c('ivw', 'wmedian', 'maxlik', 'egger'))
	dev.off()
}