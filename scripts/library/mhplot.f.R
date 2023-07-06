mhplot <-function( trait, gwas, chr="CHR", pos="POS", pval="P", cutoff=5e-8, 
	ylim_t=0, maf=NA, maf_t=1e-4, poscon_f="NO", dibiao=NA ) {

	nchr=23	
	if (grepl("\\.gz", gwas)) {
		dat <- read.table(gzfile(gwas,'r'), header=T, as.is=T)
	} else {
		dat <- read.table(gwas, header=T, as.is=T)
	}
	if (is.na(maf)) {
		dat <- dat[,c(chr,pos,pval)]; dat$maf=0.5
	} else {
		dat <- dat[,c(chr,pos,pval, maf)]
	}
	names(dat) <- c("chr","pos","pval", "maf")
	dat <- subset(dat, pval <=10^(-ylim_t))
	
	if(file.exists(dibiao)) {
		border=read.table(dibiao, header=F)
		names(border) <- c("chr","pos")
		border$maf=0.5
		border$pval=1
		dat <- rbind(dat,border)
	}	
	if (file.exists(poscon_f)) {
		poscon <- read.table(poscon_f, header=T,as.is=T)
		poscon <- poscon[,c("CHR","POS")]
		names(poscon) <- c("chr","pos")
	} else {
		poscon <- data.frame(chr=NA, pos=NA)[numeric(0),]
	}
	
	dat$chr[dat$chr=="X"] <- 23; dat$chr <- as.numeric(dat$chr)
	dat$pval <- as.numeric(dat$pval)
	dat <- subset(dat, !is.na(pval) & maf > maf_t & pval >0 & !is.na(pos) & chr %in% c(1:nchr))
	dat$maf =ifelse(dat$maf>0.5, 1-dat$maf, dat$maf)
	dat$logP<- -log10(dat$pval) 
	
	cumlen <- 0
	phy.max <- tapply(dat$pos, dat$chr, max,na.rm=T)
	for(i in 1:nchr){
		dat[dat$chr==i,"loc"]<- dat[dat$chr==i,"pos"]+cumlen
		poscon[poscon$chr==i,"loc"]<- poscon[poscon$chr==i,"pos"]+cumlen
		cumlen<-cumlen+phy.max[i]
	}
	phy.median <- tapply(dat$loc, dat$chr, median, na.rm=T)

	logP.max <- ceiling(max(dat$logP))
	dat$logP_NEW <- dat$logP
	dat$logP_NEW[dat$logP_NEW>10] <- 10+ (dat$logP_NEW[dat$logP_NEW>10]-10) * 10/(logP.max-10)
	ylim = ceiling(max(as.numeric(dat$logP_NEW)))

	dat$col = ifelse((dat$pval<=cutoff & dat$maf<=0.01), "red", ifelse(dat$chr %in% seq(2,24,2), "grey30", "grey60"))
	
	plot(dat$loc, dat[,"logP_NEW"], type="n", xaxt="n", frame.plot=F, main=paste(trait, "(N=",nrow(dat),")"), cex=1, cex.main=2, cex.axis=2, cex.lab=2, xlab="",ylab="-log10(P)", cex.lab=2, font=2, xlim=c(min(dat$loc,na.rm=T),max(dat$loc,na.rm=T)), ylim=c(ylim_t,ylim), axes=F)
	axis(side=1, at=phy.median, labels=1:nchr, font=2, cex.axis=2)
	if (logP.max > -log10(cutoff)) {
		abline(-log10(cutoff),0,col="red",lty=2)
	}
	if (logP.max <=10) {
		axis(2, font=2, cex.axis=2)
	} else {
	#	abline(10,0,col="green",lty=2)
		axis(2, at=c(0,5,10,ylim), labels=c(0,5,10,logP.max), font=2, cex.axis=2)
	}
	for(i in 1:nchr){
	#	if(i %in% seq(2,24,2)) col="grey" else col="darkgrey"
		dat1= subset(dat, chr==i) 
		points(dat1$loc, dat1$logP_NEW, col=dat1$col, pch=20, cex=0.8)
		abline(v=poscon[poscon$chr==i,"loc"], lwd=0.01, lty=1, col="green")
	}
}
