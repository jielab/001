setwd("D:/")
pacman::p_load(tidyverse, dplyr, plotly, topr, CMplot, PheGWAS)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multiple gwas manhattan plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir="D:/data/gwas/main/"
labels=c("height", "chunk", "leg", "chunk_ratio")
for (label in labels) {
	datmp <- read.table(paste0(dir, 'me.', label, '.gz'), header=TRUE) %>%
		dplyr::select(SNP, CHR, POS, EA, NEA, BETA, SE, P) %>%
		mutate(P=ifelse(P<1e-50, 1e-50, P)) %>% filter(P<1e-07, P!=0) %>% na.omit() 
		print(paste(label, nrow(datmp), ncol(datmp)))
	eval(parse(text=paste0('dat_', label, '.topr <- datmp %>% rename(REF=NEA, ALT=EA)')))
	eval(parse(text=paste0('dat_', label, '.cmpt <- datmp %>% select(SNP, CHR, POS, P) %>% rename(P.', label, '=P)') ))
	eval(parse(text=paste0('dat_', label, '.pheg <- datmp %>% rename(BP=POS, rsid=SNP, A1=EA, A2=NEA, beta=BETA, se=SE, P.value=P)') ))
}
list(mget(paste0('dat_', labels))) ???
eval(parse(text=paste0('dat_topr <- list(', toString(paste0('dat_', labels, '.topr')), ')') ))
	png("topr.png", w=2000, h=1000, res=128); topr::manhattan(dat_topr, ntop=2, legend_labels=labels, color=c('grey','blue','orange','green')); dev.off()
eval(parse(text=paste0('dat_cmpt <- list(', toString(paste0('dat_', labels, '.cmpt')), ')') ))
	dat_cmpt <- Reduce(function(x,y) merge(x,y,by="SNP",all=T), dat_cmpt)
	dat_cmpt$CHR <- apply(subset(dat_cmpt, select=grepl("CHR", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
	dat_cmpt$POS <- apply(subset(dat_cmpt, select=grepl("POS", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
	dat_cmpt <- subset(dat_cmpt, select=c('SNP','CHR','POS', grep('^P\\.', names(dat_cmpt),value=TRUE)))
	CMplot(dat_cmpt, plot.type="m", multracks=TRUE, cex=0.5, amplify=FALSE, file.output=TRUE, file="jpg", file.name="", width=20, height=4, dpi=300) 
dat_pheg <- list(dat_AFR.pheg, dat_EAS.pheg, dat_EUR.pheg, dat_HIS.pheg, dat_HIS.pheg) 
	dat <- PheGWAS::processphegwas(dat_pheg, labels) %>% na.omit()
	landscape(dat, chromosome=19)
	landscape(dat, sliceval=10, chromosome=1, bpdivision=1e+06, betaplot=TRUE)
