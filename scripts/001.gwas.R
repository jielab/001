# for dat in `ls *gz | sed 's/\.gz//g'`; do echo $dat; zcat $dat.gz | awk 'NR==1 || $NF<1e-03' > $dat.1e-3; done
pacman::p_load(dplyr, tidyverse, plotly, qqman, topr, CMplot) # PheGWAS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multiple gwas manhattan plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dir="D:/data/gwas/main/clean"
dat <- read.table(paste0(dir, '/a12.1e-3'), header=TRUE) %>% 
	filter(P<1e-03, EAF>0.001 & EAF <0.999) %>% 
	mutate(CHR=ifelse(CHR=="X",23,CHR), CHR=as.numeric(CHR), P=ifelse(P<1e-50,1e-50,P))
	qqman::manhattan(dat, chr="CHR", bp="POS", p="P", snp="SNP", col = c("blue4", "orange3"))
labels=c("bald12", "bald13", "bald14", "bald1234")
for (label in labels) {
	datmp <- read.table(paste0(dir, '/', label, '.1e-3'), header=TRUE) %>%
		filter(P<1e-03, P!=0) %>% na.omit() %>% 
		dplyr::select(SNP, CHR, POS, EA, NEA, BETA, SE, P) %>%
		mutate(P=ifelse(P<1e-50, 1e-50, P)) 
	print(paste(label, nrow(datmp), ncol(datmp)))
	eval(parse(text=paste0('dat_', label, '.cmpt <- datmp %>% select(SNP, CHR, POS, P) %>% rename(P.', label, '=P)') ))
	assign(paste0('dat_', label, '.topr'), datmp %>% rename(REF=NEA, ALT=EA))
}
dat_cmpt <- Reduce(function(x, y) merge(x, y, by ="SNP", all = TRUE), lapply(labels, function(label) get(paste0("dat_", label, ".cmpt"))))
dat_cmpt$CHR <- apply(subset(dat_cmpt, select=grepl("CHR", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
dat_cmpt$POS <- apply(subset(dat_cmpt, select=grepl("POS", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
dat_cmpt <- subset(dat_cmpt, select=c('SNP','CHR','POS', grep('^P\\.', names(dat_cmpt),value=TRUE)))
CMplot(dat_cmpt, plot.type="m", multracks=TRUE, cex=0.2, amplify=FALSE, file.output=TRUE, file="jpg", file.name="", width=20, height=4, dpi=300) 
dat_topr <- list(mget(paste0('dat_', labels, '.topr')))
eval(parse(text=paste0('dat_topr <- list(', toString(paste0('dat_', labels, '.topr')), ')') ))
png("topr.png", w=2000, h=1000, res=128); topr::manhattan(dat_topr, ntop=2, legend_labels=labels, color=c('grey','blue','orange','green')); dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# locuszoomr
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
