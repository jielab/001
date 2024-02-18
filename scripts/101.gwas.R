#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Multiple gwas manhattan plot
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
traits=c("height", "leg_ratio", "chunk", "leg")
for (trait in traits) {
	datmp <- read.table(paste0(dir, 'me.', trait, '.a.gz'), header=TRUE)[,-c(6,7)] %>% 
		mutate(P=ifelse(P<1e-50, 1e-50, P)) %>% filter(P<1e-07, P!=0) %>% na.omit() 
		print(paste(trait, nrow(datmp), ncol(datmp)))
	eval(parse(text=paste0('dat_', trait, '.topr <- datmp %>% rename(REF=NEA, ALT=EA)')))
	eval(parse(text=paste0('dat_', trait, '.cmpt <- datmp %>% select(SNP, CHR, POS, P) %>% rename(P.', trait, '=P)') ))
	eval(parse(text=paste0('dat_', trait, '.pheg <- datmp %>% select(CHR, POS, SNP, EA, NEA, BETA, SE, P) %>% rename(BP=POS, rsid=SNP, A1=EA, A2=NEA, beta=BETA, se=SE, P.value=P)') ))
}
dat_topr <- list(dat_height.topr, dat_leg_ratio.topr, dat_chunk.topr, dat_leg.topr); # eval(parse(text=paste0( 'dat_', traits, '.topr'   ) ))
	png("topr.png", w=2000, h=1000, res=128); topr::manhattan(dat_topr, ntop=1, legend_labels=traits, color=c('grey','blue','orange','green')); dev.off()
dat_cmpt <- Reduce(function(x,y) merge(x,y,by="SNP",all=T), list(dat_height.cmpt, dat_leg_ratio.cmpt, dat_chunk.cmpt, dat_leg.cmpt)) 
	dat_cmpt$CHR <- apply(subset(dat_cmpt, select=grepl("CHR", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
	dat_cmpt$POS <- apply(subset(dat_cmpt, select=grepl("POS", names(dat_cmpt))), 1, FUN=min, na.rm=TRUE)
	dat_cmpt <- subset(dat_cmpt, select=c('SNP','CHR','POS', grep('^P\\.', names(dat_cmpt),value=TRUE)))
	CMplot(dat_cmpt, plot.type="m", multracks=TRUE, cex=0.5, amplify=FALSE, file.output=TRUE, file="jpg", file.name="", width=20, height=4, dpi=300) 
dat_pheg <- list(dat_height.pheg, dat_leg_ratio.pheg, dat_chunk.pheg, dat_leg.pheg)
	dat_pheg <- PheGWAS::processphegwas(dat_pheg, traits)
	landscape(dat_pheg, chromosome=19)