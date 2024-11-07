# 参考 https://github.com/RitchieLab/Gene-level-statistical-colocalization
pacman::p_load(tidyverse, TwoSampleMR, coloc, gassocplot) # data(coloc_test_data)
pattern=c('^snp$|^rsid$', '^chr$|^chrom$|^chromosome$', '^bp$|^pos$|^position$|^base_pair_location$', '^ea$|^alt$|^a1$|^effect_allele$', '^nea$|^ref|^allele0$|^a2$|^other_allele$', '^eaf$|^a1freq$|^effect_allele_frequency$', '^n$|^Neff$', '^beta$|^effect$', '^se$|^standard_error', '^p$|^pval$|^p_value$')
replacement=c('SNP', 'CHR', 'POS', 'EA', 'NEA', 'EAF', 'N', 'BETA', 'SE', 'P')
rb <- function(x) (round(x,3)); rp <- function(x) (signif(x,2))

CHR=9; flank=100000; BEGIN=136130562-flank; END=136150630+flank
dat.X.file <- 'D:/data/gwas/main/clean/ABO.gz'
dat.Y.file <- 'D:/data/gwas/main/clean/x.bmi.gz'


dat.X.raw <- read.table(dat.X.file, header=T)
dat.Y.raw <- read.table(dat.Y.file, header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
dat.X <- dat.X.raw %>% filter(CHR==CHR, POS>=BEGIN, POS<=END) %>% format_data(type='exposure', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P') 
dat.Y <- dat.Y.raw %>% filter(CHR==CHR, POS>=BEGIN, POS<=END) %>% format_data(type='outcome', snp_col='SNP', chr_col='CHR', pos_col='POS', effect_allele_col='EA', other_allele_col='NEA', eaf_col='EAF', samplesize_col='N', beta_col='BETA', se_col='SE', pval_col='P') 
dat <- harmonise_data(dat.X, dat.Y, action=1) 
	dat.coloc.X <- dat %>% {list(type="quant", snp=.$SNP, position=.$pos.exposure, MAF=.$eaf.exposure, N=.$samplesize.exposure, beta=.$beta.exposure, varbeta =.$se.exposure^2, pvalues=.$pval.exposure)}
	dat.coloc.Y <- dat %>% {list(type="cc", snp=.$SNP, position=.$pos.outcome, MAF=.$eaf.outcome, N=100000, beta=.$beta.outcome, varbeta =.$se.outcome^2, pvalues=.$pval.outcome)}
fit.coloc <- coloc::coloc.abf(dat.coloc.X, dat.coloc.Y)
	t(as.data.frame(fit.coloc[["summary"]])) 
par(mfrow = c(2,1)); plot_dataset(dat.coloc.X); plot_dataset(dat.coloc.Y)
	sensitivity(fit.coloc,"H4>0.9")
