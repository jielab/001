pacman::p_load(readxl, tidyverse, TwoSampleMR)

X="" # no IEU id, use local file
M='ebi-a-GCST006867' # T2DM
Y='ebi-a-GCST90091033' # NAFLD
dat_X <- as.data.frame(read_excel("D:/Downloads/阜外ST7.LPL.xlsx")) %>% 
	mutate (
		CI_lower = as.numeric(CI_lower), 
		CI_upper = as.numeric(CI_upper),
		SE = (CI_upper - CI_lower)/(1.96*2) # 或者用 abs(BETA/qnorm(P/2))
	) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure="LPL")
dat_M4x <- extract_outcome_data(dat_X$SNP, M)
dat_M4x <- extract_outcome_data(dat_X$SNP, 'ebi-a-GCST006867', proxies=T, rsq=0.1, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(), splitsize=10000, proxy_splitsize=500)
dat_M <- extract_instruments(outcomes=M)
dat_Y4x <- extract_outcome_data(dat_X$SNP, Y)
dat_Y4m <- extract_outcome_data(dat_M$SNP, Y)
mr_X2M <- harmonise_data(dat_X, dat_M4x) %>% mr(); mr_X2M
mr_M2Y <- harmonise_data(dat_M, dat_Y4m) %>% mr(); mr_M2Y
mr_X2Y <- harmonise_data(dat_X, dat_Y4x) %>% mr(); mr_X2Y
