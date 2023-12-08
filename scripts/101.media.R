pacman::p_load(readxl, dplyr, tidyverse, TwoSampleMR, MVMR)

X = 'ebi-a-GCST90010149'  # LPL
M = 'ebi-a-GCST90018926' # T2D
Y = 'ebi-a-GCST90091033' # NAFLD 

dat_M.noClump <- extract_instruments(outcomes=M, clump=F)
dat_M.clumped <- dat_M.noClump %>% clump_data() 
dat_X.noClump <- extract_instruments(outcomes=X, clump=F)
dat_X.clumped <- dat_X.noClump %>% clump_data()

dat_M4x <- extract_outcome_data(dat_X.clumped$SNP, M)
dat_X2M <- harmonise_data(dat_X.clumped, dat_M4x)
dat_XnM.clumped <- rbind(dat_X.noClump, dat_M.noClump) %>% clump_data() 
dat_X.tmp <- extract_outcome_data(dat_XnM.clumped$SNP, X) %>% convert_outcome_to_exposure()
dat_M.tmp <- extract_outcome_data(dat_XnM.clumped$SNP, M) %>% convert_outcome_to_exposure()
dat_XnM <- rbind(dat_X.tmp, dat_M.tmp)
dat_XnM.2 <- mv_extract_exposures(c(X, M))

dat_Y4x <- extract_outcome_data(dat_X.clumped$SNP, Y)
dat_Y4m <- extract_outcome_data(dat_M.clumped$SNP, Y)
dat_Y4xm <- extract_outcome_data(dat_XnM$SNP, Y)
dat_Y4xm.2 <- extract_outcome_data(dat_XnM.2$SNP, Y)

# MVMR from TwoSampleMR Package
mv_harmonise_data(dat_XnM, dat_Y4xm) %>% mv_multiple()
mv_harmonise_data(dat_XnM.2, dat_Y4xm.2) %>% mv_multiple()

# Two-step MR 
res_X2Y <- harmonise_data(dat_X.clumped, dat_Y4x) %>% mr() %>% filter(method=='Inverse variance weighted') 
	beta_X2Y <- res_X2Y %>% pull(b); se_X2Y <- res_X2Y %>% pull(se); p_X2Y <- signif(res_X2Y %>% pull(pval),2)
res_X2M <- harmonise_data(dat_X.clumped, dat_M4x) %>% mr() %>% filter(method=='Wald ratio' | method=='Inverse variance weighted') # 1st step 三角形的上坡
	beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- res_X2M %>% pull(pval)
res_M2Y <- harmonise_data(dat_M.clumped, dat_Y4m) %>% mr() %>% filter(method=='Wald ratio' | method=='Inverse variance weighted') # 2nd step 三角形的下坡
	beta_M2Y <- res_M2Y %>% pull(b); se_M2Y <- res_M2Y %>% pull(se); p_M2Y <- res_M2Y %>% pull(pval)
beta_2step <- round(beta_X2M * beta_M2Y, 4); beta_2step # product method
	CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y, type='MC'); se_2step = CIs$SE; p_2step <- signif(2*pnorm(abs(beta_2step/se_2step), lower.tail=F), 3); p_2step			
	paste(X, M, Y, p_X2Y, p_2step)


