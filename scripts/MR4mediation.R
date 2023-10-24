pacman::p_load(tidyverse, TwoSampleMR, RMediation, MVMR)

ieu_X <- 'ukb-b-5192'; ieu_M <- 'ieu-b-25'; ieu_Y <- 'ieu-a-966' # watching TV, smoking, lung cancer
dat_X <- extract_instruments(outcomes=ieu_X)
dat_M_exp <- extract_instruments(outcomes=ieu_M)
dat_M_out <- extract_outcome_data(snps=dat_X$SNP, outcomes=ieu_M)
dat_Y <- extract_outcome_data(snps=dat_X$SNP, outcomes=ieu_Y)
dat_XY <- harmonise_data(dat_X, dat_Y)
res <- mr(dat_XY)
beta_X2Y <- res %>% filter(method=="Inverse variance weighted") %>% pull(b)
se_X2Y <- res %>% filter(method=="Inverse variance weighted") %>% pull(se)
p_X2Y=2*pnorm(abs(beta_X2Y/se_X2Y), lower.tail=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Two-step method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #1：from X to M
dat_XM <- harmonise_data(dat_X, dat_M_out)
res <- mr(dat_XM)
beta_X2M <- res %>% filter(method=="Inverse variance weighted") %>% pull(b)
se_X2M <- res %>% filter(method=="Inverse variance weighted") %>% pull(se)
# #2：from M to Y
mvmr_X_dat <- mv_extract_Xs(c(ieu_X, ieu_M), harmonise_strictness=2) 
dat_XMY <- mv_harmonise_data(dat_XM, dat_Y, harmonise_strictness=2)
res <- mv_multiple(dat_XMY)
beta_M2Y <- mv_res[["result"]][["b"]][2]
se_M2Y <- mv_res[["result"]][["se"]][2]
# X2M * M2Y
beta_X2M2Y <- beta_EM * beta_MO
CIs = medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y, type="dop")
se_X2M2Y = CIs$SE
p=2*pnorm(abs(beta_X2M2Y/se_X2M2Y), lower.tail=F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MVMR method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- format_mvmr(BXGs = mvmr[,c(3,4)], BYG = mvmr[,1], seBXGs = mvmr[,c(5,6)], seBYG = mvmr[,2], RSID = "NULL")
sres <- strength_mvmr(r_input = F.data, gencov = 0) # test for weak instruments
pres <- pleiotropy_mvmr(r_input = F.data, gencov = 0) # test for horizontal pleiotropy using conventional Q-statistic estimation
res <- ivw_mvmr(r_input = F.data) # estimate causal effects