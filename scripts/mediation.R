#install_github("WSpiller/MVMR", build_opts=c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
pacman::p_load(data.table, dplyr, tidyverse, TwoSampleMR, MendelianRandomization, MVMR, RMediation)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read and harmonize data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 如果不用ieu数据，可先下载GWAS到本地，然后用下一行代码
# dat_Y <- read.table("D:/data/gwas/pheno/y.vte.gz", sep="\t", header=T) %>% merge(subset(dat_X, select="SNP")) %>% format_data(type ="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
ieu_X='ukb-b-10787' # 站高; ukb-b-16881 坐高 
ieu_M='ieu-b-25' # Cystain-C 胱抑素C，creatinine 肌酐 
ieu_Y='ieu-a-966' # 
dat_X <- extract_instruments(outcomes=ieu_X, clump=T)
dat_M4x <- extract_outcome_data(dat_X$SNP, ieu_M) #proxies=TRUE, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(), splitsize=10000, proxy_splitsize=500)
dat_M4y <- extract_instruments(outcomes=ieu_M, clump=T)
dat_Y4x <- extract_outcome_data(dat_X$SNP, ieu_Y)
dat_Y4m <- extract_outcome_data(dat_M4y$SNP, ieu_Y)
dat_XM <- harmonise_data(dat_X, dat_M4x, action=2)
dat_MY <- harmonise_data(dat_M4y, dat_Y4m, action=2) 
dat_XY <- harmonise_data(dat_X, dat_Y4x, action=2) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# two-step MR 方法 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
res_X2Y <- TwoSampleMR::mr(dat_XY) %>% filter(method=="Inverse variance weighted") # total 总的
	beta_X2Y <- res_X2Y %>% pull(b); se_X2Y <- res_X2Y %>% pull(se); p_X2Y <- res_X2Y %>% pull(pval)
res_X2M <- TwoSampleMR::mr(dat_XM) %>% filter(method=="Inverse variance weighted") # 1st step 三角形的上坡
	beta_X2M <- res_X2M %>% pull(b); se_X2M <- res_X2M %>% pull(se); p_X2M <- res_X2M %>% pull(pval)
res_M2Y <- TwoSampleMR::mr(dat_MY) %>% filter(method=="Inverse variance weighted") # 1st step 三角形的下坡
	beta_M2Y <- res_M2Y %>% pull(b); se_M2Y <- res_M2Y %>% pull(se); p_M2Y <- res_M2Y %>% pull(pval)
beta_2step <- beta_X2M * beta_M2Y; beta_2step  # *** 核心结果，乘法 ！！！
	CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y); se_2step = CIs$SE 
	p_2step <- 2*pnorm(abs(beta_2step/se_2step), lower.tail=F); p_2step


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MVMR方法，建议不用，仅供参考 ！！  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Use TwoSampleMR
dat_X <- mv_extract_exposures(c("ieu-a-299", "ieu-a-300", "ieu-a-302")) 
dat_Y <- extract_outcome_data(dat_X$SNP, "ieu-a-7") 
dat <- mv_harmonise_data(dat_X, dat_Y)
res <- mv_multiple(dat)
# Use WSpiller/MVMR
dat_XM0 <- rbind(dat_X0, dat_M0) %>% clump_data() 
dat_X <- dat_X %>% convert_outcome_to_exposure() # extract_instruments()不允许用户提供指定的SNP，所以用了两步
dat_X3c <- dat_X %>% dplyr::select(SNP, beta.exposure, se.exposure) 
dat_M3c <- dat_M %>% dplyr::select(SNP, beta.outcome, se.outcome) %>% rename(beta.mediator=beta.outcome, se.mediator=se.outcome)  
dat_Y3c <- dat_Y %>% dplyr::select(SNP, beta.outcome, se.outcome)
dat0 <- merge(merge(dat_X3c, dat_M3c), dat_Y3c)
dat <- MVMR::format_mvmr(RSID=dat0$SNP, BXGs=dat0[,c("beta.exposure","beta.mediator")], BYG=dat0$beta.outcome, seBXGs=dat0[,c("se.exposure","se.mediator")], seBYG=dat0$se.outcome)
res_mvmr <- MVMR::ivw_mvmr(r_input=dat)
	res_strength <- strength_mvmr(r_input=dat) 
	res_pleiotropy <- MVMR::pleiotropy_mvmr(r_input=dat) 
	res_egger <- MendelianRandomization::mr_mvegger(mr_mvinput(bx=cbind(dat$betaX1, dat$betaX2), bxse=cbind(dat$sebetaX1, dat$sebetaX2), by=dat$betaYG, byse=dat$sebetaYG), orientate=1)
	res_mvmedian <- MendelianRandomization::mr_mvmedian(mr_mvinput(bx=cbind(dat$betaX1, dat$betaX2), bxse=cbind(dat$sebetaX1, dat$sebetaX2), by=dat$betaYG, byse=dat$sebetaYG), iterations=1000)
	beta_X2Y_mvmr <- res_mvmr[1,1]; se_X2Y_mvmr <- res_mvmr[1,2]; p_X2Y_mvmr <- res_mvmr[1,4]
	beta_M2Y_mvmr <- res_mvmr[2,1]; se_M2Y_mvmr <- res_mvmr[2,2]; p_X2Y_mvmr <- res_mvmr[2,4]
beta_mvmr <- beta_X2Y - beta_X2Y_mvmr # *** 核心结果，减法 ！！！
	se_mvmr <- (se_X2Y^2 + se_X2Y_mvmr^2)^0.5
	p_mvmr <- 2*pnorm(abs(beta_mvmr/se_mvmr), lower.tail=F)

