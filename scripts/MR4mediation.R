#install_github("WSpiller/MVMR", build_opts=c("--no-resave-data", "--no-manual"), build_vignettes=TRUE)
pacman::p_load(data.table, dplyr, tidyverse, TwoSampleMR, MendelianRandomization, MVMR, RMediation)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# read and harmonize data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ieu_X='ukb-b-10787' # 站高; ukb-b-16881 坐高 
ieu_M='ieu-b-25' # Cystain C  # creatinine 肌酐 
ieu_Y='36658437' # 
dat_X0 <- extract_instruments(outcomes=ieu_X, clump=F)
dat_M0 <- extract_instruments(outcomes=ieu_M, clump=F)
	dat_XM0 <- rbind(dat_X0, dat_M0)
	dat_XM <-clump_data(dat_XM0) 
dat_X <- extract_outcome_data(dat_XM$SNP, ieu_X, proxies=TRUE, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(), splitsize=10000, proxy_splitsize=500)
	dat_X <- convert_outcome_to_exposure(dat_X) # extract_instruments()不允许用户提供指定的SNP，所以用了两步
dat_M <- extract_outcome_data(dat_XM$SNP, ieu_M, proxies=TRUE, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(), splitsize=10000, proxy_splitsize=500)
dat_Y <- extract_outcome_data(dat_XM$SNP, ieu_Y, proxies=TRUE, rsq=0.8, align_alleles=1, palindromes=1, maf_threshold=0.3, access_token=ieugwasr::check_access_token(),  splitsize=10000, proxy_splitsize=500)
# 如果 Y 不在 ieu 里面，先下载GWAS到本地，然后用下面3行代码生成 dat_Y
	dat_Y0 <- read.table("D:/data/gwas/pheno/y.vte.gz", sep="\t", header=T) # 大约1千万行
	dat_Ytmp <- merge(dat_Y0, subset(dat_XM, select="SNP")) # 先merge变得很小后，再用下面的 format_data()命令
	dat_Y <- format_data(dat_Ytmp, type ="outcome", snp_col="SNP", effect_allele_col="EA", other_allele_col="NEA", beta_col="BETA", se_col="SE", pval_col="P") 
dat_XY <- harmonise_data(dat_X, dat_Y, action=2) 
dat_XM <- harmonise_data(dat_X, dat_M, action=2); dat_Mtmp <- convert_outcome_to_exposure(dat_M) # 要不然下面一行代码没法运行
dat_MY <- harmonise_data(dat_Mtmp, dat_Y, action=2) 


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
	CIs = RMediation::medci(beta_X2M, beta_M2Y, se_X2M, se_M2Y, type="dop")
	se_2step = CIs$SE # 也可简单用 (se_X2Y^2 + se_M2Y^2)^0.5
	p_2step <- 2*pnorm(abs(beta_2step/se_2step), lower.tail=F); p_2step


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# MVMR 方法 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
