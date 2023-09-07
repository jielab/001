pacman::p_load(readr, dplyr, tidyr, tidyverse, ggplot2, cowplot, TwoSampleMR, RMediation)

ieu_exposure <- 'ukb-b-5192' # 看电视
ieu_outcome <- 'ieu-a-966' # 肺癌
ieu_mediator <- 'ieu-b-25' # 吸烟
dat_exposure <- extract_instruments(outcomes = ieu_exposure)
dat_outcome <- extract_outcome_data(snps=dat_exposure$SNP, outcomes = ieu_outcome)

# !!! 将来，下面这些代码，不需要改动
# 计算方法
difference_method_PoE <- function(total_beta, total_se, direct_beta, direct_se){
	indirect_beta = total_beta -  direct_beta
	indirect_se = round(sqrt(total_se^2 + direct_se^2), 4)
	df <- data.frame(b= indirect_beta, se = indirect_se)
	df$lo_ci <- df$b - 1.96 * df$se; df$up_ci <- df$b + 1.96 * df$se; df$or <- exp(df$b); df$or_lci95 <- exp(df$lo_ci); df$or_uci95 <- exp(df$up_ci); df<-round(df,3)
	return(df)
}
product_method_Delta <- function(EM_beta, EM_se, MO_beta, MO_se){
	EO <- EM_beta * MO_beta
	CIs = medci(EM_beta, MO_beta, EM_se, MO_se, type="dop")
	df <- data.frame(b = EO, se = CIs$SE, lo_ci = CIs[["95% CI"]][1], up_ci= CIs[["95% CI"]][2])
	df$or <- exp(df$b); df$or_lci95 <- exp(df$lo_ci); df$or_uci95 <- exp(df$up_ci); df<-round(df,3)
	return(df)
}
product_method_PoE <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
	EO_beta <- EM_beta * MO_beta
	EO_se = round(sqrt(EM_se^2 + MO_se^2), 4)
	df <- data.frame(b= EO_beta, se = EO_se)
	df$lo_ci <- df$b - 1.96 * df$se; df$up_ci <- df$b + 1.96 * df$se; df$or <- exp(df$b); df$or_lci95 <- exp(df$lo_ci); df$or_uci95 <- exp(df$up_ci); df <- round(df,3)
	return(df)
}
# 合并数据并计算总的 effect
dat <- harmonise_data(dat_exposure, dat_outcome)
	res <- mr(dat)
	exposure_total_beta <- res %>% filter(method == "Inverse variance weighted") %>% pull(b)
	exposure_total_se <- res %>% filter(method == "Inverse variance weighted") %>% pull(se)
# 差异法
mvmr_exposure_dat <- mv_extract_exposures(c(ieu_exposure, ieu_mediator), clump_r2=0.001, clump_kb=10000, harmonise_strictness=2, access_token = ieugwasr::check_access_token(), find_proxies=TRUE, force_server=FALSE, pval_threshold=5e-8, pop="EUR") 
mvmr_outcome_dat <- extract_outcome_data(mvmr_exposure_dat$SNP, ieu_outcome)
mvmr_dat <- mv_harmonise_data(mvmr_exposure_dat, mvmr_outcome_dat, harmonise_strictness=2)
	mvmr_res <- mv_multiple(mvmr_dat)
	direct_beta <- mvmr_res[["result"]][["b"]][2]
	direct_se <- mvmr_res[["result"]][["se"]][2]
	indirect_effect <- difference_method_PoE(exposure_total_beta,exposure_total_se,direct_beta,direct_se); indirect_effect 
	indirect_effect[,1]/exposure_total_beta
# 乘积法
dat_outcome <- extract_outcome_data(snps = dat_exposure$SNP, outcomes = ieu_mediator) # mediator 作为 outcome
dat <- harmonise_data(dat_exposure, dat_outcome)
	res <- mr(dat)
	EM_beta <- res %>% filter(method == "Inverse variance weighted") %>% pull(b)
	EM_se <- res %>% filter(method == "Inverse variance weighted") %>% pull(se)
dat_exposure <- extract_instruments(outcomes = ieu_mediator) # mediator 作为 exposure
dat_outcome <- extract_outcome_data(snps = dat_exposure$SNP, outcomes = ieu_outcome) 
dat <- harmonise_data(dat_exposure, dat_outcome)
	res <- mr(dat)
	MO_beta_total <- res %>% filter(method == "Inverse variance weighted") %>% pull(b)
	MO_se_total <- res %>% filter(method == "Inverse variance weighted") %>% pull(se)
	MO_beta <- mvmr_res[["result"]][["b"]][1] #?
	MO_se <- mvmr_res[["result"]][["se"]][1] #?
product_method_Delta(EM_beta, EM_se, MO_beta_total, MO_se_total)
product_method_Delta(EM_beta, EM_se, MO_beta, MO_se)
product_method_PoE(EM_beta, EM_se, MO_beta_total, MO_se_total)
product_method_PoE(EM_beta, EM_se, MO_beta, MO_se)