pacman::p_load(readr, dplyr, tidyr, tidyverse, ggplot2, cowplot, TwoSampleMR, RMediation)

ieu_exposure <- 'ukb-b-5192' # 看电视
ieu_outcome <- 'ieu-a-966' # 肺癌
ieu_mediator <- 'ieu-b-25' # 吸烟
dat_exposure <- extract_instruments(outcomes = ieu_exposure)
dat_outcome <- extract_outcome_data(snps=dat_exposure$SNP, outcomes = ieu_outcome)

# 计算方法
difference_method <- function(total_beta, total_se, direct_beta, direct_se){
	indirect_beta = total_beta - direct_beta # 差异法
	indirect_se = sqrt(total_se^2 + direct_se^2)
	df <- data.frame(beta=indirect_beta, se=indirect_se); df$p=2*pnorm(abs(df$beta/df$se), lower.tail=F)
	# lo_ci <- beta-1.96*se, up_ci <- beta+1.96*se, or <- exp(beta), or_lci95 <- exp(lo_ci), or_uci95 <- exp(up_ci); 
	return(df)
}
product_method <- function(EM_beta, EM_se, MO_beta, MO_se){
	EO_beta <- EM_beta * MO_beta # 乘积法
#	EO_se = sqrt(EM_se^2 + MO_se^2) # 这是一种简单的算法，建议用下一行的算法
	CIs = medci(EM_beta, MO_beta, EM_se, MO_se, type="dop"); EO_se = CIs$SE # lo_ci = CIs[["95% CI"]][1], up_ci= CIs[["95% CI"]][2]
	df <- data.frame(beta=EO_beta, se=EO_se); df$p=2*pnorm(abs(df$beta/df$se), lower.tail=F)
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
	indirect_effect <- difference_method(exposure_total_beta, exposure_total_se, direct_beta, direct_se); indirect_effect 
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
	MO_beta <- mvmr_res[["result"]][["b"]][1]
	MO_se <- mvmr_res[["result"]][["se"]][1]
product_method(EM_beta, EM_se, MO_beta, MO_se) # 或者 (EM_beta, EM_se, MO_beta_total, MO_se_total)