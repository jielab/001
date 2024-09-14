setwd("D:/")
pacman::p_load(data.table, tidyverse, lubridate, poLCA, BioAge)

dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 生物年龄 (PMID: 37080981)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#train=kdm_calc(biomarkers=var_kdm); kdm=kdm_calc(NHANES4, biomarkers=var_kdm, fit=train$fit, s_ba2=train$fit$s_ba2); data=kdm$data
var_kdm = c("fev1","sbp","bb_TC","bb_HBA1C","bb_ALB","bb_CRE","bb_CRP","bb_ALP","bb_BUN") # "fev","sbp","totchol","hba1c","albumin", "creat","lncrp","alp","bun
	kdm_male   = phe0 %>% filter (sex ==1) %>% kdm_calc(biomarkers=var_kdm, fit=kdm$fit$male,   s_ba2 = kdm$fit$male$s_b2)
	kdm_female = phe0 %>% filter (sex ==0) %>% kdm_calc(biomarkers=var_kdm, fit=kdm$fit$female, s_ba2 = kdm$fit$female$s_b2)
	kdm = rbind(kdm_female$data, kdm_male$data)
var_phenoage = c("bb_ALB","bc_LYMPH_pt","bc_MSCV","bb_GLU","bc_RDW","bb_CRE","bb_CRP","bb_ALP","bc_WBC") # "albumin_gL","lymph","mcv","glucose_mmol","rdw","creat_umol","lncrp","alp","wbc"
	phenoage_ukb = phe0 %>% phenoage_calc(biomarkers=var_phenoage, orig=TRUE) 
	phenoage_ukb = phenoage_ukb$data


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 社会经济地位（SES）
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
dat <- phe0 %>% dplyr::select(eid, age, sex, ethnic_cat, emp, income, edu) %>% filter(ethnic_cat=="White") %>% drop_na() #%>% dplyr::sample_n(10000) 
sink("ses.log")
SES_LCA_list <- list()
for (n_class in 2:6) { 
	SES_LCA_list[[n_class]] <- poLCA(cbind(edu_cat, emp_cat, income_cat) ~1, data=dat, nclass=n_class, maxiter=10000, tol=1e-6, nrep=2, graphs=TRUE, probs.start=NULL)
}
sink()
source("D:/scripts/library/00_LCA_out.R")
LCA_out(SES_LCA_list,3,3)
dat$ses <- SES_LCA_list[[3]]$predclass
saveRDS(subset(dat, select=c(eid, ses)), file="D:/data/ukb/Rdata/ukb.ses.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lifestyle 之 “管住嘴、 迈开腿、 躺下睡”
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 健康植物性饮食指数（PMID: 36976560）
# 南京医科大学版本（https://github.com/XiangyuYe/Infection_SES）
dat0 <- phe0 %>% 
	mutate(
		num_pa = (days_pa_mod > 5) & (days_pa_vig > 1),
		yes_pa_mod = ifelse(dura_pa_mod > 150, T, F),
		yes_pa_vig = ifelse(dura_pa_vig > 75, T, F),
		mat_pa = cbind(num_pa, yes_pa_mod, yes_pa_vig),
		pa_score = rep(0, nrow(mat_pa)),
		pa_score = ifelse(rowSums(mat_pa, na.rm = TRUE) > 0, 1, pa_score),
		pa_score = ifelse(rowSums(is.na(mat_pa)) == 3, NA, pa_score)
	) %>% mutate(
		cannabis_score = ifelse(cannabis == 0, 1, ifelse(cannabis >0, 0, NA)),
		noalcohol_stat = recode(alcohol_status, "never" = 0, "previous" = 1, "current" = 2),
		noalcohol_score = ifelse(noalcohol_stat == 0, 1, ifelse(noalcohol_stat > 0, 0, NA)),
		smoke_stat = recode(smoke_status, "never" = 0, "previous" = 1, "current" = 2),
		nosmoke_score = ifelse(smoke_stat == 0, 1, ifelse(smoke_stat >0 , 0, NA)),
		nosmoke_year = age_visit - age_smoke_quit,
		nosmoke_score = ifelse(nosmoke_year > 30, 1, nosmoke_score)
	) %>% mutate(
		fruit_score = rowSums(cbind(fruit_fresh, fruit_dried/5), na.rm = TRUE) >= 4,
		fruit_score = ifelse(is.na(fruit_fresh) & is.na(fruit_dried), NA, fruit_score),
		veg_score = rowSums(cbind(veg_cooked/3, veg_salad/3), na.rm = TRUE) >= 4,
		veg_score = ifelse(is.na(veg_cooked) & is.na(veg_salad) , NA, veg_score),
		fish_score = (fish_oily >= 3 | fish >= 3) | (fish_oily == 2 & fish == 2),
		pmeat_score = meat_processed <= 2,
		npmeat_score = rowSums(cbind(beef, lamb, pork)) <= 3 & !(beef >= 3 | lamb >= 3 | pork >= 3),  # label1 =0.5, then 1 + 2 OR 1 * 3
		npmeat_score = ifelse(is.na(beef) & is.na(lamb) & is.na(pork), NA, npmeat_score),
		bread_grain = ifelse(bread_type==3, bread, 0),
		cereal_grain = ifelse(cereal_type %in% c(1, 2, 3), cereal, 0),
		grains_score = rowSums(cbind(bread_grain, cereal_grain), na.rm = T) >= 3,
		grains_score = ifelse(is.na(bread_grain) & is.na(cereal_grain) , NA, grains_score),
		diet_mat = data.frame(fruit_score, veg_score, fish_score, pmeat_score, npmeat_score, grains_score),
		diet_score = rowSums(diet_mat) >= 4,
		diet_score = ifelse(rowSums(diet_mat, na.rm = TRUE) >= 4, TRUE, ifelse(rowSums(is.na(diet_mat)) + rowSums(diet_mat, na.rm = TRUE) < 4, FALSE, NA))
	) %>% mutate(
		chrono_score = ifelse(sleep_chronotype < 3, 1, 0), # Ref PMID:31848595
		duration_score = ifelse(sleep_duration <= 8 & sleep_duration >= 7, 1, 0), # Sleep duration 7-8h
		insomnia_score = ifelse(sleep_insomnia == 1, 1, 0), # 1: never or rarely insomnia
		snoring_score = ifelse(sleep_snoring == 2, 1, 0), # 2: No complaints of snoring
		narcolepsy_score = ifelse(sleep_narcolepsy < 2, 1, 0), # no frequently narcolepsy 0	Never/rarely 1	Sometimes
		sleep_mat = data.frame(chrono_score, duration_score, insomnia_score, snoring_score, narcolepsy_score),
		sleep_score = rowSums(sleep_mat) >= 4,
		sleep_score = ifelse(rowSums(sleep_mat, na.rm = TRUE) >= 4, TRUE, sleep_score),
		sleep_score = ifelse(rowSums(is.na(sleep_mat)) + rowSums(sleep_mat, na.rm = TRUE) < 4, FALSE, sleep_score)
	)
dat <- dat0 %>% dplyr::select(eid, smoke_score, pa_score, diet_score, alcohol_score, sleep_score, cannabis_score) *1
dat <- dat  %>% 
	mutate(
		tt4_score = rowSums(dat[,2:5], na.rm=T),
		lf4_score =ifelse(tt4_score >=3, "3-4", 
			ifelse(tt4_score ==2, "2", 
			ifelse(tt4_score <2, "0-1", NA))),
		tt6_score = rowSums(dat[,2:7], na.rm=T),
		tt6_score =ifelse(tt6_score >=4 | rowSums(dat[,2:6], na.rm=T)==3, "3", 
			ifelse(tt6_score >=2 & (tt6_score + rowSums(is.na(dat))) < 4, "2", 
			ifelse(tt6_score ==0, "1", "??")))
	)
saveRDS(life_factor_df, file = "D:/data/ukb/Rdata/ukb.le8.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# walk-VTE文章版本
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% 
	mutate( #🏮diet_score
		veg_cooked = ifelse(veg_cooked == 50, NA, veg_cooked),
		veg_salad = ifelse(veg_salad == 50, NA, veg_salad),
		vegetable = ifelse(veg_cooked + veg_salad >= 4, 1, 0),	#at least 4 tablespoons/day
		fruit_fresh = ifelse(fruit_fresh == 50, NA, fruit_fresh),
		fruit_dried = ifelse(fruit_dried == 100, NA, fruit_dried), fruit = ifelse(fruit_fresh + fruit_dried >= 3, 1, 0), #at least 3 pieces/day
		fish_oily = case_when(fish_oily %in% c(0, 1) ~ 0, fish_oily == 2 ~ 1, fish_oily == 3 ~ 2, fish_oily == 4 ~ 5,  fish_oily == 5 ~ 7),
		fish = case_when(fish %in% c(0, 1) ~ 0, fish == 2 ~ 1, fish == 3 ~ 2, fish == 4 ~ 5,  fish == 5 ~ 7),
		fish_all = ifelse(fish_oily + fish >= 2, 1, 0), #at least twice per week(>= 2 ~ 1)	
		beef = case_when(beef %in% c(0, 1) ~ 0, beef == 2 ~ 1, beef == 3 ~ 2, beef == 4 ~ 5,  beef == 5 ~ 7),
		lamb = case_when(lamb %in% c(0, 1) ~ 0, lamb == 2 ~ 1, lamb == 3 ~ 2, lamb == 4 ~ 5,  lamb == 5 ~ 7),
		pork = case_when(pork %in% c(0, 1) ~ 0, pork == 2 ~ 1, pork == 3 ~ 2, pork == 4 ~ 5,  pork == 5 ~ 7),
		red_meat = ifelse(beef + lamb + pork <= 2, 1, 0), #no more than twice per week(<= 2 ~ 1)
		meat_processed = case_when(meat_processed %in% c(0, 1) ~ 0, meat_processed == 2 ~ 1, meat_processed == 3 ~ 2, meat_processed == 4 ~ 5,  meat_processed == 5 ~ 7),meat_processed = ifelse(meat_processed < 2, 1, 0), #less than twice per week(< 2 ~ 1)
		diet_score = vegetable + fruit + fish_all + red_meat + meat_processed, diet_score = factor(diet_score, levels = c(0,1,2,3,4,5))) %>%
	mutate( #🏮IPAQ_activity
		IPAQ_activity = factor(IPAQ_activity, levels = c(0, 1, 2), labels = c("low", "moderate", "high"))) %>% 
	mutate(
		medication = ifelse(!is.na(med_male), med_male, ifelse(!is.na(med_female), med_female, NA))) %>%
	mutate( #🏮history_cvd
		#先提取 hf(I50)、af(I48) 的icd10数据
		icdDate_cad2=as.Date(ifelse(!is.na(icdDate_cad), as.character(icdDate_cad), ifelse(!is.na(srdYear_cad), (paste0(srdYear_cad,"-07-01")), ifelse(!is.na(srdAge_cad) & srdAge_cad >0, paste0(birth_year+srdAge_cad,"-07-01"), NA)))),
		history_cad = case_when(heart_dr %in% c(1, 2) ~ 1, !is.na(icdDate_cad2) ~ 1, is.na(heart_dr) & is.na(icdDate_cad2) ~ NA_real_, TRUE ~ 0),
		icdDate_stroke2=as.Date(ifelse(!is.na(icdDate_stroke), as.character(icdDate_stroke), ifelse(!is.na(srdYear_stroke), (paste0(srdYear_stroke,"-07-01")), ifelse(!is.na(srdAge_stroke) & srdAge_stroke >0, paste0(birth_year+srdAge_stroke,"-07-01"), NA)))),
		history_stroke = case_when(heart_dr == 3 ~ 1, !is.na(icdDate_stroke2) ~ 1, is.na(heart_dr) & is.na(icdDate_stroke2) ~ NA_real_, TRUE ~ 0),
		icdDate_hf2=as.Date(ifelse(!is.na(icdDate_hf), as.character(icdDate_hf), ifelse(!is.na(srdYear_hf), (paste0(srdYear_hf,"-07-01")), ifelse(!is.na(srdAge_hf) & srdAge_hf >0, paste0(birth_year+srdAge_hf,"-07-01"), NA)))),				
		history_hf = case_when(!is.na(icdDate_hf2) ~ 1, is.na(icdDate_hf2) ~ 0),
		icdDate_af2=as.Date(ifelse(!is.na(icdDate_af), as.character(icdDate_af), ifelse(!is.na(srdYear_af), (paste0(srdYear_af,"-07-01")), ifelse(!is.na(srdAge_af) & srdAge_af >0, paste0(birth_year+srdAge_af,"-07-01"), NA)))),				
		history_af = case_when(!is.na(icdDate_af2) ~ 1, is.na(icdDate_af2) ~ 0),
		history_cvd = case_when(history_cad == 1 ~ 1, history_stroke == 1 ~ 1, history_hf == 1 ~ 1, history_af == 1 ~ 1, is.na(history_cad) & is.na(history_stroke) & is.na(history_hf) & is.na(history_af) ~ NA_real_, TRUE ~ 0), history_cvd = factor(history_cvd, levels = c(0, 1), labels = c("no", "yes")))%>%
	mutate( #🏮history_dm
		icdDate_dm2=as.Date(ifelse(!is.na(icdDate_dm), as.character(icdDate_dm), ifelse(!is.na(srdYear_dm), (paste0(srdYear_dm,"-07-01")), ifelse(!is.na(srdAge_dm) & srdAge_dm >0, paste0(birth_year+srdAge_dm,"-07-01"), NA)))),
		history_dm = case_when(dm_dr == 1 ~ 1, medication == 3 ~ 1, !is.na(icdDate_dm2) ~ 1, is.na(dm_dr) & is.na(medication) & is.na(icdDate_dm2) ~ NA_real_, TRUE ~ 0), history_dm = factor(history_dm, levels = c(0, 1), labels = c("no", "yes"))) %>%
	mutate( #🏮history_ht
		#先提取 dbp(4079)
		icdDate_ht2=as.Date(ifelse(!is.na(icdDate_ht), as.character(icdDate_ht), ifelse(!is.na(srdYear_ht), (paste0(srdYear_ht,"-07-01")), ifelse(!is.na(srdAge_ht) & srdAge_ht >0, paste0(birth_year+srdAge_ht,"-07-01"), NA)))),
		history_ht = case_when(heart_dr == 4 ~ 1, medication == 2 ~ 1, sbp >= 140 | dbp >= 90 ~ 1, !is.na(icdDate_ht2) ~ 1, is.na(heart_dr) & is.na(medication) & is.na(sbp) & is.na(dbp) & is.na(icdDate_ht2) ~ NA_real_, TRUE ~ 0), history_ht = factor(history_ht, levels = c(0, 1), labels = c("no", "yes")))%>%
	mutate( #🏮leg pain
		walk_chestpain = ifelse(chestpain_walk_normal == 1 | chestpain_walk_hurry == 1, 1, 0),
		walk_legpain = ifelse(legpain_walk_normal == 1 | legpain_walk_hurry == 1, 1, 0),
		leg_surgery = ifelse(leg_arteries == 1 | leg_amputation %in% c(1:3), 1, 0)) %>% rowwise() %>%
	mutate( #🏮sedentary_time
		tv_time = ifelse(tv_time > 24, NA, tv_time),
		computer_time = ifelse(computer_time > 24, NA, computer_time),
		driving_time = ifelse(driving_time > 24, NA, driving_time),
		sedentary_time = sum(c(tv_time, computer_time, driving_time), na.rm = TRUE),
		sedentary_time = ifelse(sedentary_time > 24, 24, sedentary_time)) #exceeded 24h/day, total sedentary time was recoded to 24h/day
saveRDS(dat, file = "D:/data/ukb/Rdata/ukb.life.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 空气和噪音污染，待完成
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 复旦公卫 Predicting particulate matter, nitrogen dioxide, and ozone across Great Britain with high spatiotemporal resolution based on random forest models
air_no2, air_nox, 
noise = rowMeans(cbind(noise_day, noise_evening+5, noise_night+10), na.rm = T),
traffic, dist2road = 1/dist2road, log_traffic = log(traffic), logi_dist2road = -log(dist2road)
