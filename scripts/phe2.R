setwd("D:/")
pacman::p_load(data.table, tidyverse, lubridate, poLCA, BioAge)

phe0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")


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
saveRDS(life_factor_df, file = "D:/data/ukb/Rdata/ukb.life.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# walk-VTE文章版本
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- phe0 %>% 
	mutate( #🏮diet_score
		veg_cooked = ifelse(veg_cooked < 0 | veg_cooked == 50, NA, veg_cooked),
		veg_salad = ifelse(veg_salad < 0 | veg_salad == 50, NA, veg_salad),
		vegetable = ifelse(veg_cooked + veg_salad >= 4, 1, 0),	#at least 4 tablespoons/day
		fruit_fresh = ifelse(fruit_fresh < 0 | fruit_fresh == 50, NA, fruit_fresh),
		fruit_dried = ifelse(fruit_dried < 0 | fruit_dried == 100, NA, fruit_dried), fruit = ifelse(fruit_fresh + fruit_dried >= 3, 1, 0), #at least 3 pieces/day
		fish_oily = case_when(fish_oily < 0 ~ NA_real_, fish_oily %in% c(0, 1) ~ 0, fish_oily == 2 ~ 1, fish_oily == 3 ~ 2, fish_oily == 4 ~ 5,  fish_oily == 5 ~ 7),
		fish = case_when(fish < 0 ~ NA_real_, fish %in% c(0, 1) ~ 0, fish == 2 ~ 1, fish == 3 ~ 2, fish == 4 ~ 5,  fish == 5 ~ 7),
		fish = ifelse(fish_oily + fish >= 2, 1, 0), #at least twice per week(>= 2 ~ 1)	
		beef = case_when(beef < 0 ~ NA_real_, beef %in% c(0, 1) ~ 0, beef == 2 ~ 1, beef == 3 ~ 2, beef == 4 ~ 5,  beef == 5 ~ 7),
		lamb = case_when(lamb < 0 ~ NA_real_, lamb %in% c(0, 1) ~ 0, lamb == 2 ~ 1, lamb == 3 ~ 2, lamb == 4 ~ 5,  lamb == 5 ~ 7),
		pork = case_when(pork < 0 ~ NA_real_, pork %in% c(0, 1) ~ 0, pork == 2 ~ 1, pork == 3 ~ 2, pork == 4 ~ 5,  pork == 5 ~ 7),
		red_meat = ifelse(beef + lamb + pork <= 2, 1, 0), #no more than twice per week(<= 2 ~ 1)
		meat_processed = case_when(meat_processed < 0 ~ NA_real_, meat_processed %in% c(0, 1) ~ 0, meat_processed == 2 ~ 1, meat_processed == 3 ~ 2, meat_processed == 4 ~ 5,  meat_processed == 5 ~ 7),meat_processed = ifelse(meat_processed < 2, 1, 0), #less than twice per week(< 2 ~ 1)
		diet_score = vegetable + fruit + fish + red_meat + meat_processed, diet_score = factor(diet_score, levels = c(0,1,2,3,4,5))) %>%
	mutate( #🏮IPAQ_activity
		IPAQ_activity = factor(IPAQ_activity, levels = c(0, 1, 2), labels = c("low", "moderate", "high"))) %>% rowwise() %>%
	mutate( #🏮sedentary_time
		tv_time = ifelse(tv_time < 0 | tv_time > 24, NA, tv_time),
		computer_time = ifelse(computer_time < 0 | computer_time > 24, NA, computer_time),
		driving_time = ifelse(driving_time < 0 | driving_time > 24, NA, driving_time),
		sedentary_time = sum(c(tv_time, computer_time, driving_time), na.rm = TRUE),
		sedentary_time = ifelse(sedentary_time > 24, 24, sedentary_time)) %>% #exceeded 24h/day, total sedentary time was recoded to 24h/day
	mutate(
		med_male = ifelse(med_male < 0, NA, med_male),
		med_female = ifelse(med_female < 0, NA, med_female),
		medication = ifelse(!is.na(med_male), med_male, ifelse(!is.na(med_female), med_female, NA)),
		cholesterol_lowering_medication = ifelse(medication == 1, 1, ifelse(is.na(medication), NA, 0)),
		ht_medication = ifelse(medication == 2, 1, ifelse(is.na(medication), NA, 0)),
		diabetes_medication = ifelse(medication == 3, 1, ifelse(is.na(medication), NA, 0))) %>%
	mutate( #🏮history_cvd
		non_cancer = ifelse(non_cancer < 0 | non_cancer == 99999, NA, non_cancer),
		heart_dr = ifelse(heart_dr < 0, NA, heart_dr),
		cad_yes = ifelse(is.na(icdDate_cad), 0, 1),
		history_cad = case_when(non_cancer %in% c(1074, 1075) ~ 1, heart_dr %in% c(1, 2) ~ 1,cad_yes == 1 ~ 1,is.na(non_cancer) & is.na(heart_dr) & is.na(cad_yes) ~ NA_real_, TRUE ~ 0),
		stroke_yes = ifelse(is.na(icdDate_stroke), 0, 1),
		history_stroke = case_when(non_cancer %in% c(1081, 1082, 1086, 1491) ~ 1,heart_dr == 3 ~ 1,stroke_yes == 1 ~ 1,is.na(non_cancer) & is.na(heart_dr) & is.na(stroke_yes) ~ NA_real_, TRUE ~ 0),
		hf_yes = ifelse(is.na(icdDate_hf), 0, 1),
		history_hf = case_when(non_cancer == 1076 ~ 1,hf_yes == 1 ~ 1,is.na(non_cancer) & is.na(hf_yes) ~ NA_real_, TRUE ~ 0),
		af_yes = ifelse(is.na(icdDate_af), 0, 1),
		history_af = case_when(non_cancer %in% c(1471, 1483, 1583) ~ 1, af_yes == 1 ~ 1, is.na(non_cancer) & is.na(af_yes) ~ NA_real_, TRUE ~ 0),
		history_cvd = ifelse(history_cad == 1 | history_stroke == 1 | history_hf == 1 | history_af == 1, 1, 0),
		history_cvd = factor(history_cvd, levels = c(0, 1), labels = c("no", "yes")),)%>%
	mutate( #🏮history_dm 
		dm_yes = ifelse(is.na(icdDate_dm), 0, 1),
		non_cancer = ifelse(non_cancer < 0 | non_cancer == 99999, NA, non_cancer),
		dm_dr = ifelse(dm_dr < 0, NA, dm_dr),
		history_diabetes = case_when(non_cancer %in% c(1220, 1222, 1223) ~ 1, dm_dr == 1 ~ 1, medication == 3 ~ 1, dm_yes == 1 ~ 1, is.na(non_cancer) & is.na(dm_dr) & is.na(medication) & is.na(dm_yes) ~ NA_real_, TRUE ~ 0 ),
		history_diabetes = factor(history_diabetes, levels = c(0, 1), labels = c("no", "yes")),) %>%
	mutate( #🏮history_ht
		ht_yes = ifelse(is.na(icdDate_ht), 0, 1),
		non_cancer = ifelse(non_cancer < 0 | non_cancer == 99999, NA, non_cancer),
		heart_dr = ifelse(heart_dr < 0, NA, heart_dr),
		history_ht = case_when(non_cancer == 1065 ~ 1,heart_dr == 4 ~ 1, medication == 2 ~ 1, sbp >= 140 | Diastolic_blood_pressure >= 90 ~ 1, ht_yes == 1 ~ 1, is.na(non_cancer) & is.na(heart_dr) & is.na(medication) & is.na(sbp) & is.na(Diastolic_blood_pressure) & is.na(ht_yes) ~ NA_real_, TRUE ~ 0),
		history_ht = factor(history_ht, levels = c(0, 1), labels = c("no", "yes")))%>%
	mutate( #🏮leg pain
		chestpain_walk_normal = ifelse(chestpain_walk_normal < 0, NA, chestpain_walk_normal),chestpain_walk_hurry = ifelse(chestpain_walk_hurry < 0, NA, chestpain_walk_hurry),chestpain_walk = ifelse(chestpain_walk_normal == 1 | chestpain_walk_hurry == 1, 1, 0),
		legpain_walk_normal = ifelse(legpain_walk_normal < 0, NA, legpain_walk_normal),legpain_walk_hurry = ifelse(legpain_walk_hurry < 0, NA, legpain_walk_hurry),legpain_walk = ifelse(legpain_walk_normal == 1 | legpain_walk_hurry == 1, 1, 0),
		leg_arteries = ifelse(leg_arteries < 0, NA, leg_arteries), leg_amputation = ifelse(leg_amputation < 0, NA, leg_amputation),
		leg_surgery = ifelse(leg_arteries == 1 | leg_amputation %in% c(1:3), 1, 0)
	)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 空气和噪音污染，待完成
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 复旦公卫 Predicting particulate matter, nitrogen dioxide, and ozone across Great Britain with high spatiotemporal resolution based on random forest models
air_no2, air_nox, 
noise = rowMeans(cbind(noise_day, noise_evening+5, noise_night+10), na.rm = T),
traffic, dist2road = 1/dist2road, log_traffic = log(traffic), logi_dist2road = -log(dist2road)
