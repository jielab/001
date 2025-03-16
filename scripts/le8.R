pacman::p_load(data.table, tidyverse, lubridate, purrr)
# 源代码（https://github.com/XiangyuYe/Infection_SES）
# 健康植物性饮食指数（PMID: 36976560）

dir0="D:"
indir=paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- read.table(paste0(dir0, '/data/ukb/phe/rap/le8.tab.gz'), sep="\t", fill=TRUE, header=TRUE)
	str(dat0, list.len=800); sapply(dat0, class) 
	dat.nN <- dat0[, !sapply(dat0, is.numeric), drop = FALSE]
	cols.combo <- sapply(dat.nN, function(x) {any(grepl("\\|", as.character(x)), na.rm = TRUE)})
	print(names.combo <- names(cols.combo)[cols.combo])
	# dat0[names.combo] <- lapply(dat0[cols_with_sep], function(x) {as.numeric(sapply(as.character(x), splitMean))})


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DASH🥢
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- dat0 %>% mutate( # 不正常饮食原因 [20085]; Typical diet yesterday [100020]
	across(matches("i[0-4]$") & !matches("p20085_i|p100020_i"),  ~ ifelse(splitMatch(get(paste0("p20085_", sub(".*_(i[0-4])$", "\\1", cur_column()))), c(3,4,6)), NA, .))
)
table(dat0$p20106_i1, useNA="always") # 🚩

names(dat0 %>% select_if(~ any(. %in% c(555, 200, 300, 400, 500, 600)))) # 🚩
	dat0 <- dat0 %>% mutate(across(where(is.numeric), ~ case_when(.==555 ~0.5, .==444 ~0.25, .==200 ~2, .==300 ~3, .==400 ~4, .==500 ~5, .==600 ~6, TRUE ~ .)))
	check <- sapply(dat0, function(x) any(is.nan(x))); names(dat0)[check] # 🚩看是否有NaN 

dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans2(select(., starts_with("p100002_")))
)

veg_fields <- c("p104060_|mixveg", "p104070_|vegpiece", "p104080_|Coleslaw", "p104090_|Side", "p104100_|Avocado", "p104130_|Beetroot", "p104140_|Broccoli", "p104150_|Butternut", "p104160_|Cabbage", "p104170_|Carrot", "p104180_|Cauliflower", "p104190_|Celery", "p104200_|Courgette", "p104210_|Cucumber", "p104220_|Garlic", "p104230_|Leek", "p104240_|Lettuce", "p104250_|Mushroom", "p104260_|Onion", "p104270_|Parsnip", "p104290_|pepper", "p104300_|Spinach", "p104310_|Sprouts", "p104320_|Sweetcorn", "p104340_|Freshtomato", "p104350_|Tinnedtomato", "p104360_|Turnip", "p104370_|Watercress", "p104380_|Otherveg")
dat <- bulk_rowMeans(dat, veg_fields) # 蔬菜 🥦 
dat <- dat %>% mutate( 
	veg = rowSums2(select(., sub(".*\\|", "", veg_fields))),
	veg.new= ifelse((rowSums2(across(starts_with("p103990_i"), ~ !is.na(.))) > 0 & is.na(veg)), 0, veg)
)

fruit_fields <- c("p104410_|Stewedfruit", "p104420_|Prune", "p104430_|Dried", "p104440_|Mixedfruit", "p104450_|Apple", "p104460_|Banana", "p104470_|Berry", "p104480_|Cherry", "p104490_|Grapefruit", "p104500_|Grape", "p104510_|Mango", "p104520_|Melon", "p104530_|Orange", "p104540_|Satsuma", "p104550_|Peach", "p104560_|Pear", "p104570_|Pineapple", "p104580_|Plum", "p104590_|Otherfruit", "p100190_|Orangejuice", "p100200_|Grapefruitjuice", "p100210_|Purefruitjuice")
dat <- bulk_rowMeans(dat, fruit_fields) # 水果🍓
dat <- dat %>% mutate( 
	fruit = rowSums2(select(., sub(".*\\|", "", fruit_fields))),
	fruit.new= ifelse((rowSums2(across(starts_with("p104400_i"), ~ !is.na(.))) > 0 & is.na(fruit)), 0, fruit)
) 

nut_fields <- c("p102410_|Saltpeanut", "p102420_|Unsaltpeanut", "p102430_|Saltnut", "p102440_|Unsaltnut", "p102450_|Seeds", "p103270_|Tofu", "p104000_|Bakedbean", "p104010_|Pulses", "p104110_|Broadbean", "p104120_|Greenbean", "p104280_|Pea")
dat <- bulk_rowMeans(dat, nut_fields) # 坚果🌲
dat <- dat %>% mutate( 
	nut = rowSums2(select(., sub(".*\\|", "", nut_fields))),
	nut.new= ifelse((rowSums2(across(starts_with("p102400_i"), ~ !is.na(.))) > 0 & is.na(nut)), 0, nut)
)

porriage_fields <- c("p100770_|Porridge", "p100800_|Muesli", "p100810_|Oat", "p100820_|Sweetcereal", "p100830_|Plaincereal", "p100840_|Brancereal", "p100850_|Wholewheat", "p100860_|Othercereal", "p102740_|Brownrice", "p102780_|Othergrain", "p102770_|Couscous")
# 主食🍚🍞: 面包类型[20091]; 面包摄入量[100950]; 1:white; 2:mixed; 3:wholemeal; 4:seeded
staples <- c("slicedbread", "baguette", "bap", "breadroll", "wholemealpasta", "crispbread", "oatcakes", "otherbread")
fields <- c("p20091|p100950", "p20092|p101020", "p20093|p101090", "p20094|p101160", "p102720", "p101250", "p101260", "p101270")
for (i in 1:4) {
  field <- strsplit(fields[i], "\\|")[[1]]; field.old <- field[1]; field.new <- field[2]
  dat <- dat %>% mutate(!!staples[i] := rowMeans(across(starts_with(paste0(field.old, "_i")),  ~ ifelse(.==3, get(sub(field.old, field.new, cur_column())), 0)), na.rm=TRUE))
}
for (i in 5:8) { dat <- dat %>% mutate(!!staples[i] := rowMeans2(select(., starts_with(paste0(fields[i], "_"))))) }
dat <- bulk_rowMeans(dat, porriage_fields)
dat <- dat %>% mutate(
	grain = rowSums2(select(., c(sub(".*\\|", "", porriage_fields), all_of(staples)))),
	grain.new= ifelse((rowSums2(across(starts_with("p100760_i"), ~ !is.na(.))) > 0 & is.na(grain)), 0, grain)
)

dat <- dat %>% mutate( # 🐄🥛 酸奶类型[20106]; 酸奶摄入量[102090]
	across(matches("p20106_i[0-4]$"),  ~ ifelse(splitMatch(.x, 210), get(sub("p20106", "p102090", cur_column())), ifelse(splitMatch(.x, 211), 0, NA)), .names = "yogurt_{.col}")
)
dat <- dat %>% mutate( yogurt = rowMeans2(select(., starts_with("yogurt_"))) ) 
dat <- dat %>% mutate(across(paste0("p100920_i", 0:4), ~ ifelse(is.na(.), NA, ifelse(. %in% c(2102, 2103), get(sub("p100920", "p100520", cur_column())), 0)), .names = "milk_{.col}"))
vars <- c("milk|milk_", "hardcheese|p102810_", "cheesespread|p102850_", "cottagecheese|p102870_")
for (var in vars) {
  var_split <- strsplit(var, "\\|")[[1]]; new_var <- var_split[1]; prefix <- var_split[2]
  dat <- dat %>% mutate(!!new_var := rowMeans2(select(., starts_with(prefix))))
}
dat <- dat %>% mutate( 
    cheese = rowSums2(select(., hardcheese, cheesespread, cottagecheese)),
	cheese = ifelse((rowSums2(across(starts_with("p102800_i"), ~ !is.na(.))) > 0 & is.na(cheese)), 0, cheese)
) 
dat <- dat %>% mutate(
    dairy.new = rowSums2(select(., milk, yogurt, cheese))
)

sugar_fields <- c("p100170_|Fizzydrink", "p100180_|Squash", "p100220_|Fruitsmootie") # "p100370_|sugarcoffee", "p100380_|artisugarcoffee", "p100490_|sugartea", "p100500_|artisugartea", "p100540_|Lowhotchocolate", "p100550_|Hotchocolate")
dat <- bulk_rowMeans(dat, sugar_fields) # 🍬🥤
dat <- dat %>% mutate( 
	sugar = rowSums2(select(., sub(".*\\|", "", sugar_fields))),
	sugar.new= ifelse((rowSums2(across(starts_with("p100020_i"), ~ !is.na(.))) > 0 & is.na(sugar)), 0, sugar)
) 

meat_fields <- c("p103010_|Sausage", "p103020_|Beef", "p103030_|Pork", "p103040_|Lamb", "p103070_|Bacon", "p103080_|Ham", "p103090_|Liver", "p102970_|Scotchegg")
dat <- bulk_rowMeans(dat, meat_fields) # 红肉🥩
dat <- dat %>% mutate(
	meat = rowSums2(select(., sub(".*\\|", "", meat_fields))),
	meat.new= ifelse((rowSums2(across(starts_with("p103000_i"), ~ !is.na(.))) > 0 & is.na(meat)), 0, meat)
)
 
dat <- dat %>% mutate( # 盐 ⛵
	sodium.new = p30530_i0,
	sodium.new2 = p26052_i0
)

# 👨‍👩‍👧‍👦
vars.dash <- paste0(c("veg", "fruit", "nut", "grain", "dairy", "sugar", "meat", "sodium"), ".new")
dash <- dat %>% select(eid, all_of(vars.dash))
adjust_breaks <- function(breaks) { # Ensure unique breaks by slightly nudging duplicates
	for (j in 2:length(breaks)) {if (breaks[j] <= breaks[j - 1]) {breaks[j] <- breaks[j-1] + .Machine$double.eps}}
	return(breaks)
}
for (i in vars.dash) {
	print(i)
	print(quin_breaks <- quantile(dash[[i]], probs = seq(0, 1, 0.2), na.rm = TRUE))
	print(quin_breaks <- adjust_breaks(quin_breaks))
	dash[[sub(".new", ".quin", i)]] <- cut(dash[[i]], breaks=quin_breaks, include.lowest=TRUE, labels=1:5, right=TRUE) %>% as.numeric()
}
dash <- dash %>% mutate(across(c(sugar.quin, meat.quin, sodium.quin), ~ ifelse(!is.na(.), 6 - ., NA)))
dash <- dash %>% mutate(dashscore = rowSums2(select(., sub(".new", ".quin", vars.dash))))
dash$dashscore[rowSums(is.na(dash[vars.dash])) > 0] <- NA # ⚡一下让30万人变成NA 🍉
dash <- dash %>% mutate(dash_pts = cut(dashscore, breaks=c(-Inf, 17, 21, 26, 31, Inf), labels=c(0, 25, 50, 80, 100), right=FALSE) %>% as.character() %>% as.numeric())
saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))
lapply(dash[-1], function(x) table(x, useNA = "always")) # 🙋‍


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LE8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds"))
srd <- read.table(paste0(indir,"/rap/srd.tab"), header=T)
srd.drug <- srd %>% mutate( # 💊
	drug.dm=as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1220", "1221", "1222", "1223"))) > 0),
	drug.hpt=as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1065","1072"))) > 0)
	) %>% select(eid, drug.dm, drug.hpt)
le8 <- merge(phe, srd.drug, by="eid") %>% mutate( # drug.cmd[6177] drug.cmd2[6153]
	drug.lipid = ifelse(rowSums2(across(starts_with("drug.cmd"), ~ . %like% "1")) > 0, 1, 0)
)

le8 <- le8 %>% mutate(
	# 运动 🏃‍ MET minutes per week for moderate [22038] or vigorous [22039] activity
	pa_mins = pa_mod + pa_vig * 2,
	pa_pts = cut(pa_mins, breaks = c(0, 1, 30, 60, 90, 120, 150, Inf), labels=c(0, 20, 40, 60, 80, 90, 100), right=FALSE) %>% as.character() %>% as.numeric(),
	# 吸烟 🚬 smoke_quit_age [p2897], smoke_history [p1249], smoke_in_house [p1259]
	smoke_quit_age = ifelse(smoke_quit_age %in% c(-3, -1), NA, smoke_quit_age), 
	smoke_history = ifelse(smoke_history==-3, NA, smoke_history), 
	smoke_in_house = ifelse(smoke_in_house==-3, NA, smoke_in_house),
	smoke_quit_year = ifelse(smoke_quit_age > 0 & age>=smoke_quit_age, age-smoke_quit_age, NA), 
	smoke_cat = case_when(smoke_status==0 ~100, smoke_status==1 & smoke_quit_year>=5 ~75, smoke_history==3 ~75, smoke_status==1 & smoke_quit_year>=1 & smoke_quit_year<5 ~50, smoke_history==2 ~50, smoke_status==1 & smoke_quit_year<1 ~25, smoke_status==2 ~0),
	smoke_pts = ifelse((smoke_in_house==1 | smoke_in_house==2) & smoke_cat>0, smoke_cat-20, smoke_cat),
	smoke_pts = ifelse(!is.na(smoke_in_house) & (smoke_in_house==1 | smoke_in_house==2) & smoke_cat>0, smoke_cat-20, smoke_cat),
	# 睡眠🛏 sleep_duration [p1160]
	sleep_duration = ifelse(sleep_duration %in% c(-1, -3), NA, sleep_duration),
	sleep_pts = cut(sleep_duration, breaks = c(0, 4, 5, 6, 7, 9, 10, 24), labels=c(0, 20, 40, 70, 100, 90, 40), right=FALSE) %>% as.character() %>% as.numeric(),
	# BMI 🎈
	bmi_pts = cut(bmi, breaks=c(0, 25, 30, 35, 40, Inf), labels=c(100, 70, 30, 15, 0), right=FALSE) %>% as.character() %>% as.numeric(),
	# 血脂 🍟
	nonhdl = pmax((bb_TC-bb_HDL)*38.67, 0),
	nonhdl_cat = cut(nonhdl, breaks = c(0, 130, 160, 190, 220, Inf), labels = c(100, 60, 40, 20, 0), right=FALSE) %>% as.character() %>% as.numeric(),
	nonhdl_pts = ifelse(drug.lipid==1 & nonhdl_cat>0, nonhdl_cat-20, nonhdl_cat),
	# 血糖 🍬 dm_dr[2443], dm_gestational[4041]
	hba1c_n = bb_HBA1C*0.0915 + 2.15,
	his.dm = ifelse(drug.dm==1, 1, 0),
	hba1c_pts = case_when(his.dm==0 & hba1c_n>0 & hba1c_n<5.7 ~ 100, his.dm==0 & hba1c_n>=5.7 & hba1c_n<6.5 ~ 60, hba1c_n>0 & hba1c_n<7 ~ 40, hba1c_n>=7 & hba1c_n<8 ~ 30, hba1c_n>=8 & hba1c_n<9 ~ 20, hba1c_n>=9 & hba1c_n<10 ~ 10, hba1c_n>=10 ~ 0),
	# 血压 🦆 sbp_man[93], dbp_man[94], dbp_auto[4079], sbp_auto[4080]
	sbp_auto = rowMeans2(select(., c("sbp_auto.i0.a0", "sbp_auto.i0.a1"))),
	sbp_manual = rowMeans2(select(., c("sbp_man.i0.a0", "sbp_man.i0.a1"))),
	dbp_auto = rowMeans2(select(., c("dbp_auto.i0.a0", "dbp_auto.i0.a1"))),
	dbp_manual = rowMeans2(select(., c("dbp_man.i0.a0", "dbp_man.i0.a1")))
)
le8 <- le8 %>% mutate(
	sbp = rowMeans2(select(., c("sbp_auto", "sbp_manual"))),
	dbp = rowMeans2(select(., c("dbp_auto", "dbp_manual"))),
	bp_cat = case_when(sbp>=160 | dbp>=100 ~ 0, (sbp>=140 & sbp<160) | (dbp>=90 & dbp<100) ~ 25, (sbp>= 130 & sbp< 140) | (dbp>= 80 & dbp< 90) ~ 50, (sbp>=120 & sbp<130) & (dbp>=0 & dbp<80) ~ 75, sbp>0 & sbp<120 & dbp>0 & dbp<80 ~ 100),
	bp_pts=ifelse(drug.hpt==1 & bp_cat>0, bp_cat-20, bp_cat)
)

dash <- readRDS(paste0(indir,"/Rdata/ukb.dash.rds")) # 🥢
le8 <- merge(le8, dash, by="eid")
le8 <- le8 %>% rowwise() %>% # 🎇 汇总
	mutate( LE8score = mean(c(dash_pts, pa_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, bp_pts), na.rm = FALSE)) %>% 
	ungroup() %>% mutate(CVH_cat = cut(LE8score, breaks=c(0, 50, 80, Inf), labels=c(0, 1, 2), right=FALSE) %>% as.character() %>% as.numeric()) %>% 
	dplyr::select(eid, all_of(vars.dash), LE8score, CVH_cat, dash_pts, pa_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, bp_pts)
saveRDS(le8, paste0(indir,"/Rdata/ukb.le8.rds"))
lapply(le8[-(1:9)], function(x) table(x, useNA="always"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 社会经济地位（SES）# 💁
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds"))
dat <- phe %>% filter(ethnic_cat=="White") %>% dplyr::select(eid, age, sex, emp_cat, income_cat, edu_cat) %>% drop_na()
library(poLCA); set.seed(12345); sink("ses.log")
SES_LCA_list <- list()
for (n_class in 2:6) { 
	SES_LCA_list[[n_class]] <- poLCA(cbind(edu_cat, emp_cat, income_cat) ~1, data=dat, nclass=n_class, maxiter=10000, tol=1e-6, nrep=2, graphs=TRUE, probs.start=NULL) 
}
sink()
source("D:/scripts/library/00_LCA_out.R")
LCA_out(SES_LCA_list,3,3)
dat$ses <- SES_LCA_list[[3]]$predclass
life_factor_dff$SES <- SES_LCA_list[[3]]$predclass #
life_factor_dff <- life_factor_dff[,c(1:8,25,9,11,13:22,23,24,10,12)]
saveRDS(subset(dat, select=c(eid, ses)), file="D:/data/ukb/Rdata/ukb.ses.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 生物年龄 (PMID: 37080981) 🈲
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
# 空气和噪音污染，待完成
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
noise = rowMeans(cbind(noise_day, noise_evening+5, noise_night+10), na.rm = T),
traffic, dist2road = 1/dist2road, log_traffic = log(traffic), logi_dist2road = -log(dist2road)
