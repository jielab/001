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
	saveRDS(dat0, paste0(indir,"/Rdata/ukb.dash0.rds")) # ⚪

dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans2(select(., starts_with("p100002_")))
)

veg_fields <- c("p104060_|mixveg", "p104070_|vegpiece", "p104080_|Coleslaw", "p104090_|Side", "p104100_|Avocado", "p104130_|Beetroot", "p104140_|Broccoli", "p104150_|Butternut", "p104160_|Cabbage", "p104170_|Carrot", "p104180_|Cauliflower", "p104190_|Celery", "p104200_|Courgette", "p104210_|Cucumber", "p104220_|Garlic", "p104230_|Leek", "p104240_|Lettuce", "p104250_|Mushroom", "p104260_|Onion", "p104270_|Parsnip", "p104290_|pepper", "p104300_|Spinach", "p104310_|Sprouts", "p104320_|Sweetcorn", "p104340_|Freshtomato", "p104350_|Tinnedtomato", "p104360_|Turnip", "p104370_|Watercress", "p104380_|Otherveg")
dat <- bulk_rowMeans(dat, veg_fields) # 蔬菜 🥦 
dat <- dat %>% mutate( # 103990: eat any vegetables yesterday?
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
lapply(dash[-1], function(x) summary(x)) # 🙋‍
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
# lifestyle 南京医科大学版本 💁
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 健康植物性饮食指数（PMID: 36976560）
# 南京医科大学版本（https://github.com/XiangyuYe/Infection_SES）
fields.habit <- c("p189_i0", "p738_i0", "p6138_i0", "p6138-0.1", "p6138-0.2", "p6138-0.3", "p6138-0.4", "p6138-0.5", "p6142_i0", "p6142-0.1", "p6142-0.2", "p6142-0.3", "p6142-0.4", "p6142-0.5", "p6142-0.6", "p20116_i0", "p2897_i0", "p884_i0", "p894_i0", "p904_i0", "p914_i0", "p1309_i0", "p1319_i0", "p1289_i0", "p1299_i0", "p1329_i0", "p1339_i0", "p1349_i0", "p1369_i0", "p1379_i0", "p1389_i0", "p1438_i0", "p1448_i0", "p1458_i0", "p1468_i0", "p20117_i0", "p1180_i0", "p1160_i0", "p1200_i0", "p1210_i0", "p1220_i0", "p20453_i0") 
fields.base <- c("21001_i0", "p22001_i0", "p21003_i0", "p31_i0", "p53_i0", "p53-i1", "p53_i2", "p53_i3", "p21000_i0", "p54_i0") 
fields.all <- c("eid", fields.habit, fields.base)
dat <- read.table(paste0(dir0, '/data/ukb/phe/rap/njmu.tab.gz'), sep="\t", fill=TRUE, header=TRUE)

BMI_grp <- cut(dat$p21001_i0, breaks = c(-Inf, 18.5, 30, Inf),  labels = c(0, 1, 2), right = FALSE)

nosmoke_cnd <- ifelse(is.na(dat$p20116_i0), NA, ifelse(dat$p20116_i0 == 0, 1, ifelse(dat$p20116_i0 > 0, 0, NA)))
	nosmoke_year <- dat$p21003_i0 - dat$p2897_i0
	nosmoke_cnd[which(nosmoke_year > 30)] <- 1

cannabis_cnd <- ifelse(is.na(dat$p20453_i0), NA, ifelse(dat$p20453_i0 == 0, 1, ifelse(dat$p20453_i0 > 0, 0, NA)))

noalcohol_cnd <- ifelse(is.na(dat$p20117_i0), NA,  ifelse(dat$p20117_i0 == 0, 1, ifelse(dat$p20117_i0 > 0, 0, NA)))

num_moderate <- ifelse(dat$p884_i0 < 0, NA, dat$p884_i0)
	num_vigorous <- ifelse(dat$p904_i0 < 0, NA, dat$p904_i0)
	num_activity <- num_moderate > 5 & num_vigorous > 1
	time_moderate <- ifelse(dat$p894_i0 > 150, T, F)
	time_vigorous <- ifelse(dat$p914_i0 > 75, T, F)
	act_mat <- cbind(num_activity, time_moderate, time_vigorous)
	activity_cnd <- rep(0, nrow(act_mat))
	activity_cnd[rowSums(act_mat, na.rm = T) > 0] <- 1
	activity_cnd[rowSums(is.na(act_mat)) == 3] <- NA

diet_numcol <- c(23:36)
	diet_df <- apply(dat[, diet_numcol], 2, function(a) {a[a == -10] <- 0; a[a<0] <- NA; return(a)}) %>% as.data.frame()
	fruit_mat <- cbind(diet_df$p1309_i0, diet_df$p1319_i0/5) 
	fruit_cnd <- rowSums(fruit_mat) >= 4
	fruit_cnd[rowSums(fruit_mat,na.rm = T) >= 4] <- T
	veg_mat <- cbind(diet_df$p1289_i0/3, diet_df$p1299_i0/3)
	veg_cnd <- rowSums(veg_mat) >= 4
	veg_cnd[rowSums(veg_mat,na.rm = T) >= 4] <- T
	fish_mat <- cbind(diet_df$p1329_i0, diet_df$p1339_i0) 
	fish_cnd <- fish_mat[,1] >= 3 | fish_mat[,2] >= 3
	fish_cnd[fish_mat[,1] ==2 & fish_mat[,2] == 2] <- T
	pmeat_cnd <- diet_df$p1349_i0 <= 2
	npmeat_mat <- cbind(diet_df$p1369_i0, diet_df$p1379_i0, diet_df$p1389_i0)
	npmeat_cnd <- rowSums(npmeat_mat, na.rm = T) <= 3
	npmeat_cnd[npmeat_mat[,1] >=3 | npmeat_mat[,2] >=3 | npmeat_mat[,3] >=3] <- F
	grains_mat <- matrix(0, nrow(dat), 2)
	bread_idx <- which(dat$p1448_i0 == 3)
	grains_mat[bread_idx, 1] <- dat$p1438_i0[bread_idx]
	cereal_idx <- which(dat$p1468_i0 %in% c(1, 2, 3))
	grains_mat[cereal_idx, 2] <- dat$p1458_i0[cereal_idx]
	grains_cnd <- rowSums(grains_mat, na.rm = T) >= 3

diet_mat <- data.frame(fruit_cnd, veg_cnd, fish_cnd, pmeat_cnd, npmeat_cnd, grains_cnd)
	diet_cnd <- rowSums(diet_mat) >= 4
	diet_cnd[rowSums(diet_mat,na.rm = T) >= 4] <- T
diet_cnd[rowSums(is.na(diet_mat)) + rowSums(diet_mat,na.rm = T) < 4] <- F

sleep_df <- apply(dat[, c(38:42)], 2, function(a) {a[a<0] <- NA; return(a)}) %>% as.data.frame()
	chrono_cnd <- ifelse(is.na(sleep_df$p1180_i0), NA, ifelse(sleep_df$p1180_i0 < 3, 1, 0) )
	duration_cnd <- ifelse(is.na(sleep_df$p1160_i0), NA, ifelse(sleep_df$p1160_i0 <= 8 & sleep_df$p1160_i0 >= 7, 1, 0) )
	insomnia_cnd <- ifelse(is.na(sleep_df$p1200_i0),NA, ifelse(sleep_df$p1200_i0 == 1, 1, 0))
	snoring_cnd <- ifelse(is.na(sleep_df$p1210_i0),NA, ifelse(sleep_df$p1210_i0 == 2, 1, 0))
	narcolepsy_cnd <- ifelse(is.na(sleep_df$p1220_i0),NA, ifelse(sleep_df$p1220_i0 < 2, 1, 0))
	sleep_mat <- data.frame(chrono_cnd, duration_cnd, insomnia_cnd, snoring_cnd, narcolepsy_cnd)
	sleep_cnd <- rowSums(sleep_mat) >= 4
	sleep_cnd[rowSums(sleep_mat,na.rm = T) >= 4] <- T
	sleep_cnd[rowSums(is.na(sleep_mat)) + rowSums(sleep_mat,na.rm = T) < 4] <- F

lf_df <- data.frame(Nosmoke = nosmoke_cnd, Activity = activity_cnd, Diet = diet_cnd, Noalcohol = noalcohol_cnd, Sleep = sleep_cnd, Nocannabis = cannabis_cnd)
lifescore1 <- rep(NA, nrow(lf_df))
	lifescore1[rowSums(lf_df, na.rm = T) >= 4] <- 3 
	lifescore1[rowSums(lf_df[,c(1:5)], na.rm = T) == 3] <- 3
	lifescore1[rowSums(lf_df, na.rm = T) >= 2 & (rowSums(lf_df, na.rm = T) + rowSums(is.na(lf_df))) < 4 ] <- 2 
	lifescore1[(rowSums(lf_df, na.rm = T) + rowSums(is.na(lf_df))) < 2 | rowSums(lf_df, na.rm = T) == 0] <- 1 
	lifescore2 <- rep(2, nrow(lf_df))
	lifescore2[rowSums(lf_df, na.rm = T) >= 4] <- 3 
	lifescore2[rowSums(lf_df[,c(1:5)], na.rm = T) >= 3] <- 3 
	lifescore2[rowSums(lf_df == 0, na.rm = T) >= 4] <- 1
	lifescore3 <- rep(2, nrow(lf_df))
	lifescore3[rowSums(lf_df, na.rm = T) >= 4] <- 3 
	lifescore3[rowSums(lf_df[,c(1:5)], na.rm = T) >= 3] <- 3 
	lifescore3[rowSums(lf_df == 0, na.rm = T) >= 5] <- 1
	str_split(dat$p53_i0, "-",simplify = T)[,1] ???


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# lifestyle 南京医科大学版本 【修改中】 💁
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- nku %>% 
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
dat <- dat %>% dplyr::select(eid, smoke_score, pa_score, diet_score, alcohol_score, sleep_score, cannabis_score) *1
dat <- dat %>% 
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
saveRDS(life_factor_df, file = "D:/data/ukb/Rdata/ukb.xian.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 社会经济地位（SES）# 💁华中科大
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
