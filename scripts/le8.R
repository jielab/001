pacman::p_load(data.table, tidyverse, lubridate, purrr)

dir0 = "D:"
indir = paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, '/scripts/f/phe.f.R'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 数据QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fn <- read.table(paste0(indir,"/common/ukb.le8.lst"), sep = "\t", header = FALSE, flush = TRUE)[,1:2]
	fn$V1[duplicated(fn$V1)]; fn$V2[duplicated(fn$V2)] #check duplication
phe0 <- read.delim(paste0(indir,"/rap/le8.tab.gz"), sep = "\t", header = TRUE) 
	phe0 <- phe0 %>% mutate(across(grep("date", names(phe0), value = TRUE), ~as.Date(.x)))
	setdiff(fn$V1, names(phe0)); setdiff(names(phe0), fn$V1)
	names(phe0) <- sapply(names(phe0), replace_code, sep = "_", mapping_df = fn) %>% as.character() # 🏮

#📍处理特殊的负数值
neg_4na <- c(-1, -2, -3, -7, -10, -121, -818)
num_cols <- names(phe0)[vapply(phe0, function(x) is.numeric(x) || inherits(x, "integer64"), logical(1))]
	neg_list <- lapply(phe0[num_cols], function(x) {vals <- sort(unique(x[is.finite(x) & x < 0])); head(vals, 5)}); neg_list[lengths(neg_list) > 0] 
	neg_cols <- names(phe0)[sapply(phe0, function(x) { nv <- unique(x[x < 0 & !is.na(x)]); length(nv) > 0 && isTRUE(all(nv %in% neg_4na)) })]
	phe0 <- phe0 %>% mutate(across(all_of(neg_cols), ~ifelse(.x == -10, 0.5, ifelse(.x < 0, NA, .x)))) 
	# 跑完后再用 neg_list检查一遍🏮

#📍处理特殊的正数值
names(phe0 %>% select_if(~ any(. %in% c(555, 200, 300, 400, 500, 600))))
	phe0 <- phe0 %>% mutate(across(where(is.numeric), ~ case_when(.==555 ~0.5, .==444 ~0.25, .==200 ~2, .==300 ~3, .==400 ~4, .==500 ~5, .==600 ~6, TRUE ~ .)))

# 🚩查找NaN 
	check <- sapply(phe0, function(x) any(is.nan(x))); names(phe0)[check]
		
# 🚩查找复合型变量
phe0 %>% select(where(~ any(str_detect(as.character(.x), fixed("|")), na.rm = TRUE))) %>% names() 
# p6144 不吃 , ⚡p20085 不正常饮食, 🥛type_milk, type_yogurt, 🍞 type_bread, type_baguette, type_bap, type_roll


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LE6 湘雅版本（参考文献 PMID: 40285703)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 地中海饮食 🥢
diet <- phe0 %>% select(matches("^(eid$|med\\.)"))
names(diet) <- gsub("med\\.|_i0$", "", names(diet)); diet[diet < 0] <- NA # ⛸
diet <- diet %>% mutate(
	fruit.sco = ifelse((fruit_fresh + fruit_dried/5) >= 3, 1, 0), # ≥3 servings/day
	veg.sco = ifelse((veg_cooked/3 + veg_raw/3) >= 3, 1, 0), # ≥3 servings/day
	fish.sco = ifelse((fish_oily + fish_other) >= 2, 1, 0), # ≥2 times/week
	meat_processed.sco = ifelse(meat_processed <= 1, 1, 0), # ≤ once a week
	meat_red.sco = ifelse((beef + lamb + pork) <= 2, 1, 0), # ≤2 times/week
	dairy.sco = ifelse(cheese <= 4 & milk_type %in% c(2, 3), 1, 0), # cheese ≤4 times/week
	# oil.sco = ifelse((bread_slice > 2 * 2 * 7 & spread_type %in% c(1,3)), 0, 1),
	oil.sco = ifelse( (((spread_type %in% c(0,2)) | (spread_type == 3 & butter_type %in% c(2, 6, 7, 8))) & bread_slice <= 2 * 2 * 7), 1, 0), # ≤2 servings/day	
	bread_slice.good = ifelse(bread_type == 3, bread_slice, 0),
	bread_slice.bad = ifelse(bread_type %in% c(1, 2, 4), bread_slice, 0),
	cereal_bowl.good = ifelse(cereal_type %in% c(1, 3, 4), cereal_bowl, 0),
	cereal_bowl.bad = ifelse(cereal_type %in% c(2, 5), cereal_bowl, 0),
	grain_whole.sco = ifelse((bread_slice.good +  cereal_bowl.good >= 3 * 7), 1, 0), # ≥3 servings/day
	grain_refined.sco = ifelse((bread_slice.bad +  cereal_bowl.bad <= 2 * 7), 1, 0), # ≤2 servings/day
	sugar.sco = ifelse(never_eat %like% "4", 1, 0) # Don’t drink
) %>% mutate(
	diet.sum = rowSums(across(matches("\\.sco$")), na.rm = TRUE),
	diet.pts = cut(diet.sum, breaks = c(-Inf, 2, 4, 6, 8, Inf), labels = c(0, 25, 50, 80, 100), right = FALSE) %>% as.character() %>% as.numeric(),
)
table(diet$diet.pts)
lapply(diet %>% select(matches("sco$")), table)

# 其他7项 👨‍👩‍👧‍👦
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds")) %>% mutate( 
	drug.lipid = ifelse(rowSums2(across(starts_with("drug.cmd"), ~ . %like% "1")) > 0, 1, 0)
)
srd <- read.table(paste0(indir,"/rap/srd.tab.gz"), header = TRUE) %>% mutate( # 💊
	drug.dm = as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1220", "1221", "1222", "1223"))) > 0),
	drug.hpt = as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1065","1072"))) > 0)
	) %>% select(eid, drug.dm, drug.hpt)
le7 <- merge(phe, srd, by = "eid") 
le7 <- le7 %>% mutate(
	# 运动 🏃‍ 
	pa_mins = dura_pa_mod + dura_pa_vig, # pa_mins = pa_mod + pa_vig * 2 # 北医版本
	pa.pts = cut(pa_mins, breaks = c(0, 1, 30, 60, 90, 120, 150, Inf), labels = c(0, 20, 40, 60, 80, 90, 100), right = FALSE) %>% as.character() %>% as.numeric(),
	# 吸烟 🚬
	smoke_quit_years = ifelse(smoke_quit_age > 0 & age >= smoke_quit_age, age - smoke_quit_age, NA), 
	smoke.pts = ifelse(is.na(smoke_status), NA, ifelse(smoke_status == 0, 100, ifelse(smoke_status == 2, 0, ifelse(smoke_quit_years >=5, 75, ifelse((smoke_quit_years <5 & smoke_quit_years >=1), 50, 25))))),
	smoke.pts = ifelse((smoke.pts > 0 & (smoke_2nd_home >0 | smoke_2nd_outside >0)), smoke.pts - 20, smoke.pts),
	# 睡眠 🛏
	sleep.pts = cut(sleep_duration, breaks = c(0, 4, 5, 6, 7, 9, 10, 24), labels = c(0, 20, 40, 70, 100, 90, 40), right = FALSE) %>% as.character() %>% as.numeric(),
	# BMI 🎈
	bmi.pts = cut(bmi, breaks = c(0, 25, 30, 35, 40, Inf), labels = c(100, 70, 30, 15, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	# 血脂 🍟
	nonhdl = pmax((bb_TC - bb_HDL)*38.67, 0),
	nonhdl_cat = cut(nonhdl, breaks = c(0, 130, 160, 190, 220, Inf), labels = c(100, 60, 40, 20, 0), right = FALSE) %>% as.character() %>% as.numeric(),
	nonhdl.pts = ifelse(drug.lipid == 1 & nonhdl_cat > 0, nonhdl_cat - 20, nonhdl_cat),
	# 血糖 🍬
	hba1c_n = bb_HBA1C * 0.0915 + 2.15,
	his.dm = ifelse(drug.dm == 1, 1, 0),
	hba1c.pts = case_when(his.dm == 0 & hba1c_n > 0 & hba1c_n <5.7 ~ 100, his.dm == 0 & hba1c_n >= 5.7 & hba1c_n < 6.5 ~ 60, hba1c_n > 0 & hba1c_n < 7 ~ 40, hba1c_n >= 7 & hba1c_n < 8 ~ 30, hba1c_n >= 8 & hba1c_n < 9 ~ 20, hba1c_n >= 9 & hba1c_n < 10 ~ 10, hba1c_n >= 10 ~ 0),
	# 血压 🦆 sbp_man[93], dbp_man[94], dbp_auto[4079], sbp_auto[4080]
    sbp = rowMeans(across(matches("^sbp_.*_i0.*")), na.rm = TRUE),
    dbp = rowMeans(across(matches("^dbp_.*_i0.*")), na.rm = TRUE),
	bp.pts = case_when(sbp >= 160 | dbp >= 100 ~ 0, (sbp >= 140 & sbp < 160) | (dbp >= 90 & dbp < 100) ~ 25, (sbp >= 130 & sbp < 140) | (dbp >= 80 & dbp < 90) ~ 50, (sbp >= 120 & sbp < 130) & (dbp >= 0 & dbp < 80) ~ 75, sbp >0 & sbp <120 & dbp >0 & dbp <80 ~ 100),
	bp.pts = ifelse(drug.hpt == 1 & bp.pts > 0, bp.pts - 20, bp.pts)
)
le8 <- merge(le7, diet, by = "eid") %>% select(matches("^eid$|\\.pts$")) %>% mutate( # 🎇 汇总
	le8.sco = rowMeans(across(matches("\\.pts$")), na.rm = TRUE),
	cvh.c = factor(cut(le8.sco, breaks = c(0, 50, 80, Inf), labels = c(0, 1, 2), right = FALSE), levels = 0:2, labels = c("unhealthy", "average", "healthy"))
	)
	lapply(le8 %>% select(matches("\\.pts$")), table)
	saveRDS(le8, paste0(indir,"/Rdata/ukb.le8.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LE8 北医版本（参考文献 PMID： 39265807）
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
vars.le8 <- c("dash", "pa", "smoke", "sleep", "bmi", "nonhdl", "hba1c", "bp")
phe <- phe0 # %>% mutate(energy.sco = rowMeans2(select(., starts_with("energy_"))))

#📍处理不正常饮食特殊情况
phe <- phe %>% mutate( # abnormal饮食原因 [20085]; typical饮食 [100020]
	across(matches("i[0-4]$") & !matches("abnormal_i|typical_i"), ~ ifelse(grepl("(^|\\|)(3|4|6)(\\||$)", 
		paste0(get(paste0("abnormal_", sub(".*_(i[0-4])$", "\\1", cur_column()))), get(paste0("typical_", sub(".*_(i[0-4])$", "\\1", cur_column()))))), NA, .))
)

#📍处理面包和奶制品特殊情况🐄🥛
phe <- phe %>% mutate( # ❗type_bread: 1 white; ❗type_milk: 2102 semiskimmed, 2103 skimmed; ❗type_yogurt: 210 low-fat, 211 full-fat
	across(matches("type_bread_i[0-4]$"),  ~ ifelse(grepl("3", .x), get(sub("type_bread", "bread", cur_column())), 0)),
	across(matches("type_milk_i[0-4]$"),   ~ ifelse(.x %in% c(2102, 2103), get(sub("type_milk", "milk", cur_column())), 0)),
	across(matches("type_yogurt_i[0-4]$"), ~ ifelse(grepl("210", .x), get(sub("type_yogurt", "yogurt", cur_column())), ifelse(grepl("211", .x), 0, NA)))
)

diets <- c("veg", "fruit", "nut", "bread", "porri", "milk", "yogurt", "cheese", "sugar", "meat", "sodium")
  veg_fields	<- c("mixveg", "vegpiece", "coleslaw", "side", "avocado", "beetroot", "broccoli", "butternut", "cabbage", "carrot", "cauliflower", "celery", "courgette", "cucumber", "garlic", "leek", "lettuce", "mushroom", "onion", "parsnip", "pepper", "spinach", "sprouts", "sweetcorn", "freshtomato", "tinnedtomato", "turnip", "watercress", "otherveg")
  fruit_fields	<- c("stewedfruit", "prune", "dried", "mixedfruit", "apple", "banana", "berry", "cherry", "grapefruit", "grape", "mango", "melon", "orange", "satsuma", "peach", "pear", "pineapple", "plum", "otherfruit", "orangejuice", "grapefruitjuice", "purefruitjuice")
  nut_fields	<- c("saltpeanut", "unsaltpeanut", "saltnut", "unsaltnut", "seeds", "tofu", "bakedbean", "pulses", "broadbean", "greenbean", "pea")
  bread_fields 	<- c("bread", "baguette", "bap", "roll", "pasta", "crispbread", "oatcake", "other_bread")
  porri_fields	<- c("porridge", "muesli", "oat", "sweetcereal", "plaincereal", "brancereal", "wholewheat", "othercereal", "brownrice", "othergrain", "couscous")
  milk_fields	<- c("milk")
  yogurt_fields	<- c("yogurt")
  cheese_fields	<- c("hardcheese", "cheesespread", "cottagecheese")
  sugar_fields	<- c("fizzydrink", "squash", "fruitsmootie", "sugarcoffee", "artisugarcoffee", "sugartea", "artisugartea", "lowhotchocolate", "hotchocolate")
  meat_fields	<- c("sausage", "beef", "pork", "lamb", "bacon", "ham", "liver", "scotchegg")
  sodium_fields	<- c("sodium", "sodium_urine")

for (diet in diets) { # .mi 表示多个 *_i 的平均值mean
	phe <- bulk_rowMeans(phe, get(paste0(diet, "_fields"))) #下面第二个paste0之前必须要用get🏮
	phe <- phe %>% mutate(!!paste0(diet, ".sco") := rowSums2(select(., paste0(get(paste0(diet, "_fields")), ".mi"))))
}

# 把 🐄🥛三合一
grep("\\.sco$", names(phe), value = TRUE)
	vars.dairy = c("milk", "yogurt", "cheese")
	phe <- phe %>% mutate(dairy.sco = rowSums2(select(., paste0(vars.dairy, ".sco"))))
	diets <- append(diets[!(diets %in% vars.dairy)], "dairy", after=5)

# 📍处理有_overall变量的特殊情况
# 原代码: meat.new= ifelse(rowSums2(across(starts_with("p103000_i"), ~ !is.na(.))) > 0 & is.na(meat)), 0, meat) 
vars <- c("veg", "fruit", "nut", "meat") # ❓"cereal", "cheese" 
for (f in vars) {
	f.sco <- paste0(f, ".sco"); f.1 <- grep(paste0("^", f, "_overall"), names(phe))
	phe[[f.sco]] <- ifelse(rowSums(!is.na(phe[f.1])) > 0 &is.na(phe[[f.sco]]), 0, phe[[f.sco]])
}

vars.dash <- paste0(diets, ".sco")
	dash <- phe %>% select(eid, all_of(vars.dash))
	lapply(dash, function(x) summary(x)) # 🙋‍

# 📍处理特殊情况
adjust_breaks <- function(breaks) { # Ensure unique breaks by slightly nudging duplicates
	for (j in 2:length(breaks)) {if (breaks[j] <= breaks[j - 1]) {breaks[j] <- breaks[j-1] + .Machine$double.eps}}
	return(breaks)
}
for (i in vars.dash) {
	print(i)
	print(quin_breaks <- quantile(dash[[i]], probs = seq(0, 1, 0.2), na.rm = TRUE))
	print(quin_breaks <- adjust_breaks(quin_breaks))
	dash[[sub(".sco", ".quin", i)]] <- cut(dash[[i]], breaks = quin_breaks, include.lowest = TRUE, labels = 1:5, right = TRUE) %>% as.numeric()
}
dash <- dash %>% mutate(across(c(sugar.quin, meat.quin, sodium.quin), ~ ifelse(!is.na(.), 6 - ., NA)))
	dash <- dash %>% mutate(dash.sco = rowSums2(select(., sub(".sco", ".quin", vars.dash))))
	dash$dash.sco[rowSums(is.na(dash[vars.dash])) > 0] <- NA # ⚡一下让30万人变成NA 🍉
	dash <- dash %>% mutate(dash.pts = cut(dash.sco, breaks = c(-Inf, 17, 21, 26, 31, Inf), labels = c(0, 25, 50, 80, 100), right = FALSE) %>% 
		as.character() %>% as.numeric()) # 🏮 如果省略掉as.character()，就会得到 1, 2, 3, 4, 5，而不是 0, 25, 50, 80, 100
	lapply(dash[2:9], function(x) summary(x))
	
saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))

# 其他7项 👨‍👩‍👧‍👦
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds"))
srd <- read.table(paste0(indir,"/rap/srd.tab.gz"), header=T)
srd <- srd %>% mutate( # 💊
	drug.dm=as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1220", "1221", "1222", "1223"))) > 0),
	drug.hpt=as.integer(rowSums2(across(starts_with("p20002_i0_"),  ~ . %in% c("1065","1072"))) > 0)
	) %>% select(eid, drug.dm, drug.hpt)
le7 <- merge(phe, srd, by="eid") %>% mutate( # drug.cmd[6177] drug.cmd2[6153]
	drug.lipid = ifelse(rowSums2(across(starts_with("drug.cmd"), ~ . %like% "1")) > 0, 1, 0)
)
le7 <- le7 %>% mutate(
	# 运动 🏃‍ MET minutes per week for moderate [22038] or vigorous [22039] activity
	pa_mins = pa_moderate + pa_vigorous * 2,
	pa.pts = cut(pa_mins, breaks = c(0, 1, 30, 60, 90, 120, 150, Inf), labels=c(0, 20, 40, 60, 80, 90, 100), right=FALSE) %>% as.character() %>% as.numeric(),
	# 吸烟 🚬 smoke_quit_age [p2897], smoke_past [p1249], smoke_at_home [p1259]
	smoke_quit_age = ifelse(smoke_quit_age %in% c(-3, -1), NA, smoke_quit_age), 
	smoke_quit_year = ifelse(smoke_quit_age > 0 & age >= smoke_quit_age, age - smoke_quit_age, NA), 
	smoke_cat = case_when(smoke_status==0 ~100, smoke_status==1 & smoke_quit_year>=5 ~75, smoke_past==3 ~75, smoke_status==1 & smoke_quit_year>=1 & smoke_quit_year<5 ~50, smoke_past==2 ~50, smoke_status==1 & smoke_quit_year<1 ~25, smoke_status==2 ~0),
	smoke.pts = ifelse(smoke_at_home %in% 1:2 & smoke_cat > 0, smoke_cat - 20, smoke_cat),
	# 睡眠🛏 sleep_duration [p1160]
	sleep_duration = ifelse(sleep_duration %in% c(-1, -3), NA, sleep_duration),
	sleep.pts = cut(sleep_duration, breaks = c(0, 4, 5, 6, 7, 9, 10, 24), labels=c(0, 20, 40, 70, 100, 90, 40), right=FALSE) %>% as.character() %>% as.numeric(),
	# BMI 🎈
	bmi.pts = cut(bmi, breaks=c(0, 25, 30, 35, 40, Inf), labels=c(100, 70, 30, 15, 0), right=FALSE) %>% as.character() %>% as.numeric(),
	# 血脂 🍟
	nonhdl = pmax((bb_TC-bb_HDL)*38.67, 0),
	nonhdl_cat = cut(nonhdl, breaks = c(0, 130, 160, 190, 220, Inf), labels = c(100, 60, 40, 20, 0), right=FALSE) %>% as.character() %>% as.numeric(),
	nonhdl.pts = ifelse(drug.lipid==1 & nonhdl_cat>0, nonhdl_cat-20, nonhdl_cat),
	# 血糖 🍬 dm_dr[2443], dm_gestational[4041]
	hba1c_n = bb_HBA1C*0.0915 + 2.15,
	his.dm = ifelse(drug.dm==1, 1, 0),
	hba1c.pts = case_when(his.dm==0 & hba1c_n>0 & hba1c_n<5.7 ~ 100, his.dm==0 & hba1c_n>=5.7 & hba1c_n<6.5 ~ 60, hba1c_n>0 & hba1c_n<7 ~ 40, hba1c_n>=7 & hba1c_n<8 ~ 30, hba1c_n>=8 & hba1c_n<9 ~ 20, hba1c_n>=9 & hba1c_n<10 ~ 10, hba1c_n>=10 ~ 0),
	# 血压 🦆 sbp_man[93], dbp_man[94], dbp_auto[4079], sbp_auto[4080]
	sbp_auto = rowMeans2(select(., c("sbp_auto_i0_a0", "sbp_auto_i0_a1"))),
	sbp_manual = rowMeans2(select(., c("sbp_man_i0_a0", "sbp_man_i0_a1"))),
	dbp_auto = rowMeans2(select(., c("dbp_auto_i0_a0", "dbp_auto_i0_a1"))),
	dbp_manual = rowMeans2(select(., c("dbp_man_i0_a0", "dbp_man_i0_a1")))
)
le7 <- le7 %>% mutate(
	sbp = rowMeans2(select(., c("sbp_auto", "sbp_manual"))),
	dbp = rowMeans2(select(., c("dbp_auto", "dbp_manual"))),
	bp_cat = case_when(sbp>=160 | dbp>=100 ~ 0, (sbp>=140 & sbp<160) | (dbp>=90 & dbp<100) ~ 25, (sbp>= 130 & sbp< 140) | (dbp>= 80 & dbp< 90) ~ 50, (sbp>=120 & sbp<130) & (dbp>=0 & dbp<80) ~ 75, sbp>0 & sbp<120 & dbp>0 & dbp<80 ~ 100),
	bp.pts=ifelse(drug.hpt==1 & bp_cat>0, bp_cat-20, bp_cat)
)
saveRDS(le7, paste0(indir,"/Rdata/ukb.pku.le7.rds"))
 
le8 <- le7 %>% merge(dash, by="eid") %>% select(eid, all_of(paste0(vars.le8, ".pts"))) %>% mutate(
	le8.sco = rowMeans(across(all_of(paste0(vars.le8, ".pts"))), na.rm = TRUE),
	cvh.c = cut(le8.sco, breaks=c(0, 50, 80, Inf), labels=c(0, 1, 2), right=FALSE)
) 
	
lapply(le8[grep("\\.pts", names(le8))], function(x) table(x, useNA="always"))
names(le8)[-1] <- gsub("$", ".pku", names(le8)[-1])
saveRDS(le8, paste0(indir,"/Rdata/ukb.pku.le8.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LE6 南医版本（PMID: 36717939）
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://github.com/XiangyuYe/Infection_SES
dat0 <- read.table(paste0(dir0, '/data/ukb/phe/rap/njmu.tab.gz'), sep = "\t", fill = TRUE, header = TRUE)
dat <- dat0 %>% mutate(across(everything(), ~ ifelse(. < 0, NA, .)))

dat <- dat %>% mutate( #🚬🍺
	nosmoke_sco = ifelse(p20116_i0 == 0, 1, ifelse(p20116_i0 > 0, 0, NA)),
	nosmoke_year = p21003_i0 - p2897_i0,
	nosmoke_sco[which(nosmoke_year > 30)] <- 1,
	cannabis_sco = ifelse(p20453 == 0, 1, ifelse(p20453 > 0, 0, NA)),
	noalcohol_sco = ifelse(p20117_i0 == 0, 1, ifelse(p20117_i0 > 0, 0, NA))
)

dat <- dat %>% mutate( # 🥦
	fruit_sco = rowSums(cbind(p1309_i0, p1319_i0 / 5), na.rm = TRUE) >= 4,
	veg_sco = rowSums(cbind(p1289_i0 / 3, p1299_i0 / 3), na.rm = TRUE) >= 4,
	fish_sco = (p1329_i0 >= 3 | p1339_i0 >= 3) | (p1329_i0 == 2 & p1339_i0 == 2),
	pmeat_sco = p1349_i0 <= 2,
	npmeat_sco = rowSums(cbind(p1369_i0, p1379_i0, p1389_i0), na.rm = TRUE) <= 3 & !(p1369_i0 >= 3 | p1379_i0 >= 3 | p1389_i0 >= 3),
	grains_sco = rowSums(cbind(ifelse(p1448_i0 == 3, p1438_i0, 0), ifelse(p1468_i0 %in% c(1, 2, 3), p1458_i0, 0)), na.rm = TRUE) >= 3
)
diet_fileds <- c("fruit_sco", "veg_sco", "fish_sco", "pmeat_sco", "npmeat_sco", "grains_sco")
dat <- dat %>% mutate(
  diet_sum = rowSums(select(., all_of(diet_fileds)), na.rm = TRUE),
  diet_sco = ifelse(diet_sum >= 4, TRUE, ifelse(rowSums(is.na(select(., all_of(diet_fileds)))) + diet_sum < 4, FALSE, NA))
) 

dat <- dat %>% mutate( # 🏃‍
	num_activity = p884_i0 > 5 & p904_i0 > 1,
	time_moderate = p894_i0 > 150,
	time_vigorous = p914_i0 > 75,
	mat_activity = rowSums(cbind(num_activity, time_moderate, time_vigorous), na.rm = TRUE),
	activity_sco = ifelse(mat_activity > 0, 1, ifelse(mat_activity == 3, NA, 0))
)

dat <- dat %>% mutate( # 🛏
	chrono_sco = ifelse(!is.na(p1180_i0) & p1180_i0 < 3, 1, 0),
	duration_sco = ifelse(!is.na(p1160_i0) & p1160_i0 >= 7 & p1160_i0 <= 8, 1, 0),
	insomnia_sco = ifelse(!is.na(p1200_i0) & p1200_i0 == 1, 1, 0),
	snoring_sco = ifelse(!is.na(p1210_i0) & p1210_i0 == 2, 1, 0),
	narcolepsy_sco = ifelse(!is.na(p1220_i0) & p1220_i0 < 2, 1, 0)
) 
sleep_mat <- select(dat, chrono_sco, duration_sco, insomnia_sco, snoring_sco, narcolepsy_sco)
dat <- dat %>% mutate(
	sleep_sum = rowSums(sleep_mat, na.rm = TRUE),
	sleep_sco = ifelse(sleep_sum >= 4, TRUE, ifelse(rowSums(is.na(sleep_mat)) + sleep_sum < 4, FALSE, NA))
)

le_fileds <- c("nosmoke_sco", "cannabis_sco", "noalcohol_sco", "diet_sco", "activity_sco", "sleep_sco")
dat <- dat %>% mutate( # 🎇 汇总
	le_sum = rowSums(select(., all_of(le_fileds)), na.rm = TRUE),
	missing_cnt = rowSums(is.na(select(., all_of(le_fileds)))),
	le6 = ifelse(le_sum >= 4 | le_sum == 3, 3, ifelse(le_sum >= 2 & (le_sum + missing_cnt) < 4, 2, ifelse((le_sum + missing_cnt) < 2 | le_sum == 0, 1, NA)))
)
saveRDS(dat, paste0(indir,"/Rdata/ukb.nku.le6.rds"))
