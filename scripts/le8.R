pacman::p_load(tidyverse, lubridate, purrr)

dir0="D:"
indir=paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, "/scripts/ukb/le8.dash.fields"))
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- read.table(paste0(dir0, '/data/ukb/phe/rap/raw/le8.pku.raw.gz'), sep="\t", fill=TRUE, header=TRUE) 
dat0 <- dat0 %>% mutate(
	across(where(is.numeric), ~ case_when(. == 555 ~ 0.5, . == 444 ~ 0.25, . == 200 ~ 2, . == 300 ~ 3, . == 400 ~ 4, . == 500 ~ 5, . == 600 ~ 6, TRUE ~ .))
)
# str(dat0, list.len=800); sapply(dat0, class);  

dat0 <- dat0 %>% mutate(
  across(.cols = matches("i0$") & !matches("p20085_i0") & !matches("p100020_i0"), .fns = ~ ifelse(p20085_i0 %in% c(3, 4, 6), NA, .)),
  across(.cols = matches("i1$") & !matches("p20085_i1") & !matches("p100020_i1"), .fns = ~ ifelse(p20085_i1 %in% c(3, 4, 6), NA, .)),
  across(.cols = matches("i2$") & !matches("p20085_i2") & !matches("p100020_i2"), .fns = ~ ifelse(p20085_i2 %in% c(3, 4, 6), NA, .)),
  across(.cols = matches("i3$") & !matches("p20085_i3") & !matches("p100020_i3"), .fns = ~ ifelse(p20085_i3 %in% c(3, 4, 6), NA, .)),
  across(.cols = matches("i4$") & !matches("p20085_i4") & !matches("p100020_i4"), .fns = ~ ifelse(p20085_i4 %in% c(3, 4, 6), NA, .))
)
# lapply(dat0, function(x) any(is.nan(x))) # 🏮看是否有NaN 


dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans2(select(., starts_with("p100002_")))
)

dat <- bulk_rowMeans(dat, veg_fields_names) # 蔬菜 🥦 
dat <- dat %>% mutate( 
	vegetable = rowSums2(select(., sub(".*\\|", "", veg_fields_names))),
	vegetablenew= ifelse((rowSums(across(starts_with("p103990_i"), ~ !is.na(.))) > 0 & is.na(vegetable)), 0, vegetable)
)

dat <- bulk_rowMeans(dat, fruit_fields_names) # 水果🍓
dat <- dat %>% mutate( 
	fruit = rowSums2(select(., sub(".*\\|", "", fruit_fields_names))),
	fruitnew= ifelse((rowSums(across(starts_with("p104400_i"), ~ !is.na(.))) > 0 & is.na(fruit)), 0, fruit)
) 

dat <- bulk_rowMeans(dat, nut_fields_names) # 坚果🌲
dat <- dat %>% mutate( 
	nut = rowSums2(select(., sub(".*\\|", "", nut_fields_names))),
	nutnew= ifelse((rowSums(across(starts_with("p102400_i"), ~ !is.na(.))) > 0 & is.na(nut)), 0, nut)
)

dat <- bulk_rowMeans(dat, porriage_fields_names) # 主食🍚🍞
dat <- dat %>% mutate( 
    Slicedbread = rowMeans2(replace_if_equal(dat, slicedbread_fields, 3)),
    Baguette = rowMeans2(replace_if_equal(dat, baguette_fields, 3)),
    Bap = rowMeans2(replace_if_equal(dat, bap_fields, 3)),
    Breadroll = rowMeans2(replace_if_equal(dat, breadroll_fields, 3)),
    Wholemealpasta = rowMeans2(select(., starts_with("p102720_"))),
    Crispbread = rowMeans2(select(., starts_with("p101250_"))),
    Oatcakes = rowMeans2(select(., starts_with("p101260_"))),
    Otherbread = rowMeans2(select(., starts_with("p101270_")))
)
dat <- dat %>% mutate(
	grain = rowSums2(select(., c(sub(".*\\|", "", porriage_fields_names), "Slicedbread", "Baguette", "Bap", "Breadroll", "Wholemealpasta", "Crispbread", "Oatcakes", "Otherbread"))),
	grainnew= ifelse((rowSums(across(starts_with("p100760_i"), ~ !is.na(.)), na.rm=TRUE) > 0 & is.na(grain)), 0, grain)
)

dat <- dat %>% mutate( # 🐄🥛
    yogurt_0 = ifelse(p20106_i0 == 210, p102090_i0, ifelse(p20106_i0 == 211, 0, NA)),
    yogurt_1 = ifelse(p20106_i1 == 210, p102090_i1, ifelse(p20106_i1 == 211, 0, NA)),
    yogurt_2 = ifelse(p20106_i2 == 210, p102090_i2, ifelse(p20106_i2 == 211, 0, NA)),
    yogurt_3 = ifelse(p20106_i3 == 210, p102090_i3, ifelse(p20106_i3 == 211, 0, NA)),
    yogurt_4 = ifelse(p20106_i4 == 210, p102090_i4, ifelse(p20106_i4 == 211, 0, NA))
)
dat <- dat %>% mutate(
    Yogurt = rowMeans2(select(., starts_with("yogurt_"))),
    milk_0 = ifelse(p100920_i0 %in% c(2102, 2103), p100520_i0, 0),
    milk_1 = ifelse(p100920_i1 %in% c(2102, 2103), p100520_i1, 0),
    milk_2 = ifelse(p100920_i2 %in% c(2102, 2103), p100520_i2, 0),
    milk_3 = ifelse(p100920_i3 %in% c(2102, 2103), p100520_i3, 0),
    milk_4 = ifelse(p100920_i4 %in% c(2102, 2103), p100520_i4, 0),
    milk_0 = ifelse(is.na(p100920_i0), NA, milk_0),
    milk_1 = ifelse(is.na(p100920_i1), NA, milk_1),
    milk_2 = ifelse(is.na(p100920_i2), NA, milk_2),
    milk_3 = ifelse(is.na(p100920_i3), NA, milk_3),
    milk_4 = ifelse(is.na(p100920_i4), NA, milk_4)
)
dat <- dat %>% mutate(
    Milk = rowMeans2(select(., starts_with("milk_"))),
    hardcheese = rowMeans2(select(., starts_with("p102810_"))),
    cheesespread = rowMeans2(select(., starts_with("p102850_"))),
    cottagecheese = rowMeans2(select(., starts_with("p102870_")))
)
dat <- dat %>% mutate( 
    cheese = rowSums2(select(., hardcheese, cheesespread, cottagecheese)),
	cheesenew= ifelse((rowSums(across(starts_with("p102800_i"), ~ !is.na(.))) > 0 & is.na(cheese)), 0, cheese)
) 
dat <- dat %>% mutate(
    lowfatdairy = rowSums2(select(., Milk, Yogurt, cheesenew))
)

dat <- bulk_rowMeans(dat, sugar_fields_names) # 🍬🥤
dat <- dat %>% mutate( 
	sugar = rowSums2(select(., sub(".*\\|", "", sugar_fields_names)))
) 
dat <- dat %>% mutate(
#	sugarnew= ifelse(any(!is.na(c("p100020_i0", "p100020_i1", "p100020_i2", "p100020_i3", "p100020_i4")) && is.na(sugarnew)), 0, sugar) 
	sugarnew= ifelse((rowSums2(across(starts_with("p100020_i"), ~ !is.na(.))) > 0 & is.na(sugar)), 0, sugar)
) 

dat <- bulk_rowMeans(dat, meat_fields_names) # 红肉🥩
dat <- dat %>% mutate(
	meat = rowSums2(select(., sub(".*\\|", "", meat_fields_names))),
	meatnew= ifelse((rowSums(across(starts_with("p103000_i"), ~ !is.na(.))) > 0 & is.na(meat)), 0, meat)
)
 
dat <- dat %>% mutate( # 盐 ⛵
	sodium = p30530_i0
)

dash <- dat %>% select(eid, vegetablenew, fruitnew, nutnew, grainnew, lowfatdairy, sugarnew, meatnew, sodium)
saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.semi.rds"))
write.table(dash, "ukb.le8.dash.R.tsv", append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)

QUIN <- c("vegetablenew", "fruitnew", "nutnew", "grainnew", "lowfatdairy", "sugarnew", "meatnew", "sodium")
adjust_breaks <- function(breaks) { # Ensure unique breaks by slightly nudging duplicates
  for (j in 2:length(breaks)) {if (breaks[j] <= breaks[j - 1]) {breaks[j] <- breaks[j - 1] + .Machine$double.eps}}
  return(breaks)
}
for (i in QUIN) {
	quin_boundaries <- quantile(dash[[i]], probs = seq(0, 1, 0.2), na.rm = TRUE)
	quin_boundaries <- adjust_breaks(quin_boundaries)
	new_col_name <- paste0("quin", sub("new", "", i))
	dash <- dash %>% mutate(!!sym(new_col_name) := as.numeric(cut(dash[[i]], breaks=quin_boundaries, include.lowest = TRUE, labels = 1:5, right = TRUE)))
}
# 🏮 上面的代码，最好改成下面这样
#for (item in QUIN) {  dash <- dash %>% mutate(!!paste0("quin", i) := ntile(.data[[i]], 5))  }

dash <- dash %>% mutate(
    quinsugar = ifelse(!is.na(quinsugar), 6 - quinsugar, NA),
    quinmeat = ifelse(!is.na(quinmeat), 6 - quinmeat, NA),
    quinsodium = ifelse(!is.na(quinsodium), 6 - quinsodium, NA)
)

dash <- dash %>% mutate(dashscore = rowSums2(select(., quinvegetable, quinfruit, quinnut, quingrain, quinlowfatdairy, quinsugar, quinmeat, quinsodium)))
dash$dashscore[apply(is.na(select(dash, vegetablenew, fruitnew, nutnew, grainnew, lowfatdairy, sugarnew, meatnew, sodium)), 1, any)] <- NA
# dash <- dash %>% arrange(dashscore)
dash <- dash %>% mutate( dash_pts = case_when((dashscore >= 0 & dashscore < 17) ~ 0, (dashscore >= 17 & dashscore < 21) ~ 25, (dashscore >= 21 & dashscore < 26) ~ 50, (dashscore >= 26 & dashscore < 31) ~ 80, (dashscore >= 31) ~ 100, TRUE ~ NA_real_ ))
saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# LE8
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phe <- readRDS(paste0(indir,"/Rdata/ukb.phe.rds"))
srd <- read.table(paste0(indir,"/rap/srd.tab"), header=T)
srd.drug <- srd %>% mutate( # 💊
	drug.dm=as.integer(rowSums(across(starts_with("p20002_i0_"),  ~ . %in% c("1220", "1221", "1222", "1223")), na.rm = TRUE) > 0),
	drug.hpt=as.integer(rowSums(across(starts_with("p20002_i0_"),  ~ . %in% c("1065","1072")), na.rm = TRUE) > 0)
	) %>% select(eid, drug.dm, drug.hpt)
le8 <- merge(phe, srd.drug, by="eid") %>% mutate( 
	drug.lipid = ifelse(rowSums(across(starts_with("drug.cmd"), ~ . %like% "1")) > 0, 1, 0)
)
le8 <- le8 %>% mutate(
	# 运动 🏃‍ MET minutes per week for moderate [22038] or vigorous [22039] activity
	pa_mins = pa_mod + pa_vig * 2, 
	pa_pts = case_when(pa_mins >=150 ~ 100, (pa_mins>=120 & pa_mins<150) ~ 90, (pa_mins>=90 & pa_mins<120) ~ 80, (pa_mins >=60 & pa_mins<90) ~ 60, (pa_mins>=30 & pa_mins<60) ~ 40, (pa_mins>= 1 & pa_mins< 30) ~ 20, pa_mins==0 ~ 0),
	# 吸烟 🚬 smoke_quit_age [p2897], smoke_history [p1249], smoke_in_house [p1259]
	smoke_quit_age = ifelse(smoke_quit_age %in% c(-3, -1), NA, smoke_quit_age), 
	smoke_history = ifelse(smoke_history==-3, NA, smoke_history), 
	smoke_in_house = ifelse(smoke_in_house==-3, NA, smoke_in_house),
	smoke_quit_year = ifelse(smoke_quit_age > 0 & age>=smoke_quit_age, age-smoke_quit_age, NA), 
	smoke_cat = case_when(smoke_status==0 ~ 100, smoke_status==1 & smoke_quit_year>=5 ~ 75, smoke_history==3 ~ 75, smoke_status==1 & smoke_quit_year>=1 & smoke_quit_year<5 ~ 50, smoke_history==2 ~ 50, smoke_status==1 & smoke_quit_year<1 ~ 25, smoke_status==2 ~ 0),
	smoke_pts = ifelse((smoke_in_house==1 | smoke_in_house==2) & smoke_cat>0, smoke_cat-20, smoke_cat),
	# 睡眠🛏 sleep_duration [p1160]
	sleep_duration = ifelse(sleep_duration %in% c(-1, -3), NA, sleep_duration),
	sleep_pts = case_when(sleep_duration >= 7 & sleep_duration < 9 ~ 100, sleep_duration >= 9 & sleep_duration < 10 ~ 90, sleep_duration >= 6 & sleep_duration < 7 ~ 70, (sleep_duration >= 5 & sleep_duration < 6) | sleep_duration >= 10 ~ 40, sleep_duration >= 4 & sleep_duration < 5 ~ 20, sleep_duration >= 0 & sleep_duration < 4 ~ 0),
	# BMI 🎈
	bmi_pts = case_when(bmi>0 & bmi<25 ~ 100, bmi>=25 & bmi<30 ~ 70, bmi>=30 & bmi<35 ~ 30, bmi>= 35 & bmi<40 ~ 15, bmi>=40 ~ 0),
	# 血脂 🍟  drug.cmd[6177] drug.cmd2[6153]
	nonhdl = pmax((bb_TC-bb_HDL)*38.67, 0),
	nonhdl_cat = case_when(nonhdl>= 0 & nonhdl<130 ~ 100, nonhdl>=130 & nonhdl<160 ~ 60, nonhdl>=160 & nonhdl<190 ~ 40, nonhdl>=190 & nonhdl<220 ~ 20, nonhdl>=220 ~ 0),
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

dash <- readRDS(paste0(indir,"/Rdata/ukb.dash.rds"))[,c("eid", "dash_pts")] # 🍚
le8 <- merge(le8, dash, by="eid")
le8 <- le8 %>% rowwise() %>% # 🎇 汇总
	mutate( LE8score = mean(c(dash_pts, pa_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, bp_pts), na.rm = FALSE)) %>% 
	ungroup() %>% mutate(CVH_cat = case_when(LE8score < 50 & LE8score >= 0 ~ 0, LE8score < 80 & LE8score >= 50 ~ 1, LE8score >= 80 ~ 2 )) %>%
	select(eid, LE8score, CVH_cat, dash_pts, pa_pts, smoke_pts, sleep_pts, bmi_pts, nonhdl_pts, hba1c_pts, bp_pts)
saveRDS(le8, paste0(indir,"/Rdata/ukb.le8.rds"))