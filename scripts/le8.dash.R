pacman::p_load(tidyverse, lubridate, purrr)

dir0="D:"
indir=paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, "/scripts/ukb/le8.dash.fields"))
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- read.table("D:/github/001/scripts/le8/le8.pku.raw.gz", sep="\t", fill=TRUE, header=TRUE) 
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
for (i in QUIN) {
  quin_boundaries <- quantile(dash[[i]], probs = c(0, 0.2, 0.4, 0.6, 0.8, 1), na.rm = TRUE)
  quin_boundaries[2] <- ifelse(quin_boundaries[2] == 0, quin_boundaries[2] + 0.000001, quin_boundaries[2])
  quin_boundaries[3] <- ifelse(quin_boundaries[3] == 0, quin_boundaries[3] + 0.000002, quin_boundaries[3])
  quin_boundaries[4] <- ifelse(quin_boundaries[4] == 0, quin_boundaries[4] + 0.000003, quin_boundaries[4])
  quin_boundaries[5] <- ifelse(quin_boundaries[5] == 0, quin_boundaries[5] + 0.000004, quin_boundaries[5])
  new_col_name <- paste0("quin", gsub("new", "", i)) 
  dash <- dash %>% mutate(!!sym(new_col_name) := as.numeric(cut(dash[[i]],breaks=quin_boundaries, include.lowest=TRUE,labels=c(1, 2, 3, 4, 5),right=TRUE)))
}
# 🏮 上面的代码，最好改成下面这样
#for (item in QUIN) { dash <- dash %>% mutate(!!paste0("quin", item) := ntile(.data[[item]], 5)) }

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

