# 说明：在mutate()里面生成的新变量
pacman::p_load(tidyverse, lubridate, purrr)
source("D:/scripts/ukb/le8.dash.fields")

indir="D:/data/ukb/phe"

pick_first_number <- function(x) {
	as.numeric(strsplit(x, "\\|")[[1]][1]) # for strings such as "1|2|3"
}
rowSums2 <- function(x, na.rm = TRUE) { # 🏮 如果所有的row都是NA，那么结果也是NA，而不是0
  row_sums <- rowSums(x, na.rm = na.rm)
  row_sums[row_sums == 0 & rowSums(is.na(x)) == ncol(x)] <- NA
  return(row_sums)
}
bulk_rowMeans <- function(dat, fields_names) {
	dat %>% mutate( purrr::imap_dfc(fields_names, function(x, name) {
	prefix <- strsplit(x, "\\|")[[1]][1]
	new_col <- rowMeans(select(dat, starts_with(prefix)), na.rm = TRUE)
	setNames(as.data.frame(new_col), strsplit(x, "\\|")[[1]][2])
	})) %>% mutate_all(~ ifelse(is.nan(.), NA, .)) # 🏮 non-numeric 的数据，会被转换成 NaN，而不是NA
}
replace_if_equal <- function(dat, var_list, num) {
  sapply(var_list, function(pair) {
    fields <- strsplit(pair, "\\|")[[1]]
    ifelse(dat[[fields[1]]] == num, dat[[fields[2]]], 0)
  })
}


dat0 <- read.table(paste0(indir,"/rap/le8.pku.tab.gz"), sep="\t", header=TRUE, flush=TRUE) %>% mutate(
	across(where(is.numeric), ~ case_when(. == 555 ~ 0.5, . == 444 ~ 0.25, . == 200 ~ 2, . == 300 ~ 3, . == 400 ~ 4, . == 500 ~ 5, . == 600 ~ 6, TRUE ~ .))) %>% mutate(
	across(c(p104000_i0:p104590_i0),  ~ ifelse(p20085_i0 %in% c(3, 4, 6), NA, .)), 
	across(c(p104000_i1:p102970_i1),  ~ ifelse(p20085_i1 %in% c(3, 4, 6), NA, .)), 
	across(c(p104000_i2:p102970_i2),  ~ ifelse(p20085_i2 %in% c(3, 4, 6), NA, .)), 
	across(c(p104000_i3:p102970_i3),  ~ ifelse(p20085_i3 %in% c(3, 4, 6), NA, .)), 
	across(c(p104000_i4:p102970_i4),  ~ ifelse(p20085_i4 %in% c(3, 4, 6), NA, .))
)

dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans(select(., starts_with("p100002_")), na.rm = TRUE)
)

dat <- bulk_rowMeans(dat, veg_fields_names) %>% mutate( # 蔬菜 🥦 
	Vegetables = rowSums2(select(., sub(".*\\|", "", veg_fields_names)), na.rm = TRUE),
	Vegetablesnew = ifelse(rowSums2(is.na(select(., starts_with("p103990_")))) > 0 & is.na(Vegetables), 0, Vegetables)
) 

dat <- bulk_rowMeans(dat, fruit_fields_names) %>% mutate( # 水果🍓
	Fruits = rowSums2(select(., sub(".*\\|", "", fruit_fields_names)), na.rm = TRUE),
	Fruitsnew = ifelse(rowSums2(is.na(select(., starts_with("p104400_")))) > 0 & is.na(Fruits), 0, Fruits)
) 

dat <- bulk_rowMeans(dat, nut_fields_names) %>% mutate( # 坚果🌲
	nut = rowSums2(select(., sub(".*\\|", "", nut_fields_names)), na.rm = TRUE),
	nutnew = ifelse(rowSums2(is.na(select(., starts_with("p102400_")))) > 0 & is.na(nut), 0, nut)
)

dat <- bulk_rowMeans(dat, porriage_fields_names) # 主食🍚🍞
dat <- dat %>% mutate( 
    Slicedbread = rowMeans(replace_if_equal(dat, slicedbread_fields, 3), na.rm = TRUE),
    Baguette = rowMeans(replace_if_equal(dat, baguette_fields, 3), na.rm = TRUE),
    Bap = rowMeans(replace_if_equal(dat, bap_fields, 3), na.rm = TRUE),
    Breadroll = rowMeans(replace_if_equal(dat, breadroll_fields, 3), na.rm = TRUE),
    Wholemealpasta = rowMeans(select(., starts_with("p102720_")), na.rm = TRUE),
    Crispbread = rowMeans(select(., starts_with("p101250_")), na.rm = TRUE),
    Oatcakes = rowMeans(select(., starts_with("p101260_")), na.rm = TRUE),
    Otherbread = rowMeans(select(., starts_with("p101270_")), na.rm = TRUE)
) %>% mutate(
	grains = rowSums2(select(., c(sub(".*\\|", "", veg_fields_names), "Slicedbread", "Baguette", "Bap", "Breadroll", "Wholemealpasta", "Crispbread", "Oatcakes", "Otherbread")), na.rm = TRUE),
	grainsnew = ifelse(rowSums2(is.na(select(., starts_with("p100760_")))) > 0 & is.na(grains), 0, grains)
)

dat <- dat %>% mutate( # 🐄🥛
    yogurt_0 = ifelse(p20106_i0 == 210, p102090_i0, ifelse(p20106_i0 == 211, 0, NA)),
    yogurt_1 = ifelse(p20106_i1 == 210, p102090_i1, ifelse(p20106_i1 == 211, 0, NA)),
    yogurt_2 = ifelse(p20106_i2 == 210, p102090_i2, ifelse(p20106_i2 == 211, 0, NA)),
    yogurt_3 = ifelse(p20106_i3 == 210, p102090_i3, ifelse(p20106_i3 == 211, 0, NA)),
    yogurt_4 = ifelse(p20106_i4 == 210, p102090_i4, ifelse(p20106_i4 == 211, 0, NA)),
    Yogurt = rowMeans(select(., starts_with("yogurt_")), na.rm = TRUE),
    milk_0 = ifelse(p100920_i0 %in% c(2102, 2103), p100520_i0, 0),
    milk_1 = ifelse(p100920_i1 %in% c(2102, 2103), p100520_i1, 0),
    milk_2 = ifelse(p100920_i2 %in% c(2102, 2103), p100520_i2, 0),
    milk_3 = ifelse(p100920_i3 %in% c(2102, 2103), p100520_i3, 0),
    milk_4 = ifelse(p100920_i4 %in% c(2102, 2103), p100520_i4, 0),
    milk_0 = ifelse(is.na(p100920_i0), NA, milk_0),
    milk_1 = ifelse(is.na(p100920_i1), NA, milk_1),
    milk_2 = ifelse(is.na(p100920_i2), NA, milk_2),
    milk_3 = ifelse(is.na(p100920_i3), NA, milk_3),
    milk_4 = ifelse(is.na(p100920_i4), NA, milk_4),
    Milk = rowMeans(select(., starts_with("milk_")), na.rm = TRUE),
    hardcheese = rowMeans(select(., starts_with("p102810_")), na.rm = TRUE),
    cheesespread = rowMeans(select(., starts_with("p102850_")), na.rm = TRUE),
    Cottagecheese = rowMeans(select(., starts_with("p102870_")), na.rm = TRUE)
) %>% mutate(
    cheese = rowSums2(select(., hardcheese, cheesespread, Cottagecheese), na.rm = TRUE),
    cheesenew = ifelse((!is.na(p102800_i0) | !is.na(p102800_i1) | !is.na(p102800_i2) |  !is.na(p102800_i3) | !is.na(p102800_i4)) & is.na(cheese), 0, cheese)
) %>% mutate(
    lowfatdairy = rowSums2(select(., Milk, Yogurt, cheesenew), na.rm = TRUE)
)

dat <- bulk_rowMeans(dat, sugar_fields_names) %>% mutate( # 🍬🥤
	sugar = rowSums2(select(., sub(".*\\|", "", sugar_fields_names)), na.rm = TRUE),
	sugarnew = ifelse(rowSums2(is.na(select(., starts_with("p100020_")))) > 0 & is.na(sugar), 0, sugar)
) 

dat <- bulk_rowMeans(dat, meat_fields_names) %>% mutate( # 红肉🥩
	meat = rowSums2(select(., sub(".*\\|", "", meat_fields_names)), na.rm = TRUE),
	meatnew = ifelse(rowSums2(is.na(select(., starts_with("p103000_")))) > 0 & is.na(meat), 0, meat),
) 
 
dat <- dat %>% mutate( # 盐 ⛵
	Sodium = p30530_i0
)

dash <- dat %>% select(eid, Vegetablesnew, Fruitsnew, nutnew, grainsnew, lowfatdairy, sugarnew, meatnew, Sodium)

dash <- dash %>% mutate(
	quinVegetables = ntile(Vegetablesnew, 5),
	quinFruits = ntile(Fruitsnew, 5),
	quinnut = ntile(nutnew, 5),
	quingrains = ntile(grainsnew, 5),
	quinlowfatdairy = ntile(lowfatdairy, 5),
	quinsugar = ntile(sugarnew, 5),
	quinmeat = ntile(meatnew, 5),
	quinSodium = ntile(Sodium, 5),
	quinsugar = ifelse(!is.na(quinsugar), 6 - quinsugar, NA),
	quinmeat = ifelse(!is.na(quinmeat), 6 - quinmeat, NA),
	quinSodium_new = ifelse(!is.na(quinSodium), 6 - quinSodium, NA)
)

dash <- dash %>% mutate(
	dashscore = rowSums2(select(., c("quinVegetables", "quinFruits", "quinnut", "quingrains", "quinlowfatdairy", "quinsugar", "quinmeat", "quinSodium")), na.rm = TRUE),
	dashscore = ifelse(is.na(Vegetablesnew) | is.na(Fruitsnew) | is.na(nutnew) | is.na(grainsnew) | is.na(lowfatdairy) | is.na(sugarnew) | is.na(meatnew) | is.na(Sodium), NA, dashscore),
	dashscore_pctile = cut(dashscore, breaks = c(-Inf, unique(quantile(dashscore, probs = c(0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.93, 0.94, 0.95), na.rm = TRUE)), Inf), labels = FALSE, right = TRUE),
	dashpts = case_when(dashscore > 0 & dashscore < 17 ~ 0, dashscore >= 17 & dashscore < 21 ~ 25, dashscore >= 21 & dashscore < 26 ~ 50, dashscore >= 26 & dashscore < 31 ~ 80, dashscore >= 31 ~ 100, TRUE ~ NA_real_)
)

saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))