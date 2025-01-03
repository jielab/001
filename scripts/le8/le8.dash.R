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
	across(where(is.numeric), ~ case_when(. == 555 ~ 0.5, . == 444 ~ 0.25, . == 200 ~ 2, . == 300 ~ 3, . == 400 ~ 4, . == 500 ~ 5, . == 600 ~ 6, TRUE ~ .))
) 
	
dat0 <- dat0 %>% mutate(
	across(.cols = matches("i0$") & !matches("p20085_i0"),.fns = ~ ifelse(p20085_i0 %in% c(3, 4, 6), NA, .)),
	across(.cols = matches("i1$") & !matches("p20085_i1"),.fns = ~ ifelse(p20085_i1 %in% c(3, 4, 6), NA, .)),
	across(.cols = matches("i2$") & !matches("p20085_i2"),.fns = ~ ifelse(p20085_i2 %in% c(3, 4, 6), NA, .)),
	across(.cols = matches("i3$") & !matches("p20085_i3"),.fns = ~ ifelse(p20085_i3 %in% c(3, 4, 6), NA, .)),
	across(.cols = matches("i4$") & !matches("p20085_i4"),.fns = ~ ifelse(p20085_i4 %in% c(3, 4, 6), NA, .))
)

dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans(select(., starts_with("p100002_")), na.rm = TRUE)
)

dat <- bulk_rowMeans(dat, veg_fields_names) # 蔬菜 🥦 
dat <- dat %>% mutate( 
	vegetable = rowSums2(select(., sub(".*\\|", "", veg_fields_names)), na.rm = TRUE),
	vegetablenew= ifelse((rowSums(across(starts_with("p103990_i"), ~ !is.na(.))) > 0 & is.na(vegetable)), 0, vegetable)
)

dat <- bulk_rowMeans(dat, fruit_fields_names) # 水果🍓
dat <- dat %>% mutate( 
	fruit = rowSums2(select(., sub(".*\\|", "", fruit_fields_names)), na.rm = TRUE),
	fruitnew= ifelse((rowSums(across(starts_with("p104400_i"), ~ !is.na(.))) > 0 & is.na(fruit)), 0, fruit)
) 

dat <- bulk_rowMeans(dat, nut_fields_names) # 坚果🌲
dat <- dat %>% mutate( 
	nut = rowSums2(select(., sub(".*\\|", "", nut_fields_names)), na.rm = TRUE),
	nutnew= ifelse((rowSums(across(starts_with("p102400_i"), ~ !is.na(.))) > 0 & is.na(nut)), 0, nut)
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
) 
dat <- dat %>% mutate(
	grain = rowSums2(select(., c(sub(".*\\|", "", veg_fields_names), "Slicedbread", "Baguette", "Bap", "Breadroll", "Wholemealpasta", "Crispbread", "Oatcakes", "Otherbread")), na.rm = TRUE),
	grainnew= ifelse((rowSums(across(starts_with("p100760_i"), ~ !is.na(.))) > 0 & is.na(grain)), 0, grain)
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
) 
dat <- dat %>% mutate( 
    cheese = rowSums2(select(., hardcheese, cheesespread, Cottagecheese), na.rm = TRUE),
	cheesenew= ifelse((rowSums(across(starts_with("p102800_i"), ~ !is.na(.))) > 0 & is.na(cheese)), 0, cheese)
) 
dat <- dat %>% mutate(
    lowfatdairy = rowSums2(select(., Milk, Yogurt, cheesenew), na.rm = TRUE)
)

dat <- bulk_rowMeans(dat, sugar_fields_names) # 🍬🥤
dat <- dat %>% mutate( 
	sugar = rowSums2(select(., sub(".*\\|", "", sugar_fields_names)), na.rm = TRUE),
	sugarnew= ifelse((rowSums(across(starts_with("p100020_i"), ~ !is.na(.))) > 0 & is.na(sugar)), 0, sugar)
) 

dat <- bulk_rowMeans(dat, meat_fields_names) # 红肉🥩
dat <- dat %>% mutate(
	meat = rowSums2(select(., sub(".*\\|", "", meat_fields_names)), na.rm = TRUE),
	meatnew= ifelse((rowSums(across(starts_with("p103000_i"), ~ !is.na(.))) > 0 & is.na(meat)), 0, meat)
)
 
dat <- dat %>% mutate( # 盐 ⛵
	sodium = p30530_i0
)

dash <- dat %>% select(eid, vegetablenew, fruitnew, nutnew, grainnew, lowfatdairy, sugarnew, meatnew, sodium)

dash <- dash %>% mutate(
	quinvegetable = ntile(vegetablenew, 5),
	quinfruit = ntile(fruitnew, 5),
	quinnut = ntile(nutnew, 5),
	quingrain = ntile(grainnew, 5),
	quinlowfatdairy = ntile(lowfatdairy, 5),
	quinsugar = ntile(sugarnew, 5),
	quinmeat = ntile(meatnew, 5),
	quinsodium = ntile(sodium, 5),
	quinsugar = ifelse(!is.na(quinsugar), 6 - quinsugar, NA),
	quinmeat = ifelse(!is.na(quinmeat), 6 - quinmeat, NA),
	quinsodium_new = ifelse(!is.na(quinsodium), 6 - quinsodium, NA)
)

dash <- dash %>% mutate(
	dashscore = rowSums2(select(., c("quinvegetable", "quinfruit", "quinnut", "quingrain", "quinlowfatdairy", "quinsugar", "quinmeat", "quinsodium")), na.rm = TRUE),
	dashscore = ifelse(is.na(vegetablenew) | is.na(fruitnew) | is.na(nutnew) | is.na(grainnew) | is.na(lowfatdairy) | is.na(sugarnew) | is.na(meatnew) | is.na(sodium), NA, dashscore),
	dashscore_pctile = cut(dashscore, breaks = c(-Inf, unique(quantile(dashscore, probs = c(0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.93, 0.94, 0.95), na.rm = TRUE)), Inf), labels = FALSE, right = TRUE),
	dashpts = case_when(dashscore > 0 & dashscore < 17 ~ 0, dashscore >= 17 & dashscore < 21 ~ 25, dashscore >= 21 & dashscore < 26 ~ 50, dashscore >= 26 & dashscore < 31 ~ 80, dashscore >= 31 ~ 100, TRUE ~ NA_real_)
)

saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))