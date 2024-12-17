# 说明：凡是在 mutate()里面生成的新变量，在新的mutate()里面才能使用

setwd("D:/")
pacman::p_load(tidyverse, lubridate)
indir="D:/data/ukb/phe"

dat0 <- read.table(paste0(indir,"/rap/le8.tab"), sep="\t", header=TRUE, flush=TRUE) %>% 
	mutate(across(where(is.numeric), ~ case_when(. == 555 ~ 0.5, . == 444 ~ 0.25, . == 200 ~ 2, . == 300 ~ 3, . == 400 ~ 4, . == 500 ~ 5, . == 600 ~ 6, TRUE ~ .))) %>%
	mutate(across(c(p104000_i0:p104590_i0),  ~ ifelse(p20085_i0 %in% c(3, 4, 6), NA, .)), across(c(p104000_i1:p102970_i1),  ~ ifelse(p20085_i1 %in% c(3, 4, 6), NA, .)), across(c(p104000_i2:p102970_i2),  ~ ifelse(p20085_i2 %in% c(3, 4, 6), NA, .)), across(c(p104000_i3:p102970_i3),  ~ ifelse(p20085_i3 %in% c(3, 4, 6), NA, .)), across(c(p104000_i4:p102970_i4),  ~ ifelse(p20085_i4 %in% c(3, 4, 6), NA, .)))

dat <- dat0 %>% mutate( # 热量 🌋
	energy_total = rowMeans(select(., starts_with("p100002_")), na.rm = TRUE)
)

dat <- dat %>% mutate( # 蔬菜 🥦
    mixveg = rowMeans(select(., starts_with("p104060_")), na.rm = TRUE),
    vegpiece = rowMeans(select(., starts_with("p104070_")), na.rm = TRUE),
    Coleslaw = rowMeans(select(., starts_with("p104080_")), na.rm = TRUE),
    Side = rowMeans(select(., starts_with("p104090_")), na.rm = TRUE),
    Avocado = rowMeans(select(., starts_with("p104100_")), na.rm = TRUE),
    Beetroot = rowMeans(select(., starts_with("p104130_")), na.rm = TRUE),
    Broccoli = rowMeans(select(., starts_with("p104140_")), na.rm = TRUE),
    Butternut = rowMeans(select(., starts_with("p104150_")), na.rm = TRUE),
    Cabbage = rowMeans(select(., starts_with("p104160_")), na.rm = TRUE),
    Carrot = rowMeans(select(., starts_with("p104170_")), na.rm = TRUE),
    Cauliflower = rowMeans(select(., starts_with("p104180_")), na.rm = TRUE),
    Celery = rowMeans(select(., starts_with("p104190_")), na.rm = TRUE),
    Courgette = rowMeans(select(., starts_with("p104200_")), na.rm = TRUE),
    Cucumber = rowMeans(select(., starts_with("p104210_")), na.rm = TRUE),
    Garlic = rowMeans(select(., starts_with("p104220_")), na.rm = TRUE),
    Leek = rowMeans(select(., starts_with("p104230_")), na.rm = TRUE),
    Lettuce = rowMeans(select(., starts_with("p104240_")), na.rm = TRUE),
    Mushroom = rowMeans(select(., starts_with("p104250_")), na.rm = TRUE),
    Onion = rowMeans(select(., starts_with("p104260_")), na.rm = TRUE),
    Parsnip = rowMeans(select(., starts_with("p104270_")), na.rm = TRUE),
    pepper = rowMeans(select(., starts_with("p104290_")), na.rm = TRUE),
    Spinach = rowMeans(select(., starts_with("p104300_")), na.rm = TRUE),
    Sprouts = rowMeans(select(., starts_with("p104310_")), na.rm = TRUE),
    Sweetcorn = rowMeans(select(., starts_with("p104320_")), na.rm = TRUE),
    Freshtomato = rowMeans(select(., starts_with("p104340_")), na.rm = TRUE),
    Tinnedtomato = rowMeans(select(., starts_with("p104350_")), na.rm = TRUE),
    Turnip = rowMeans(select(., starts_with("p104360_")), na.rm = TRUE),
    Watercress = rowMeans(select(., starts_with("p104370_")), na.rm = TRUE),
    Otherveg = rowMeans(select(., starts_with("p104380_")), na.rm = TRUE)
  )
dat <- dat %>% mutate( 
	Vegetables = rowSums(select(., c("mixveg", "vegpiece", "Coleslaw", "Side", "Avocado", "Beetroot", "Broccoli", "Butternut", "Cabbage", "Carrot", "Cauliflower", "Celery", "Courgette", "Cucumber", "Garlic", "Leek", "Lettuce", "Mushroom", "Onion", "Parsnip", "pepper", "Spinach", "Sprouts", "Sweetcorn", "Freshtomato", "Tinnedtomato", "Turnip", "Watercress", "Otherveg")), na.rm = TRUE),
	Vegetablesnew = ifelse(rowSums(is.na(select(., starts_with("p103990_")))) > 0 & is.na(Vegetables), 0, Vegetables)
) 
 
dat <- dat %>% mutate( # 水果🍓
    Stewedfruit = rowMeans(select(., starts_with("p104410_")), na.rm = TRUE),
    Prune = rowMeans(select(., starts_with("p104420_")), na.rm = TRUE),
    Dried = rowMeans(select(., starts_with("p104430_")), na.rm = TRUE),
    Mixedfruit = rowMeans(select(., starts_with("p104440_")), na.rm = TRUE),
    Apple = rowMeans(select(., starts_with("p104450_")), na.rm = TRUE),
    Banana = rowMeans(select(., starts_with("p104460_")), na.rm = TRUE),
    Berry = rowMeans(select(., starts_with("p104470_")), na.rm = TRUE),
    Cherry = rowMeans(select(., starts_with("p104480_")), na.rm = TRUE),
	Grapefruit = rowMeans(select(., starts_with("p104490_")), na.rm = TRUE),
    Grape = rowMeans(select(., starts_with("p104500_")), na.rm = TRUE),
    Mango = rowMeans(select(., starts_with("p104510_")), na.rm = TRUE),
    Melon = rowMeans(select(., starts_with("p104520_")), na.rm = TRUE),
    Orange = rowMeans(select(., starts_with("p104530_")), na.rm = TRUE),
    Satsuma = rowMeans(select(., starts_with("p104540_")), na.rm = TRUE),
    Peach = rowMeans(select(., starts_with("p104550_")), na.rm = TRUE),
    Pear = rowMeans(select(., starts_with("p104560_")), na.rm = TRUE),
    Pineapple = rowMeans(select(., starts_with("p104570_")), na.rm = TRUE),
    Plum = rowMeans(select(., starts_with("p104580_")), na.rm = TRUE),
    Otherfruit = rowMeans(select(., starts_with("p104590_")), na.rm = TRUE),
    Orangejuice = rowMeans(select(., starts_with("p100190_")), na.rm = TRUE),
	Grapefruitjuicejuice = rowMeans(select(., p100200_i0, p100200_i1, p100200_i2, p100200_i3, p100200_i4), na.rm = TRUE),
    Purefruitjuice = rowMeans(select(., starts_with("p100210_")), na.rm = TRUE)
)
dat <- dat %>% mutate(Fruits = rowSums(select(., c("Stewedfruit", "Prune", "Dried", "Mixedfruit", "Apple",  "Banana", "Berry", "Cherry", "Grapefruit", "Grape", "Mango",  "Melon", "Orange", "Satsuma", "Peach", "Pear", "Pineapple", "Plum", "Otherfruit", "Orangejuice", "Grapefruitjuice", "Purefruitjuice")), na.rm = TRUE),
	Fruitsnew = ifelse(rowSums(is.na(select(., starts_with("p104400_")))) > 0 & is.na(Fruits), 0, Fruits)
) 
 
dat <- dat %>% mutate( # 坚果🌲
    Saltpeanut = rowMeans(select(., starts_with("p102410_")), na.rm = TRUE),
    Unsaltpeanut = rowMeans(select(., starts_with("p102420_")), na.rm = TRUE),
    Saltnut = rowMeans(select(., starts_with("p102430_")), na.rm = TRUE),
    Unsaltnut = rowMeans(select(., starts_with("p102440_")), na.rm = TRUE),
    Seeds = rowMeans(select(., starts_with("p102450_")), na.rm = TRUE),
    Tofu = rowMeans(select(., starts_with("p103270_")), na.rm = TRUE),
    Bakedbean = rowMeans(select(., starts_with("p104000_")), na.rm = TRUE),
    Pulses = rowMeans(select(., starts_with("p104010_")), na.rm = TRUE),
    Broadbean = rowMeans(select(., starts_with("p104110_")), na.rm = TRUE),
    Greenbean = rowMeans(select(., starts_with("p104120_")), na.rm = TRUE),
    Pea = rowMeans(select(., starts_with("p104280_")), na.rm = TRUE)
)
dat <- dat %>% mutate(
	nuts = rowSums(select(., c("Saltpeanut", "Unsaltpeanut", "Saltnut", "Unsaltnut", "Seeds", "Tofu", "Bakedbean", "Pulses", "Broadbean", "Greenbean", "Pea")), na.rm = TRUE),
    nutsnew = ifelse(rowSums(is.na(select(., starts_with("p102400_")))) > 0 & is.na(nuts), 0, nuts)
)
 
dat <- dat %>% mutate( # 粗粮🍞
    Porridge = rowMeans(select(., starts_with("p100770_")), na.rm = TRUE),
    Muesli = rowMeans(select(., starts_with("p100800_")), na.rm = TRUE),
    Oat = rowMeans(select(., starts_with("p100810_")), na.rm = TRUE),
    Sweetcereal = rowMeans(select(., starts_with("p100820_")), na.rm = TRUE),
    Plaincereal = rowMeans(select(., starts_with("p100830_")), na.rm = TRUE),
    Brancereal = rowMeans(select(., starts_with("p100840_")), na.rm = TRUE),
    Wholewheat = rowMeans(select(., starts_with("p100850_")), na.rm = TRUE),
    Othercereal = rowMeans(select(., starts_with("p100860_")), na.rm = TRUE),
    Brownrice = rowMeans(select(., starts_with("p102740_")), na.rm = TRUE),
    Othergrain = rowMeans(select(., starts_with("p102780_")), na.rm = TRUE),
    Couscous = rowMeans(select(., starts_with("p102770_")), na.rm = TRUE),
    Slicedbread = ifelse(p20091_i0 == 3, p100950_i0, 0),
    Slicedbread_0_0 = ifelse(is.na(p20091_i0), NA, Slicedbread),
    Slicedbread_1_0 = ifelse(p20091_i1 == 3, p100950_i1, 0),
    Slicedbread_1_0 = ifelse(is.na(p20091_i1), NA, Slicedbread_1_0),
    Slicedbread_2_0 = ifelse(p20091_i2 == 3, p100950_i2, 0),
    Slicedbread_2_0 = ifelse(is.na(p20091_i2), NA, Slicedbread_2_0),
    Slicedbread_3_0 = ifelse(p20091_i3 == 3, p100950_i3, 0),
    Slicedbread_3_0 = ifelse(is.na(p20091_i3), NA, Slicedbread_3_0),
    Slicedbread_4_0 = ifelse(p20091_i4 == 3, p100950_i4, 0),
    Slicedbread_4_0 = ifelse(is.na(p20091_i4), NA, Slicedbread_4_0),
    Slicedbread = rowMeans(select(., starts_with("Slicedbread_")), na.rm = TRUE),
    Baguette = ifelse(p20092_i0 == 3, p101020_i0, 0),
    Baguette_0_0 = ifelse(is.na(p20092_i0), NA, Baguette),
    Baguette_1_0 = ifelse(p20092_i1 == 3, p101020_i1, 0),
    Baguette_1_0 = ifelse(is.na(p20092_i1), NA, Baguette_1_0),
    Baguette_2_0 = ifelse(p20092_i2 == 3, p101020_i2, 0),
    Baguette_2_0 = ifelse(is.na(p20092_i2), NA, Baguette_2_0),
    Baguette_3_0 = ifelse(p20092_i3 == 3, p101020_i3, 0),
    Baguette_3_0 = ifelse(is.na(p20092_i3), NA, Baguette_3_0),
    Baguette_4_0 = ifelse(p20092_i4 == 3, p101020_i4, 0),
    Baguette_4_0 = ifelse(is.na(p20092_i4), NA, Baguette_4_0),
    Baguette = rowMeans(select(., starts_with("Baguette_")), na.rm = TRUE),
    Bap = ifelse(p20093_i0 == 3, p101090_i0, 0),
    Bap_0_0 = ifelse(is.na(p20093_i0), NA, Bap),
    Bap_1_0 = ifelse(p20093_i1 == 3, p101090_i1, 0),
    Bap_1_0 = ifelse(is.na(p20093_i1), NA, Bap_1_0),
    Bap_2_0 = ifelse(p20093_i2 == 3, p101090_i2, 0),
    Bap_2_0 = ifelse(is.na(p20093_i2), NA, Bap_2_0),
    Bap_3_0 = ifelse(p20093_i3 == 3, p101090_i3, 0),
    Bap_3_0 = ifelse(is.na(p20093_i3), NA, Bap_3_0),
    Bap_4_0 = ifelse(p20093_i4 == 3, p101090_i4, 0),
    Bap_4_0 = ifelse(is.na(p20093_i4), NA, Bap_4_0),
    Bap = rowMeans(select(., starts_with("Bap_")), na.rm = TRUE),
    Breadroll = ifelse(p20094_i0 == 3, p101160_i0, 0),
    Breadroll_0_0 = ifelse(is.na(p20094_i0), NA, Breadroll),
    Breadroll_1_0 = ifelse(p20094_i1 == 3, p101160_i1, 0),
    Breadroll_1_0 = ifelse(is.na(p20094_i1), NA, Breadroll_1_0),
    Breadroll_2_0 = ifelse(p20094_i2 == 3, p101160_i2, 0),
    Breadroll_2_0 = ifelse(is.na(p20094_i2), NA, Breadroll_2_0),
    Breadroll_3_0 = ifelse(p20094_i3 == 3, p101160_i3, 0),
    Breadroll_3_0 = ifelse(is.na(p20094_i3), NA, Breadroll_3_0),
    Breadroll_4_0 = ifelse(p20094_i4 == 3, p101160_i4, 0),
    Breadroll_4_0 = ifelse(is.na(p20094_i4), NA, Breadroll_4_0),
    Breadroll = rowMeans(select(., starts_with("Breadroll_")), na.rm = TRUE),
    Wholemealpasta = rowMeans(select(., starts_with("p102720_")), na.rm = TRUE),
    Crispbread = rowMeans(select(., starts_with("p101250_")), na.rm = TRUE),
    Oatcakes = rowMeans(select(., starts_with("p101260_")), na.rm = TRUE),
    Otherbread = rowMeans(select(., starts_with("p101270_")), na.rm = TRUE)
)
dat <- dat %>% mutate(
	grains = rowSums(select(., c("Porridge", "Muesli", "Oat", "Sweetcereal", "Plaincereal", "Brancereal", "Wholewheat", "Othercereal", "Brownrice", "Othergrain", "Couscous", "Slicedbread", "Baguette", "Bap", "Breadroll", "Wholemealpasta", "Crispbread", "Oatcakes", "Otherbread")), na.rm = TRUE),
	grainsnew = ifelse(rowSums(is.na(select(., starts_with("p100760_")))) > 0 & is.na(grains), 0, grains)
)
 
dat <- dat %>% mutate( # 🐄🥛
    yogurt = ifelse(p20106_i0 == 210, p102090_i0, 0),
    yogurt_0_0 = ifelse(is.na(p20106_i0), NA, yogurt),
    yogurt_1_0 = ifelse(p20106_i1 == 210, p102090_i1, 0),
    yogurt_1_0 = ifelse(is.na(p20106_i1), NA, yogurt_1_0),
    yogurt_2_0 = ifelse(p20106_i2 == 210, p102090_i2, 0),
    yogurt_2_0 = ifelse(is.na(p20106_i2), NA, yogurt_2_0),
    yogurt_3_0 = ifelse(p20106_i3 == 210, p102090_i3, 0),
    yogurt_3_0 = ifelse(is.na(p20106_i3), NA, yogurt_3_0),
    yogurt_4_0 = ifelse(p20106_i4 == 210, p102090_i4, 0),
    yogurt_4_0 = ifelse(is.na(p20106_i4), NA, yogurt_4_0),
    Yogurt = rowMeans(select(., starts_with("yogurt_")), na.rm = TRUE),
    Milk = ifelse(p100920_i0 %in% c(2102, 2103), p100520_i0, 0),
    milk_0_0 = ifelse(is.na(p100920_i0), NA, Milk),
    milk_1_0 = ifelse(p100920_i1 %in% c(2102, 2103), p100520_i1, 0),
    milk_1_0 = ifelse(is.na(p100920_i1), NA, milk_1_0),
    milk_2_0 = ifelse(p100920_i2 %in% c(2102, 2103), p100520_i2, 0),
    milk_2_0 = ifelse(is.na(p100920_i2), NA, milk_2_0),
    milk_3_0 = ifelse(p100920_i3 %in% c(2102, 2103), p100520_i3, 0),
    milk_3_0 = ifelse(is.na(p100920_i3), NA, milk_3_0),
    milk_4_0 = ifelse(p100920_i4 %in% c(2102, 2103), p100520_i4, 0),
    milk_4_0 = ifelse(is.na(p100920_i4), NA, milk_4_0),
    Milk = rowMeans(select(., starts_with("milk_")), na.rm = TRUE),
    hardcheese = rowMeans(select(., starts_with("p102810_")), na.rm = TRUE),
    cheesespread = rowMeans(select(., starts_with("p102850_")), na.rm = TRUE),
    Cottagecheese = rowMeans(select(., starts_with("p102870_")), na.rm = TRUE),
    cheese = hardcheese + cheesespread + Cottagecheese,
    cheesenew = ifelse(rowSums(is.na(select(., starts_with("p102800_")))) > 0 & is.na(cheese), 0, cheese),
    lowfatdairy = Milk + Yogurt + cheesenew
) 

dat <- dat %>% mutate( # 🍬🥤
    Fizzydrink = rowMeans(select(., starts_with("p100170_")), na.rm = TRUE),
    Squash = rowMeans(select(., starts_with("p100180_")), na.rm = TRUE),
    Fruitsmootie = rowMeans(select(., starts_with("p100220_")), na.rm = TRUE),
    sugarcoffee = rowMeans(select(., p100370_i0, p100370_i1, p100370_i2, p100370_i3, p100370_i4), na.rm = TRUE),
    artisugarcoffee = rowMeans(select(., p100380_i0, p100380_i1, p100380_i2, p100380_i3, p100380_i4), na.rm = TRUE),
    sugartea = rowMeans(select(., p100490_i0, p100490_i1, p100490_i2, p100490_i3, p100490_i4), na.rm = TRUE),
    artisugartea = rowMeans(select(., p100500_i0, p100500_i1, p100500_i2, p100500_i3, p100500_i4), na.rm = TRUE),
    Lowhotchocolate = rowMeans(select(., p100540_i0, p100540_i1, p100540_i2, p100540_i3, p100540_i4), na.rm = TRUE),
    Hotchocolate = rowMeans(select(., p100550_i0, p100550_i1, p100550_i2, p100550_i3, p100550_i4), na.rm = TRUE),
    sugarsweetened = Fizzydrink + Squash + Fruitsmootie,
    sugarsweetenednew = ifelse(rowSums(is.na(select(., starts_with("p100020_")))) > 0 & is.na(sugarsweetened), 0, sugarsweetened), 
) 
 
dat <- dat %>% mutate( # 红肉🥩
    Sausage = rowMeans(select(., starts_with("p103010_")), na.rm = TRUE),
    Beef = rowMeans(select(., starts_with("p103020_")), na.rm = TRUE),
    Pork = rowMeans(select(., starts_with("p103030_")), na.rm = TRUE),
    Lamb = rowMeans(select(., starts_with("p103040_")), na.rm = TRUE),
    Bacon = rowMeans(select(., starts_with("p103070_")), na.rm = TRUE),
    Ham = rowMeans(select(., starts_with("p103080_")), na.rm = TRUE),
    Liver = rowMeans(select(., starts_with("p103090_")), na.rm = TRUE),
    Scotchegg = rowMeans(select(., starts_with("p102970_")), na.rm = TRUE),
    redmeat = Sausage + Beef + Pork + Lamb + Bacon + Ham + Liver + Scotchegg,
    redmeatnew = ifelse(rowSums(is.na(select(., starts_with("p103000_")))) > 0 & is.na(redmeat), 0, redmeat),
) 
 
dat <- dat %>% mutate( # 盐 ⛵
    Sodium = p30530_i0
)

dash <- dat %>% select(eid, energy_total, Vegetablesnew, Fruitsnew, nutsnew, grainsnew, lowfatdairy, sugarsweetenednew, redmeatnew, Sodium)
dash <- dash %>%
mutate(
    quinVegetables = ntile(Vegetablesnew, 5),
    quinFruits = ntile(Fruitsnew, 5),
    quinnuts = ntile(nutsnew, 5),
    quingrains = ntile(grainsnew, 5),
    quinlowfatdairy = ntile(lowfatdairy, 5),
    quinsugarsweetened = ntile(sugarsweetenednew, 5),
    quinredmeat = ntile(redmeatnew, 5),
    quinSodium = ntile(Sodium, 5),
    quinsugar = ifelse(!is.na(quinsugarsweetened), 6 - quinsugarsweetened, NA),
    quinmeat = ifelse(!is.na(quinredmeat), 6 - quinredmeat, NA),
    quinSodium_new = ifelse(!is.na(quinSodium), 6 - quinSodium, NA)
)
dash <- dash %>% mutate(
	dashscore = rowSums(select(., c("quinVegetables", "quinFruits", "quinnuts", "quingrains", "quinlowfatdairy", "quinsugarsweetened", "quinredmeat", "quinSodium")), na.rm = TRUE),
	dashscore = ifelse(is.na(Vegetablesnew) | is.na(Fruitsnew) | is.na(nutsnew) | is.na(grainsnew) | is.na(lowfatdairy) | is.na(sugarsweetenednew) | is.na(redmeatnew) | is.na(Sodium), NA, dashscore)
	dashscore_pctile = cut(dashscore, breaks = c(-Inf, unique(quantile(dashscore, probs = c(0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.93, 0.94, 0.95), na.rm = TRUE)), Inf), labels = FALSE, right = TRUE),
	dashpts = case_when(dashscore > 0 & dashscore < 17 ~ 0, dashscore >= 17 & dashscore < 21 ~ 25, dashscore >= 21 & dashscore < 26 ~ 50, dashscore >= 26 & dashscore < 31 ~ 80, dashscore >= 31 ~ 100, TRUE ~ NA_real_)
)

saveRDS(dash, paste0(indir,"/Rdata/ukb.dash.rds"))