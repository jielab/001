# 复旦（PMID 40603580), UKB category 100090
pacman::p_load(readxl, data.table, tidyverse, lubridate, purrr, survival)

dir0 = "D:"
indir = paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, '/scripts/f/phe.f.R'))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Diet 24-hour recall QC [Prune_data.R]
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS(file = paste0(dir0, '/data/ukb/phe/Rdata/all.rds')) 
diet0 <- readRDS(paste0(indir, "/Rdata/le8.raw.rds")) %>% select(!matches("^medi\\."))
names.new <- sub("^diet(w)?_([0-9]+)", "p\\2", names(diet0)); check_dup(names.new); names(diet0) <- names.new

dura <- diet0 %>% select(eid, starts_with("p20082")) # #1 排除问卷做的太快的，选择至少有1次dietary recalls的
	q <- quantile(as.matrix(dura[, -1]), na.rm = TRUE, probs = c(0.05, 0.95))
	dura2 <- dura %>% mutate(across(-eid, ~ ifelse(.x >= q[1] & .x <= 60 * 24, .x, NA_real_)), enroll = rowSums(!is.na(across(-eid)))); table(dura2$enroll)
	vari <- diet0 %>% select(eid, matches("vari")) %>% mutate(vari = as.integer(rowSums(across(-eid, ~ .x == 3, .names = NULL), na.rm = TRUE) > 0))
	keep <- dura2 %>% select(eid, enroll) %>% inner_join(vari %>% select(eid, vari), by = "eid") %>% filter(enroll >= 1, vari == 0) %>% select(eid)
diet <- diet0 %>% filter(eid %in% keep$eid); nrow(diet) # ✅ 184501

energy <- diet0 %>% select(eid, starts_with("energy_i")) %>% # 排除可疑能量摄入的人群（男女区分) 
	left_join(dat0 %>% select(eid, sex), by = "eid") %>% semi_join(keep, by = "eid") %>% 
	mutate(energy.kJ = rowMeans(across(starts_with("energy_i")), na.rm = TRUE)) %>% 
	filter((sex == 0 & between(energy.kJ, 2092, 14644)) | (sex == 1 & between(energy.kJ, 3347.2, 17572.8))) %>% select(eid, energy.kJ)
diet <- diet %>% filter(eid %in% energy$eid); nrow(diet) # ✅ 182451

inst <- 0:4 # 🍞面包分为黑白，🥥油分为好坏
breads <- c("bread", "baguette", "bap", "roll")
for (k in inst) {
	for (f in breads) { 
		type_col <- sprintf("type_%s_i%d", f, k)  # type_bread_i0
		int_col  <- sprintf("%s_i%d", f, k)       # bread_i0
		total <- str_count(diet[[type_col]], "[1-9]")
		w_n   <- str_count(diet[[type_col]], "[15]")
		b_n   <- str_count(diet[[type_col]], "[2-4]")
		diet[[sprintf("%s.w_i%d", f, k)]] <- ifelse(!is.na(diet[[int_col]]) & total > 0 & w_n > 0, diet[[int_col]] * w_n / total, 0)
		diet[[sprintf("%s.b_i%d", f, k)]] <- ifelse(!is.na(diet[[int_col]]) & total > 0 & b_n > 0, diet[[int_col]] * b_n / total, 0)
	}
	raw_col <- sprintf("type_oil_i%d", k)
	new_col <- sprintf("oil_i%d", k)
	if (raw_col %in% names(diet)) { diet[[new_col]] <- ifelse(str_detect(diet[[raw_col]], "^35[4-7]$"), 10, 0)
	} else  {diet[[new_col]] <- NA_real_ }
}
diet <- diet %>% select(-matches(paste0("^", breads, "_i", collapse = "|")), -matches("^type"))

field <- read_excel(paste0(dir0, "/files/foods.xlsx"))
field$name <- case_when( # 把4个🍞变量重命名
	field$name %in% breads ~ paste0(field$name, ".w"),
	grepl("^type_", field$name) & sub("^type_", "", field$name) %in% breads ~ paste0(sub("^type_", "", field$name), ".b"),
	grepl("^type_", field$name) ~ sub("^type_", "", field$name), TRUE ~ field$name
)

dat1 <- diet %>% select(eid)
	for (i in field$name) {dat1 <- cbind(dat1, select(diet, starts_with(as.character(i))))}
dat2 <- dat1 %>% select(eid) # 🏮取均值
	for (i in unique(sub("_i.*", "", names(dat1)))) { dat2[[i]] <- rowMeans(dat1 %>% select(starts_with(as.character(i))), na.rm = TRUE)}
lapply(dat2, summary); length(names(dat2))
saveRDS(dat2, "Food_items.rds") 

dat.food <- readRDS("Food_items.rds") # 合并食物条目为食物组
wt_vars <- grep("26064|26152|26153|26151|26067|26138|26102|26131|26133|26096|26150|26154|26106", names(dat.food), value = TRUE) # fat不处理 : 26062, 26063, 26111, 26112
	portion <- field %>% filter(as.character(FieldID) %in% wt_vars) %>% select(FieldID, Portion) %>% tibble::deframe()
	dat.food <- dat.food %>% mutate(across(all_of(names(portion)), ~ .x / portion[cur_column()]))  
dat.grp <- dat.food %>% select(eid)
	grp.link <- split(field$name, field$Food.Group)
	for (i in names(grp.link)) {dat.grp[[i]] <- rowSums(as.matrix(dat.food[, grp.link[[i]], drop = FALSE]), na.rm = TRUE)}
lapply(dat.grp, summary)
saveRDS(dat.grp, "Food_groups.rds")


## ------------------------------------------------------------------
## 🚩 疾病变量、 协变量
## ------------------------------------------------------------------
dat.food <- readRDS("Food_items.rds")
dat.grp <- readRDS("Food_groups.rds")
diet_date <- diet %>% select(eid, starts_with("p105010_i")) %>% mutate(across(starts_with("p105010_i"), ~ as.Date(.x)))
	mask_valid <- diet %>% select(starts_with("p20082_i")) %>% as.matrix() %>% { !is.na(.) }
	date_mat <- as.matrix(diet_date[, grep("^p105010_i", names(diet_date), value = TRUE)])
	date_mat[!mask_valid] <- NA
	diet_date[, grep("^p105010_i", names(diet_date), value = TRUE)] <- date_mat
	diet_date$last_date <- apply(date_mat, 1, function(x) {if (all(is.na(x))) NA_Date_ else max(x, na.rm = TRUE)}) %>% as.Date()
	diet_date <- diet_date %>% select(eid, last_date)
dat.list <- list(dat.grp, dat.food, diet_date, energy, dat0) %>% map(~ mutate(.x, eid = as.character(eid)))
dat <- reduce(dat.list, left_join, by = "eid")

covs.basic <- 'energy.kJ age sex.c ethnic.b tdi edu.sco smoke.c pa_met_tot bmi apoe.e4' %>% strsplit(' ') %>% unlist()
cov.Y <- c("htn", "dm", "cvd", "depress")
covs <- c(covs.basic, paste0("hist_", cov.Y))
dat <- dat %>% mutate(
	ethnic.b = factor(ifelse(ethnic.c=="White", "White", "Others"), levels=c("White", "Others")),
	apoe.e4 = factor(apoe.e4, levels = c(0, 1, 2), labels = c("Non-carrier", "One allele", "Two alleles")),
)
for(v in cov.Y){
	dat[[paste0("hist_", v)]] <- factor(if_else(!is.na(dat[[paste0("icd10Date_", v)]]) & dat[[paste0("icd10Date_", v)]] <= dat$last_date, 1L, 0L))
}

library(mice) # 协变量插补：5 datasets x 5 iterations
covs_raw <- dat %>% select(eid, all_of(covs)) 
covs_x <- covs_raw %>% select(-eid)
	ini  <- mice(covs_x, maxit = 0, print = FALSE)
	meth <- ini$method
	pred <- ini$predictorMatrix
	imp <- mice(covs_x, m = 5, maxit = 5, method = meth, predictorMatrix = pred, seed = 20251126)
covs_imp_body <- complete(imp, action = 1) # 取第 1 套插补结果作为最终协变量表
	covs_imp <- cbind(eid = covs_raw$eid, covs_imp_body)
saveRDS(covs_imp, "covs_imp.rds")

# outcome
dxs <- c(ad = "Alzheimer's Disease", acd = "All Cause Dementia", htn = "Hypertension", dm = "Diabetes", cvd = "Cerebrovascular Disease", cvd_o = "Other Cardiovascular Disease", depress = "Depression")
Ys <- names(dxs)
for (Y in Ys) {
	dat[grep(paste0("^", Y, "\\.([Y]?(t2e|r2e))$"), names(dat))] <- NULL
	dat <- t2e(dat, NA, paste0("icd10Date_", Y), "birth_date", "last_date", "date_lost", "date_death", "2022-12-31", Y, "day") 
}
Y <- "acd" 
dat1 <- dat %>% mutate(
	exit_date = pmin(date_lost, date_death, as.Date("2022-12-31"), na.rm = TRUE),	#计算随访结束日期（退出/死亡/行政截尾 中的最早者）
	Y_bl = !is.na(paste0("icd10Date_", Y)) & paste0("icd10Date_", Y) <= last_date
	) %>% 
	filter(is.na(exit_date) | last_date <= exit_date)  %>%  # 只保留 baseline 早于或等于随访结束的个体
	filter(!Y_bl)  # 只留 baseline 无 dementia 的人
table(dat1[[paste0(Y, ".Yt2e")]]) # ✅ 1,987 and 696 developed ACD and AD, mean ageonset 74.6 and 75.0

saveRDS(dat1, "dat1.rds")  # 182374


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Food_Dementia.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- readRDS("dat1.rds") 
dat_covs <- readRDS("covs_imp.rds")
field <- readxl::read_excel("foods.xlsx")

covs.basic <- 'energy.kJ age sex.c ethnic.b tdi edu.sco smoke.c pa_met_tot bmi apoe.e4' %>% strsplit(' ') %>% unlist()
covs <- c(covs.basic, "hist_htn", "hist_dm", "hist_cvd", "hist_cvd_o")
dat <- dat1 %>% select(-all_of(covs)) %>% left_join(dat_covs, by = "eid") # 换imp的协变量

Y <- "acd" 
covariate <- "covs"
	time_var <- paste0(Y, ".t2e")
	event_var <- paste0(Y, ".Yt2e")

dat1 <- dat %>% mutate(across(all_of(field$name), ~ifelse(. == 0, 0, 1))) # 把食物变量二分类化
output <- data.frame(exposure = character(), class = numeric(), outcome = character(), model = character(), total = integer(), case = integer(), HR = numeric(), LCI = numeric(), UCI = numeric(), P = numeric(), res.zph  = numeric())
for (p in field$name) {
	temp <- dat1 %>% select(eid, all_of(time_var), all_of(event_var), all_of(p), all_of(get(covariate)))
		colnames(temp)[4] <- "pheno"
		temp$pheno <- scale(temp$pheno)  
		temp <- na.omit(temp)  
	total <- temp %>% group_by(pheno) %>% summarise(n = n(), .groups = "drop")
		case <- temp %>% group_by(pheno) %>% summarise(n = sum(.data[[event_var]]), .groups = "drop")  
	fit <- coxph(as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ .")), data = temp[,-1])
	fit0 <- summary(fit)
	res <- cox.zph(fit)
		coef <- as.data.frame(fit0$coefficients)
		confint <- as.data.frame(fit0$conf.int)
	output0 <- data.frame(exposure = rep(p, 2), class = 1:2, outcome = rep(Y, 2), model = rep(covariate, 2), total = total$n[1:2], case = case$n[1:2], HR = c(1, coef$`exp(coef)`[1]), LCI = c(1, confint$`lower .95`[1]), UCI = c(1, confint$`upper .95`[1]), P = c(1, coef$`Pr(>|z|)`[1]), res.zph = res$table[1,3])
		output <- rbind(output, output0)
}
FDR <- function(p, n){row <- length(p); for (i in 1:n) {p[seq(i, row, n)] <- p.adjust(p[seq(i, row, n)],method = "fdr")}; p}
output$FDR <- FDR(output$P, 2)
write.csv(output, paste0("Food_", Y, "_", covariate, ".csv"), row.names = F)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Group_Dementia.R
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
func_quin <- function(x){ # 五分位
  breaks = quantile(x, probs = seq(0, 1, 1/5))
  if(length(unique(breaks)) == 2){a <- factor(ifelse(x > 0, 1, 0))
  }else if(length(unique(breaks)) < 6){a <- cut(x, breaks = c(min(x) - 1, unique(breaks)), include.lowest = T)
  }else{a <- cut(x, breaks = unique(breaks), include.lowest = T)}
  return(a)
}
dat1 <- dat %>% mutate(across(all_of(unique(field$Food.Group)), func_quin))

# 统计（各分位的中位数 / Range / n / case）
stat <- function(x, data_q, data_num){
	temp <- data_q %>% select(eid, all_of(time_var), all_of(event_var), quantile = all_of(x)) %>% left_join(data_num %>% select(eid, Food = all_of(x)), by = "eid")
	b <- temp %>% group_by(quantile) %>% 
		summarise(Median = round(median(Food), 2), Range = paste(round(range(Food), 2), collapse = "~"),  n = n(), case = sum(!!sym(event_var)), .groups = "drop") %>% 
		mutate(pheno = x, .before = "quantile") %>% 
		mutate(`Consumption level` = paste0(Median, " (", Range, ")"), .after = "quantile")
	return(b)
}
stat_data <- data.frame()
for (i in unique(field$Food.Group)) {
	stat0 <- stat(i, dat1, dat)
	stat_data <- rbind(stat_data, stat0)
}
write.csv(stat_data, "Statistic/Group_stat.csv", row.names = F)

# cox
dat1 <- mutate(dat1, across(all_of(unique(field$Food.Group)), as.integer))  # 把 quintile 因子转成整数
func_cox <- function(component, ref) {
	temp <- dat1 %>% select(eid, all_of(time_var), all_of(event_var), all_of(component), all_of(get(covariate))) %>% rename(pheno = all_of(component))
	lv <- sort(unique(temp$pheno[!is.na(temp$pheno)]))
	ref_use <- min(ref, max(lv))
	lv <- c(ref_use, lv[lv != ref_use])
	temp$pheno <- factor(temp$pheno, levels = lv)
	n <- length(lv)
	fit <- coxph(as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ .")), data = temp[, -1])
	s <- summary(fit); z <- cox.zph(fit); coef <- as.data.frame(s$coefficients)
	ci <- as.data.frame(s$conf.int)
	k <- n - 1L
	data.frame(Ref = rep(ref, n), exposure = rep(component, n), outcome = rep(Y, n), model = rep(covariate, n), class = levels(temp$pheno), HR = c(1, coef$`exp(coef)`[1:k]), LCI = c(1, ci$`lower .95`[1:k]), HCI = c(1, ci$`upper .95`[1:k]), P = c(1, coef$`Pr(>|z|)`[1:k]), res.zph = rep(z$table[1, 3], n))
}
output <- data.frame()

for (ref in 1:5) {
	for (component in unique(field$Food.Group)) {
		output0 <- func_cox(component, ref)
		output <- rbind(output, output0)
	}
}
write.csv(output, paste0("Result/Group_", Y, "_", covariate, ".csv"), row.names = F)



========== 下面的还没弄完 🈲 =============
pacman::p_load(readxl, data.table, tidyverse, lubridate, purrr, survival)

dir0 = "D:"
indir = paste0(dir0, "/data/ukb/phe")
source(paste0(dir0, '/scripts/f/phe.f.R'))

setwd("D:/analysis/diet")
dat0 <- readRDS("dat1.rds") 
field <- readxl::read_excel("foods.xlsx")
dat <- dat0 %>% select(eid, all_of(unique(field$Food.Group)))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Pattern.R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ML assisted diet score
dat1 <- dat %>% mutate(`Berries and citrus fruits` = Berries + `Citrus fruits`) %>% 
	transmute( eid,
		`Green leafy vegetables` = ifelse(`Green leafy vegetables` > 0.25, 1, 0),
		`Olive oil` = ifelse(`Olive oil` > 0, 1, 0),
		`Berries and citrus fruits` = ifelse(`Berries and citrus fruits` > 0 & `Berries and citrus fruits` <= 2, 1, 0),
		Potatoes = ifelse(Potatoes > 0 & Potatoes <= 0.75, 1, 0),
		Eggs = ifelse(Eggs > 0 & Eggs <= 1, 1, 0),
		Poultry = ifelse(Poultry > 0 & Poultry <= 0.5, 1, 0),
		`Sweetened beverages` = ifelse(`Sweetened beverages` == 0, 1, 0)
	)
dat1$ML_score <- rowSums(dat1[ , -1])
saveRDS(dat1, "ML_score.rds")

# MIND score
data <- dat %>% select(eid, Berries, `Green leafy vegetables`, `Red/Orange vegetables`, `Other vegetables`, Legumes, Nuts, `Whole grains`, Poultry, `Fish and other seafood`, `Olive oil`, Wine, `Red meats`, `Processed meats`, `Butter and margarine`,`Sweets and desserts`,`Fried and fast foods`)
data <- data %>% mutate(`Red and processed meats` = `Red meats` + `Processed meats`, .after = Wine) %>% select(-`Red meats`, -`Processed meats`)
cheese <- dat0 %>% select(eid, all_of(c("102810", "102820", "102830", "102840", "102850", "102860", "102870", "102880", "102890", "102900", "102910")) %>% 
	mutate(Cheese = rowSums(across(-eid), na.rm = TRUE)) %>% select(eid, Cheese)
data <- data %>% left_join(cheese, by = "eid") %>% relocate(Cheese, .after = `Red and processed meats`) 
#顺序：eid, Berries, Green leafy vegetables, Red/Orange veg, Other veg, Legumes, Nuts, Whole grains, Poultry, Fish & seafood, Olive oil, Wine, Red+proc meats, Cheese, Sweets & desserts, Fried & fast foods
func_ter <- function(x){
	breaks = quantile(x,probs = seq(0,1,1/3))
	if(length(unique(breaks)) == 2){a <- factor(ifelse(x > 0, 1, 0))
	}else if(length(unique(breaks)) < 4){a <- cut(x, breaks = c(min(x)-1, unique(breaks)), include.lowest = T)
	}else{a <- cut(x,breaks = unique(breaks),include.lowest = T)}
	return(a)
}
data1 <- mutate(data, across(.cols = 2:16, func_ter)) #三分位分组
data2 <- mutate(data1, across(.cols = 2:16, as.integer))
func_for1 <- function(x){a <- ifelse(x == 1, 0, 1); return(a)}
func_for2 <- function(x){a <- ifelse(x == 1, 0,ifelse(x == 2,0.5,1)); return(a)}
func_aginst1 <- function(x){a <- ifelse(x == 1, 1, 0); return(a)} # 
func_aginst2 <- function(x){a <- ifelse(x == 1, 1, ifelse(x == 2, 0.5, 0)); return(a)}
func_wine <- function(x){a <- ifelse(x == 1, 1, ifelse(x == 0 | x > 1, 0, 0.5)); return(a)}
data3 <- data2 %>% 
	mutate(across(.cols = c(2, 6, 10), func_for1)) %>%
	mutate(across(.cols = c(3, 4, 5, 7, 8, 9, 11), func_for2)) %>%
	mutate(across(.cols = c(12, 13, 14, 15, 16), func_aginst2))
data3 <- data3 %>% mutate(Wine = func_wine(data$Wine))
data3 <- data3 %>% mutate(MIND_score = rowSums(across(-eid), na.rm = TRUE))
saveRDS(data3, "MIND_score.rds")

####
# 有益食物：简单 for1
cols_for1  <- c("Berries", "Whole grains", "Olive oil")
# 有益食物：for2（三档）
cols_for2  <- c("Green leafy vegetables", "Red/Orange vegetables", "Other vegetables", "Legumes", "Nuts", "Poultry", "Fish and other seafood")
# 有害食物：against2
cols_again2 <- c("Red and processed meats", "Cheese", "Sweets and desserts", "Fried and fast foods")
data3 <- data2 %>%
  mutate(across(all_of(cols_for1),  func_for1)) %>%
  mutate(across(all_of(cols_for2),  func_for2)) %>%
  mutate(across(all_of(cols_again2), func_aginst2))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Pattern_Dementia.R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##Associations between dietary patterns and incident dementia
library(survival)
ML_score <- fread("ML_score.csv")
MIND_score <- fread("MIND_score.csv")
covariate <- fread("covariate_mice.csv")
outcome <- fread("Dementia.csv")
data <- merge(ML_score,MIND_score,by = "eid",sort = F)
data <- merge(data,covariate,by = "eid",sort = F)
data <- merge(outcome,data,by = "eid",sort = F)
#Multivariate Cox proportional hazards model one
func_quan <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
  setnames(data,c(out,exposure),c("outcome","pheno"))
  data$pheno <- cut(data$pheno,breaks = quantile(data$pheno),include.lowest = T)
  total <- group_by(data,pheno) %>% summarise(n=n())
  case <- group_by(data,pheno) %>% summarise(n=sum(outcome))
  fit <- coxph(Surv(days,outcome) ~ .,data = data)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  n=2;a=n+2;b=n+1
  output <- data.frame(exposure=rep(exposure,a),
                       outcome=rep(out,a),
                       model=rep(paste0("model_",model),a),
                       class=c(1:a),
                       class_value=total$pheno[1:a],
                       total=total$n[1:a],
                       case=case$n[1:a],
                       HR=c(1,coef$`exp(coef)`[1:b]),
                       LCI=c(1,confint$`lower .95`[1:b]),
                       HCI=c(1,confint$`upper .95`[1:b]),
                       P=c(1,coef$`Pr(>|z|)`[1:b]))
  return(output)
}
func_linear <- function(exposure,out,model){
  if(model == 1){data <- select(data,days,all_of(out),all_of(exposure),energy:APOE4)
  }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
  setnames(data,c(out,exposure),c("outcome","pheno"))
  # data$pheno <- scale(data$pheno,scale = 0.2*7)
  data$pheno <- (data$pheno-mean(data$pheno))/1.4
  n <- nrow(data);case <- sum(data$outcome,na.rm = T)
  fit <- coxph(Surv(days,outcome) ~ .,data = data)
  fit0 <- summary(fit)
  coef <- as.data.frame(fit0$coefficients)
  confint <- as.data.frame(fit0$conf.int)
  output <- data.frame(exposure=exposure,
                       outcome=out,
                       model=paste0("model_",model),
                       class=5,
                       class_value="Per 20% increment",
                       total=n,
                       case=case,
                       HR=c(coef$`exp(coef)`[1]),
                       LCI=c(confint$`lower .95`[1]),
                       HCI=c(confint$`upper .95`[1]),
                       P=c(coef$`Pr(>|z|)`[1]))
  return(output)
}

exposure <- c("ML_score","MIND_score")
out <- c("ACD_status","AD_status")
model <- 1:2
output <- data.frame(exposure=character(),outcome=character(),model=as.character(),class=numeric(),class_value=factor(),
                     total=integer(),case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
for (i in exposure) {
  for (j in out) {
    for (k in model) {
      output0 <- func_quan(i,j,k)
      output <- rbind(output,output0)
      output0 <- func_linear(i,j,k)
      output <- rbind(output,output0)
    }
  }
}
write.csv(output,"Result/Pattern_dementia.csv",row.names = F)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 Pattern_chronic_disease.R 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##Associations between ML_score and other health-related outcomes
library(survival)
ML_score <- fread("ML_score.csv")
covariate <- fread("covariate_mice.csv")
data <- merge(ML_score,covariate,by = "eid",sort = F)
output <- data.frame(exposure=character(),outcome=character(),model=as.character(),class=numeric(),class_value=factor(),
                     total=integer(),case=integer(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
output2 <- data.frame(exposure=character(),outcome=character(),model=as.character(),HR=numeric(),LCI=numeric(),UCI=numeric(),P=numeric())
start_date <- read.csv("UKB_dietary_date_max.csv",header = T)
start_date$start_date <- as.Date(start_date$start_date)
death <- read.csv("UKB_death.csv",header = T)
death$death <- pmax(death$X40000.0.0,death$X40000.1.0,na.rm = T)
end_date <- mutate(death,end_date=as.Date(ifelse(is.na(death),"2022-12-31",death))) %>% dplyr::select("eid","end_date")
outco <- merge(start_date,end_date,by="eid")
name <- fread("Target_code.csv")
for (out_name in name$Disease_code) {
  outcome <- read.csv(paste0("chronic/Targets_RAW/",out_name,".csv"),header = T)[,c("eid","target_date","target_y")]
  outcome <- merge(outco,outcome,by="eid")
  outcome <- dplyr::select(outcome,eid,start_date,end_date,target_date) %>% mutate(target_date=as.Date(target_date)) %>%
    mutate(days=ifelse(is.na(.data$target_date),difftime(end_date,start_date,units = "days"),difftime(.data$target_date,start_date,units = "days")),
           status=ifelse(is.na(.data$target_date),0,1),.keep="unused",.after=eid)
  data1 <- merge(outcome,data,by = "eid",sort = F)
  data1 <- filter(data1,days > 0)
  exposure <- "ML_score"
  out <- "status"
  func_quan <- function(exposure,out,model,input){
    data <- select(input,days,all_of(out),all_of(exposure),energy:BMI)
    setnames(data,c(out,exposure),c("outcome","pheno"))
    data$pheno <- cut(data$pheno,breaks = quantile(data$pheno),include.lowest = T)
    total <- group_by(data,pheno) %>% summarise(n=n())
    case <- group_by(data,pheno) %>% summarise(n=sum(outcome))
    fit <- coxph(Surv(days,outcome) ~ .,data = data)
    fit0 <- summary(fit)
    coef <- as.data.frame(fit0$coefficients)
    confint <- as.data.frame(fit0$conf.int)
    n=2;a=n+2;b=n+1
    output <- data.frame(exposure=rep(exposure,a),
                         outcome=rep(name$Disease[name$Disease_code == out_name],a),
                         model=rep(paste0("model_",model),a),
                         class=c(1:a),
                         class_value=total$pheno[1:a],
                         total=total$n[1:a],
                         case=case$n[1:a],
                         HR=c(1,coef$`exp(coef)`[1:b]),
                         LCI=c(1,confint$`lower .95`[1:b]),
                         HCI=c(1,confint$`upper .95`[1:b]),
                         P=c(1,coef$`Pr(>|z|)`[1:b]))
    return(output)
  }
  for (i in exposure) {
    for (j in out) {
      for (k in model) {
        output0 <- func_quan(i,j,k,dat)
        output <- rbind(output,output0)
      }
    }
  }
  
  func_linear <- function(exposure,out,model,input){
    if(model == 1){data <- select(input,days,all_of(out),all_of(exposure),energy:BMI)
    }else if(model == 2){data <- select(data,days,all_of(out),all_of(exposure),energy:CVD)}
    setnames(data,c(out,exposure),c("outcome","pheno"))
    data$pheno <- scale(data$pheno,scale = 0.2*7)
    n <- nrow(data);case <- sum(data$outcome,na.rm = T)
    fit <- coxph(Surv(days,outcome) ~ .,data = data)
    fit0 <- summary(fit)
    coef <- as.data.frame(fit0$coefficients)
    confint <- as.data.frame(fit0$conf.int)
    output <- data.frame(exposure=exposure,
                         outcome=name$Disease[name$Disease_code == out_name],
                         model=paste0("model_",model),
                         class=5,
                         class_value="Per 20% increment",
                         total=n,
                         case=case,
                         HR=c(coef$`exp(coef)`[1]),
                         LCI=c(confint$`lower .95`[1]),
                         HCI=c(confint$`upper .95`[1]),
                         P=c(coef$`Pr(>|z|)`[1]))
    return(output)
  }
  for (i in exposure) {
    for (j in out) {
      for (k in model) {
        output0 <- func_linear(i,j,k,dat)
        output <- rbind(output,output0)
      }
    }
  }
}
output$FDR <- FDR(output$P,5)
write.csv(output,paste0("pattern_result_linear_New","",".csv"),row.names = F)
write.csv(output2,paste0("pattern_result_New","",".csv"),row.names = F)