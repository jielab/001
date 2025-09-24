pacman::p_load(data.table, tidyverse, lubridate, survival, rcssci)

dir0 = 'D:'
source(paste0(dir0, '/scripts/f/phe.f.R'))
le8.vars = c("diet", "pa", "smoke", "bmi", "nonhdl", "hba1c", "bp", "sleep")

dat0 <- readRDS(file = paste0(dir0, '/data/ukb/phe/Rdata/all.rds')) %>% mutate (
	C677T = 2 - folate.mthfr.rs1801133_G, 
	A1298C = 2 - folate.mthfr.rs1801131_T
)
snp = "C677T" # cyp2c19.rs12248560_C cyp2c19.rs4986893_G  cyp2c19.rs4244285_G # 不同人种间的差别 👨‍👩‍👧‍👦
	snp; dat0$snp = round(dat0[[snp]])
	round(prop.table(table(dat0$ethnic.c, dat0$snp), 1), 2)
	dat0 %>% group_by(ethnic.c) %>% summarise(count = n(), af = mean(snp, na.rm = TRUE) / 2) 

dat <- dat0 %>% filter(ethnic.c == "White", gen_ethnicity==1) # sex == 1, prot.yes == 1 
sum(!is.na(dat$icdDate_cmd))
	dat %>% group_by(cvh.c) %>% summarise(count = n(), mean = mean(le8.sco, na.rm = TRUE))
	hist(dat$date_attend, breaks = 'months', freq = TRUE)
	print(date_attend.med <- median(dat$date_attend, na.rm = TRUE))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🧠 以 stroke 为主的 Cox 分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
grep("death", names(dat), value = TRUE); table(substr(dat$death_icd10, 0, 1))
sum(!is.na(dat$icdDate_cad))
Y = "cad"; X = "le8.sco"; X1 = "cvh.c"; X2 = "cad.2m.score_sum" # sleep.pts
dat1 <- dat %>% mutate(
	X = dat[[X]], X1 = dat[[X1]], X2 = dat[[X2]],
	X2 = cut(X2, breaks = quantile(X2, probs = seq(0, 1, 0.2), na.rm = T), include.lowest = T, labels = paste0('q', 1:5)),
	X2 = factor(ifelse(X2 == 'q1', 'low', ifelse(X2 == 'q5', 'high', 'middle')), levels = c('low', 'middle', 'high')),
	Y_date = dat[[paste0("icdDate_", Y)]],
	Y_age = time_length(interval(birth_date, Y_date), unit = "years"),
	Y_yes = ifelse(is.na(Y_date), 0, 1), 
	follow_date = ifelse(!is.na(Y_date), Y_date, ifelse(!is.na(date_lost), date_lost, ifelse(!is.na(date_death), date_death, as.Date('2021-12-31')))),
	follow_years = (as.numeric(follow_date) - as.numeric(date_attend)) / 365.25,
	before_date = ifelse(!is.na(Y_date), Y_date, ifelse(!is.na(birth_date), birth_date, NA)),
	before_years = (as.numeric(date_attend) - as.numeric(before_date)) / 365.25
)
dat1 %>% group_by(sex.c, age.c) %>% summarise(n = n(), n_case = sum(Y_yes, na.rm = TRUE), freq = mean(Y_yes, na.rm = TRUE), .groups = "drop")

prop.table(table(dat1$bald, dat1$Y_yes), 1) # 看分布🌈，线性、非线性
	rcssci_cox(data = dat1, time = 'follow_years', y = 'Y_yes', x = 'X2', covs = c('age', 'bald'), prob = 0.1, filepath = 'D:/tmp') # 🐂
	summary(glm(Y_yes ~ age + sex + bald, data = dat1, family = "binomial"))
	subset(dat1, Y_date > date_attend) %>% nrow(); summary(dat1$date_attend)
	# datmp <- dat1 %>% dplyr::select(ethnicity, age, sex, bmi, date_attend, birth_year, birth_month, Y_date, matches("^prot_"))
	# fwrite(datmp, file = "D:/dat1.txt", append = FALSE, sep = "\t", na = NA, row.names = FALSE, quote = FALSE)

par(mar = c(5, 4, 4, 5) + 0.1) # 看分布，前、后🥚🐓❓
	myhist <- hist(dat1$Y_date, breaks = "years", freq = TRUE, main = "X", ylab = "Y (freq)")
	X.avgs <- by(dat1$X, cut(dat1$Y_date, breaks = as.Date(myhist$breaks)), function(x) mean(x, na.rm = TRUE)) # 用 x == x_val 求频率🏮
	X.sds <- by(dat1$X, cut(dat1$Y_date, breaks = as.Date(myhist$breaks)), function(x) sd(x, na.rm = TRUE))
	par(new = T)
	plot(as.Date(myhist$mids), X.avgs, xlim = range(as.Date(myhist$breaks)), ylim = range_3d(dat1$X), pch = 16, axes = FALSE, xlab = NA, ylab = NA, cex = 1.2, col = "red")
	arrows(as.Date(myhist$mids), X.avgs-X.sds, as.Date(myhist$mids), X.avgs+X.sds, angle = 90, code = 3, length = 0.05, col = "red")
	axis(side = 4); mtext(side = 4, line = 3, "X", col = "red")
	abline(v = date_attend.med, col = "blue", lwd = 3, lty = 2)
	abline(h = mean(dat1$X, na.rm = TRUE), col = "green", lwd = 3, lty = 2)

dat1 <- dat1 %>% filter(follow_years > 0) # 下面进行Survival🏊‍分析
	surv.obj <- Surv(time = dat1$Y_age, event = dat1$Y_yes) 
	survdiff(surv.obj ~ bald, data = dat1)
covs <- 'age bmi deprive PC1 PC2' %>% strsplit(' ') %>% unlist()
	fit <- coxph(as.formula(paste('surv.obj ~ X1 + X2 +', paste(covs, collapse = '+'))), data = dat1); summary(fit)
	survminer::ggforest(fit, main = '', cpositions = c(0, 0.1, 0.3), fontsize = 1.2, data = dat1)
	km.obj <- survfit(surv.obj ~ sex.c, data = dat1) # KM不能带协变量，M作为分层变量
	plot(km.obj, fun = function(x) 1-x)
	survminer::ggsurvplot(km.obj, fun = 'event', pval = TRUE, risk.table = TRUE, ncensor.plot = TRUE, ggtheme = theme_bw())

surv.obj <- Surv(time = dat1$follow_years, event = dat1$Y_yes) # ✝ 10年风险计算
	fit <- coxph(surv.obj ~ X1 + X2, data = dat1) 
	dat2 <- expand.grid(X1 = levels(dat1$X1), X2 = levels(dat1$X2))
	surv.pred <- survfit(fit, newdata = dat2)
	surv.10y <- summary(surv.pred, times = 10)
	dat2 <- dat2 %>% mutate(risk = 1 - t(surv.10y$surv), ci_lower = t(1 - surv.10y$upper), ci_upper = t(1 - surv.10y$lower))
	# library(nricens); dat1$risk <- nricens::get.risk.coxph(fit, 10) #  ✝ 10年风险计算简化版 🏮
	# dat2 <- dat1 %>% group_by(X2, X1) %>% summarise(n = n(), mean_risk = mean(risk), sd_risk = sd(risk), se = sd_risk / sqrt(n), ci_lower = mean_risk - qt(0.975, df = n - 1) * se, ci_upper = mean_risk + qt(0.975, df = n - 1) * se, .groups = "drop") %>% rename(risk = mean_risk)
	ggplot(dat2, aes(x = X2, y = risk, fill = X1)) + geom_bar(stat = 'identity', position = position_dodge(width = 0.7), width = 0.7) +
		geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, position = position_dodge(width = 0.7)) +
		geom_text(aes(label = sprintf('%.2f', risk), y = ci_upper), vjust = -1, position = position_dodge(width = 0.7), size = 2) +
		labs(x = '', y = '10-year risk (%)', title = '') + scale_fill_manual(values = rainbow(length(levels(dat1$X1))), name = '') + theme_minimal() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 💇‍ 以 bald 为主的表型分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
# 表一🪴🦂
dat %>% filter(sex == 1) %>% group_by(bald) %>% summarise(count = n(), mean = mean(bald.score_sum, na.rm = TRUE))
covs <- "age bmi X stroke.c sex_first.c sex_par.c bb_TES bb_SHBG PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()
	library(compareGroups); descrTable(formula = as.formula(paste(c("bald ~ ", paste(covs, collapse = "+")), collapse = "")), data = dat, digits = )
	library(gtsummary); dat %>% select(all_of(c('bald', covs))) %>% tbl_summary(by = bald, missing = 'no') %>% add_p()

dat1 <- dat %>% mutate(Y = bald, Y2 = stroke.c) %>% drop_na(X, Y, Y2) # 看某个❌的分布 
	group_by(dat1, Y, Y2) %>% summarise(count = n(), mean = mean(X, na.rm = TRUE))
	aggregate(X ~ Y+Y2, dat1, FUN = function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs = c(0,0.5,1))), 3)} ) # Y*C
	bp <- boxplot(X ~ Y+Y2, dat1, las = 2, col = rainbow(4), font = 2); bp$stats

# 尝试不同的回归方法☯
	fit.lm <- lm(formula = as.formula(paste(c('bald ~ ', paste(covs, collapse = '+')), collapse = '')), data = dat)
	fit.ordi <- MASS::polr(formula = as.formula(paste(c('bald.ord ~ ', paste(covs, collapse = '+')), collapse = '')), data = dat, Hess = TRUE); 
	fit.mnom <- nnet::multinom(formula = as.formula(paste(c('bald.cat ~ ', paste(covs, collapse = '+')), collapse = '')), data = dat)
	print(fit.lm.sum <- summary(fit.lm)); AIC(fit.lm)
	print(fit.ordi.sum <- summary(fit.ordi)); AIC(fit.ordi); cbind(fit.ordi.sum$coef, 2 * (1 - pnorm(abs(fit.ordi.sum$coef[,1] / fit.ordi.sum$coef[,2]))))
	print(fit.mnom.sum <- summary(fit.mnom)); AIC(fit.mnom); fit.mnom.sum$coefficients; fit.mnom.sum$standard.errors

plots <- list() # 多个GLM回归结果展示🎄
traits <- c("bald12", "bald13", "bald14"); colors <- c("red", "blue", "green") 
for (i in 1:3) {
	Y <- traits[i]
	fit <- glm(as.formula(paste(Y, '~', paste(covs, collapse = '+'))), data = dat1, family = 'binomial')
	results <- broom::tidy(fit, conf.int = TRUE) %>% filter(term != "(Intercept)") %>% mutate(term = factor(term, levels = rev(unique(term))), OR = exp(estimate), OR_CI_low = exp(conf.low), OR_CI_high = exp(conf.high), p_text = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value)))
	if (nrow(results) > 4) { results <- results[1:(nrow(results) - 4), ]} # 去掉最后的PC1-PC4
	p <- ggplot(results, aes(x = OR, y = term)) + geom_pointrange(aes(xmin = OR_CI_low, xmax = OR_CI_high), shape = 24, fill = colors[i], color = "black") +
		geom_vline(xintercept = 1, linetype = "dashed") + labs(x = "Odds Ratio (95% CI)", y = if (i == 1) "Exposure" else NULL) + ggtitle(Y) + theme_minimal() + geom_text(aes(label = paste0("OR = ", round(OR, 2), "\nP = ", p_text)), hjust = -0.2, size = 3) 
	if (i > 1) { p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) }
	plots[[i]] <- p; # plots[[i]] <- forestmodel::forest_model(fit)
}
library(patchwork); wrap_plots(plots, ncol = 3)

x = "prot_eda2r"; y= "bald" # 某个回归结果的深挖🔍
	dat1 <- dat %>% mutate(X = dat[[x]], Y = dat[[y]]) %>% drop_na(X, Y) 
	dat1 %>% group_by(Y) %>% summarise(count = n(), mean = mean(X, na.rm = TRUE))
	library(ggpubr); ggboxplot(dat1, x = y, y = x, add = 'mean_se', color = y) + 
		stat_compare_means(method = 'anova', label = 'p.format', label.x = 1.5) +
		stat_compare_means(aes(label = paste0('p = ', round(after_stat(p), 3))), method = 't.test', ref.group = '1', comparisons = list(c('1', '2'), c('1', '3'), c('1', '4'))) +
		stat_summary(fun = 'mean', geom = 'text', aes(label = round(after_stat(y), 2)), vjust = -1.5, color = 'black') +theme_minimal()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 众多因素 vivid 交互作用分析 💃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- head(dat, 10000)
	library(GGally); ggpairs(dat1, mapping = aes(color = as.factor(bald)), columns = c('age', 'bmi', 'stiffness'))
	ggpairs(subset(dat1, select = c('age', 'bmi', 'stiffness', 'bald')), upper = list(continuous = 'density', combo = 'box_no_facet'), lower = list(continuous = 'points', combo = 'dot_no_facet'))

library(ComplexHeatmap); library(circlize)
	Ys <- c('asthma', 'breast_cancer', 'cad', 'cancer', 'circulatory', 'copd', 'covid', 'lung_cancer', 'mi', 'prostate_cancer', 'ra', 'stroke', 't1dm', 't2dm', 'varicose')
	Xs <- grep('bb_', names(dat), value = TRUE)[1:10] 
	icdDate_Ys <- paste0('icdDate_', Ys); setdiff(icdDate_Ys, names(dat)); t2_Ys <- paste0('t2_', Ys)
	dat1 <- dat %>% mutate(across(.cols = all_of(icdDate_Ys), .fns = ~ as.numeric(. - date_attend), .names = 't2_{col}')) %>% 
		rename_with(~ gsub('icdDate_', '', .), .cols = starts_with('t2_icdDate_')) %>% 
		dplyr::select(all_of(c(Xs, t2_Ys)))
	cor_matrix <- sapply(t2_Ys, function(y) { 
		sapply(Xs, function(x) { 
			cor(dat1[[x]], dat1[[y]], method = "spearman", use = "pairwise.complete.obs") 
		})
	})
	Heatmap(cor_matrix, col = colorRamp2(c(-0.5, 0, 0.5), c('blue', 'white', 'red')), column_title = 'diseases onset', row_title = 'blood biochemistry', heatmap_legend_param = list(title = 'r2', at = c(-1, 0, 1), labels = c('-1', '0', '1')))

library(vivid) # 正式的vivid分析
	dat1$Y_yes <- as.factor(dat1$Y_yes); dat1$follow_years <- NULL
	fit.rforest <- randomForest(Y_yes ~ ., na.action = na.omit, data = dat1)
	# fit.rforest <- rfsrc(surv.obj ~ ., data = dat1)
	fit.vivi <- vivi(fit = fit.rforest, response = 'Y_yes', data = dat1); print(fit.vivi, digits = 1)
	pdf(paste(Y,X,'heatmap.pdf',sep = '.')); print(viviHeatmap(mat = fit.vivi)); dev.off()
	pdf(paste(Y,X,'network.pdf',sep = '.')); print(viviNetwork(mat = fit.vivi)); dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GWAS数据比较🧬
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(qqman, CMplot, topr, plotly) # PheGWAS
trait = 'ALDH5A1'
dir = "D:/data/gwas/ppp/clean"
dat <- read.table(paste0(dir, '/', trait, '.1e-03'), header = TRUE) %>% 
	filter(P < 1e-03, EAF > 0.001 & EAF < 0.999) %>% 
	mutate(CHR = ifelse(CHR == "X",23,CHR), CHR = as.numeric(CHR), P = ifelse(P<1e-50,1e-50,P))
	qqman::manhattan(dat, main = trait, chr = "CHR", bp = "POS", p = "P", snp = "SNP", col = c("blue4", "orange3"))
	
labels = c("bald12", "bald13", "bald14", "bald1234")
for (label in labels) {
	datmp <- read.table(paste0(dir, '/', label, '.1e-3'), header = TRUE) %>%
		filter(P < 1e-03, P != 0) %>% na.omit() %>% 
		dplyr::select(SNP, CHR, POS, EA, NEA, BETA, SE, P) %>%
		mutate(P = ifelse(P < 1e-50, 1e-50, P)) 
	print(paste(label, nrow(datmp), ncol(datmp)))
	eval(parse(text = paste0('dat_', label, '.cmpt <- datmp %>% select(SNP, CHR, POS, P) %>% rename(P.', label, ' = P)') ))
	assign(paste0('dat_', label, '.topr'), datmp %>% rename(REF = NEA, ALT = EA))
}
dat_cmpt <- Reduce(function(x, y) merge(x, y, by = "SNP", all = TRUE), lapply(labels, function(label) get(paste0("dat_", label, ".cmpt"))))
dat_cmpt$CHR <- apply(subset(dat_cmpt, select = grepl("CHR", names(dat_cmpt))), 1, FUN = min, na.rm = TRUE)
dat_cmpt$POS <- apply(subset(dat_cmpt, select = grepl("POS", names(dat_cmpt))), 1, FUN = min, na.rm = TRUE)
dat_cmpt <- subset(dat_cmpt, select = c('SNP','CHR','POS', grep('^P\\.', names(dat_cmpt),value = TRUE)))
CMplot(dat_cmpt, plot.type = "m", multracks = TRUE, cex = 0.2, amplify = FALSE, file.output = TRUE, file = "jpg", file.name = "", width = 20, height = 4, dpi = 300) 
dat_topr <- mget(paste0('dat_', labels, '.topr')) # 🏮
png("topr.png", w = 2000, h = 1000, res = 128); topr::manhattan(dat_topr, ntop = 2, legend_labels = labels, color = c('grey','blue','orange','green')); dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 蛋白质cisMr数据挖掘🎂
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phes <- c('bald12', 'bald13', 'bald14', 'bald134', 'bald1234')
anas <- c('phe', 'clp08', 'cml08', 'g1', 'g2')
dat.phe <- read.table(paste0(dir0, '/analysis/assoc.sum/all.assoc.txt'), header = TRUE) %>%
	mutate(X = ifelse(X %like% 'prot_', gsub('PROT_', 'prot_', toupper(X)), X), Analysis = 'phe') 
	Y = 'bald12'; valcano(Y, dat.phe[dat.phe$Y == Y & dat.phe$X %like% 'prot_', ], 'X', 'BETA', 'P', 0.05, '', '') # 👀
dat.mr <- read.table(paste0(dir0, '/analysis/assoc.sum/all.mr.log'), header = TRUE)[,c(1:4,6:8)] %>% rename(P = P.ivw) %>% 
	mutate(X_new = ifelse(Analysis %like% 'way2', Y, X), Y_new = ifelse(Analysis %like% 'way2', X, Y), X = X_new, Y = Y_new) %>% dplyr::select(-X_new, -Y_new)
dat.cis <- read.table(paste0(dir0, '/analysis/assoc.sum/all.cisMr.log'), header = TRUE)[,-5] 
dat <- rbind(dat.mr, dat.cis) %>% mutate(X = paste0('prot_', X), Y = gsub('pancre', 'pancre_cancer', Y)) %>%
	mutate(Analysis = paste(Analysis, sprintf('%02d', p_t), sep = '.')) %>% dplyr::select(-p_t) %>% 
	mutate(Analysis = recode(Analysis, 'gwIV.way1.08' = 'g1', 'gwIV.way2.08' = 'g2', 'cis.clump.06' = 'clp06', 'cis.clump.08' = 'clp08', 'Mr.08' = 'cml08', 'Mr.06' = 'cml06'))
dat <- dat[, colnames(dat.phe)] # 🏮
dat <- rbind(dat, dat.phe) 
dat <- dat %>% 
	pivot_wider(names_from = c(Y, Analysis), values_from = c(BETA, SE, P), names_glue = '{Y}.{.value}.{Analysis}'); names(dat)
	dat <- dat %>% select(-matches('06$|SE')) 
	cols <- apply(expand.grid(c('BETA','P'), anas, phes)[,c(3,1,2)], 1, function(x) paste(x, collapse = '.'))
	cols <- c('X', cols); dat <- dat[, cols] # 🏮
dat.coloc <- read.table(paste0(dir0, '/analysis/assoc.sum/all.coloc.log'), header = TRUE)[,c(1,2,8)] %>%
	mutate(X = paste0('prot_',X)) %>% 
	pivot_wider(names_from = Y, values_from = H4, names_glue = '{Y}.{.value}'); names(dat.coloc)
dat <- merge(dat, dat.coloc, by = 'X', all = TRUE) 
	valcano('bald12', dat, 'X', 'bald12.BETA.g1', 'bald12.P.g1', 0.05, '', '') # 👀
	for (i in phes) { for (j in anas) {
		assign(paste0('plot.', i,'.',j), valcano(i, dat, 'X', paste0(i,'.BETA.',j), paste0(i,'.P.',j), 0.05, '', '')) # 最后两个空白表示 label_x, label_y
	}}
	plots <- mget( outer(phes, anas, function(i, j) paste0('plot.', i, '.', j)) |> as.vector() )
	gridExtra::grid.arrange(grobs = plots, ncol = length(phes), nrow = length(anas))
	bbplot('bald12 Clump^cML BETA', dat, 'X', 'bald12.BETA.clp08', 'bald12.BETA.cml08', 'yes', 'bald12.P.clp08', 'bald12.P.cml08', 'no', 0.05)	
	bbplot('bald12 Clump^cML P', dat, 'X', 'bald12.BETA.clp08', 'bald12.BETA.cml08', 'no', 'bald12.P.clp08', 'bald12.P.cml08', 'yes', 0.05)	

for (Y in Ys) { # 画火山🌋图
	res1 <- read.table(paste0(Y,'.assoc.txt'), header = TRUE) %>% filter(X %like% 'prot_') %>%
		pivot_wider(names_from = Analysis, values_from = c(BETA, SE, P), names_glue = '{.value}.{Analysis}'); names(res1)
		pdf(paste0(Y,'.pdf'))
		valcano(Y, res1[res1$X %like% 'prot_', ], 'X', 'BETA.cox.a', 'P.cox.a', 0.05, '', '')
		valcano(Y, res1[res1$X %like% 'prot_', ], 'X', 'BETA.cox.b', 'P.cox.b', 0.05, '', '')
		dev.off()
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 💊药物靶点Prioritize⭕️
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppp <- read.table(paste0(dir0, '/data/gwas/ppp/map_3k_v1.tsv'), sep = '\t', quote = '', header = TRUE, fill = TRUE)[,c(1,3,5)]
drug <- read.table(paste0(dir0, '/files/prot-drug.txt'), sep = '\t', header = TRUE, fill = TRUE)
drug <- merge(ppp, drug, by.x = 'UniProt', by.y = 'UNIPROT_ACCESSION') %>% 
	mutate(drug = ifelse(APPROVED_DRUG_TARG_CONF != '', 'approved', ifelse(CLINICAL_DRUG_TARG_CONF != '', 'clinical', ifelse(DRUG_POTENTIAL_TARGET != '', 'potential', NA)))) %>%
	select(ID, olink_target_fullname, drug) %>% mutate(ID = paste0('prot_', ID)) %>% arrange('drug')
	names(drug) <- c('X', 'olink_target', 'drug')
	drug[duplicated(drug$X) | duplicated(drug$X, fromLast = TRUE), ]
	drug <- drug[!duplicated(drug$X), ]
dat <- merge(dat, drug, by = 'X', all.x = TRUE)
saveRDS(dat, "bald.rds")
write.table(dat, 'bald.tsv', sep = '\t', append = FALSE, quote = FALSE, row.names = FALSE)