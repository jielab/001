pacman::p_load(writexl, tidyverse, survival, flextable, patchwork, gtsummary, broom, cowplot)

dir0 <- ifelse(Sys.info()[["sysname"]] == "Windows", "D:", "/work/sph-huangj")
indir <- paste0(dir0, "/data/ukb/phe")
invisible(lapply(c("phe.f.R", "assoc.f.R", "plot.f.R", "pred.f.R"), \(f) source(file.path(dir0, "/scripts", "f", f))))
setwd2(paste0(dir0, "/analysis/maha"))

log_file <- "maha.log.txt"
if (file.exists(log_file)) file.remove(log_file)
log_con <- file(log_file, open = "wt")
sink(log_con, split = TRUE)
on.exit({ while (sink.number() > 0) sink(); close(log_con) }, add = TRUE)
options(width = 200, warn = 1)
step_header <- function(x) cat("\n", strrep("=", 88), "\n", x, "\n", strrep("=", 88), "\n", sep = "")
save_xlsx1 <- function(x, file) writexl::write_xlsx(if (is.data.frame(x)) list(data = x) else x, file)
save_plot <- function(p, file, width, height, dpi = 320, bg = "white") ggsave(file, p, width = width, height = height, dpi = dpi, bg = bg)

center.lst <- c(
	"c10003" = "Stockport", "c11001" = "Manchester", "c11002" = "Oxford", "c11003" = "Cardiff", "c11004" = "Glasgow",
	"c11005" = "Edinburgh", "c11006" = "Stoke", "c11007" = "Reading", "c11008" = "Bury", "c11009" = "Newcastle",
	"c11010" = "Leeds", "c11011" = "Bristol", "c11012" = "Barts", "c11013" = "Nottingham", "c11014" = "Sheffield",
	"c11016" = "Liverpool", "c11017" = "Middlesbrough", "c11018" = "Hounslow", "c11020" = "Croydon", "c11021" = "Birmingham",
	"c11022" = "Swansea", "c11023" = "Wrexham", "c11024" = "Cheadle", "c11025" = "Cheadle", "c11026" = "Reading",
	"c11027" = "Newcastle", "c11028" = "Bristol"
)

dx.lst <- c(
	cvd_htn = "Hypertension", cvd_cad = "Coronary artery disease", mi = "Myocardial infarction",
	cvd_stroke_i = "Ischemic stroke", cvd_stroke_ih = "Intracerebral hemorrhage", cvd_stroke_sh = "Subarachnoid hemorrhage",
	cvd_hfail = "Heart failure", cvd_afib = "Atrial fibrillation", cvd_pad = "Peripheral arterial disease",
	cvd_dvt = "Deep vein thrombosis", cvd_pe = "Pulmonary embolism", t2dm = "Type 2 diabetes",
	acd = "All-cause dementia", ad = "Alzheimer's disease", depress = "Depression",
	copd = "Chronic obstructive pulmonary disease", asthma = "Asthma",
	ckd = "Chronic kidney disease", ra = "Rheumatoid arthritis", death = "All-cause death"
)

diet.lst <- c(medito = "MEDI-Touch", medi24 = "MEDI", dash = "DASH", mind = "MIND", hpdi = "hPDI", ahei = "AHEI", phdi = "PHDI", modern = "MODERN", digm = "DIGM", maha = "MAHA")
diet.sum.cols <- paste0("diet.", names(diet.lst), ".sum")
diet.order <- unname(diet.lst)
diet.pts.cols <- paste0("diet.", names(diet.lst), ".pts")
diet.s100.cols <- paste0("diet.", names(diet.lst), ".s100")
diet.q5.cols <- paste0("diet.", names(diet.lst), ".q5")
diet.s100_q5.cols <- paste0("diet.", names(diet.lst), ".s100_q5")

vars.basic <- c("age", "sex", "center", "tdi", "PC1", "PC2")
vars.le8 <- c("diet.pts", "pa.pts", "smoke.pts", "bmi.pts", "nonhdl.pts", "hba1c.pts", "bp.pts", "sleep.pts")
names.le8 <- gsub("\\.pts$", "", vars.le8)

covs <- c("age", "sex.f", "center", "tdi", "PC1", "PC2", "bmi.pts", "bp.pts", "nonhdl.pts", "smoke.pts", "pa.pts", "sleep.pts")
diet.inc <- c("maha", "dash", "mind", "medi24")
diet.inc.pts <- paste0("diet.", diet.inc, ".pts")
diet.inc.s100 <- paste0("diet.", diet.inc, ".s100")
diet.inc.q5 <- paste0("diet.", diet.inc, ".q5")
diet.inc.s100_q5 <- paste0("diet.", diet.inc, ".s100_q5")

Y.inc <- c("cvd_cad", "cvd_stroke_i", "cvd_hfail", "t2dm", "ckd", "death")
Y.cols <- setNames(c("#B97AF7", "#F26D60", "#D98C2B", "#7CAE00", "#16B9C0", "#1F78B4"), dx.lst[Y.inc])

group_pct_th <- 0.40
mk_hml <- function(x, p = 0.40) { q1 <- quantile(x, probs = p, na.rm = TRUE, names = FALSE); q2 <- quantile(x, probs = 1 - p, na.rm = TRUE, names = FALSE); case_when(is.na(x) ~ NA_character_, x <= q1 ~ "low", x >= q2 ~ "high", TRUE ~ "middle") }
get_mode <- function(x) { x <- x[!is.na(x)]; if (length(x) == 0) NA else names(sort(table(x), decreasing = TRUE))[1] }
make_typical_value <- function(x) { x2 <- x[!is.na(x)]; if (length(x2) == 0) NA else if (is.numeric(x)) mean(x2) else if (is.factor(x)) factor(get_mode(x2), levels = levels(x)) else if (is.character(x)) get_mode(x2) else x2[1] }
var2lab <- function(x) { x0 <- gsub("^diet\\.|\\.(s100_q5|s100|sum|pts|q5|3c|hml)$", "", x); out <- unname(diet.lst[x0]); ifelse(is.na(out), x, out) }
zstd <- function(x) { x <- as.numeric(x); if (is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) return(rep(NA_real_, length(x))); as.numeric(scale(x)) }


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 数据及QC
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("数据及QC")
dat0 <- readRDS(paste0(indir, "/Rdata/all.rds"))

std_10_90 <- function(x) { x <- as.numeric(x); q <- quantile(x, c(0.1, 0.9), na.rm = TRUE, names = FALSE); if (!all(is.finite(q)) || q[1] == q[2]) rep(NA_real_, length(x)) else (x - q[1]) / (q[2] - q[1]) }
score_q5 <- function(x) { x <- as.numeric(x); out <- rep(NA_real_, length(x)); ok <- is.finite(x); out[ok] <- c(0, 25, 50, 75, 100)[dplyr::ntile(x[ok], 5)]; out }
score_0_100 <- function(x) { x <- as.numeric(x); ok <- is.finite(x); if (sum(ok) < 2 || diff(range(x[ok])) == 0) rep(NA_real_, length(x)) else { out <- rep(NA_real_, length(x)); out[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100; out } }
score_0_100_q5 <- function(x) { x <- as.numeric(x); ok <- is.finite(x); if (sum(ok) < 2 || diff(range(x[ok])) == 0) rep(NA_real_, length(x)) else { z <- rep(NA_real_, length(x)); z[ok] <- (x[ok] - min(x[ok])) / diff(range(x[ok])) * 100; as.numeric(as.character(cut(z, c(-Inf, 20, 40, 60, 80, Inf), labels = c(0, 25, 50, 75, 100), right = FALSE))) } }
domain_of <- function(Y) ifelse(grepl("^cvd_", Y) | Y == "mi", "cvd", NA)

dat <- dat0 %>%
	filter(ethnic.c == "White") %>%
	mutate(
		across(all_of(diet.sum.cols), std_10_90, .names = "{sub('\\\\.sum$', '', .col)}.pts"),
		across(all_of(diet.sum.cols), score_0_100, .names = "{sub('\\\\.sum$', '', .col)}.s100"),
		across(all_of(diet.sum.cols), score_q5, .names = "{sub('\\\\.sum$', '', .col)}.q5"),
		across(all_of(diet.sum.cols), score_0_100_q5, .names = "{sub('\\\\.sum$', '', .col)}.s100_q5"),
		across(all_of(diet.sum.cols), f3c, .names = "{sub('\\\\.sum$', '', .col)}.3c")
	)

for (Y in names(dx.lst)) {
	dat[grep(paste0("^", Y, "\\.([Y]?(t2e|r2e))$"), names(dat))] <- NULL
	dat <- t2e(dat, domain_of(Y), paste0("fod_icd10_", Y), "birth_date", "date_attend", "date_lost", "date_death", date_follow_end, Y, "year")
}

dat <- dat %>% mutate(across(any_of(diet.inc.pts), ~ mk_hml(., group_pct_th), .names = "{sub('\\\\.pts$', '', .col)}.hml"))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 表1
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("表1")
score.cols <- c("bmi.pts","bp.pts","nonhdl.pts","smoke.pts","pa.pts","sleep.pts",diet.s100.cols)
lbls <- list(
	age ~ "Age", female ~ "Female", tdi ~ "TDI", bmi ~ "BMI", energy_kcal ~ "Total Energy Intake, kcal/d", bmi.pts ~ "BMI score", bp.pts ~ "Blood pressure score", nonhdl.pts ~ "Non-HDL cholesterol score", smoke.pts ~ "Smoking score", pa.pts ~ "Physical activity score", sleep.pts ~ "Sleep score"
)
for (n in names(diet.lst)) lbls[[paste0("diet.", n, ".s100")]] <- diet.lst[[n]]
reset_gtsummary_theme(); theme_gtsummary_compact()
table1 <- dat %>% mutate(
	MAHA_Group = factor(diet.maha.3c, c("low","middle","high"), c("Low Adherence","Moderate Adherence","High Adherence")), 
	energy_kcal = energy.kJ * 0.239, female = if_else(tolower(as.character(sex.f)) == "female", 1L, 0L, missing = NA_integer_)
	) %>% select(MAHA_Group, age, female, tdi, bmi, energy_kcal, all_of(score.cols)) %>% drop_na(MAHA_Group) %>% 
	tbl_summary(by = MAHA_Group, type = list(female ~ "dichotomous", all_of(c("age","tdi","bmi","energy_kcal",score.cols)) ~ "continuous"), value = list(female ~ 1), statistic = list(all_continuous() ~ "{mean} ({sd})", female ~ "{n} ({p}%)"), digits = list(all_continuous() ~ 1, female ~ c(0,1)), label = lbls, missing = "no") %>% 
	add_overall() %>% modify_header(label = "**Baseline Characteristics**") %>% modify_footnote(all_stat_cols() ~ "Continuous variables are presented as mean (SD)") %>% bold_labels()
table1$table_styling$footnote <- NULL
table1 %>% as_flex_table() %>% flextable::set_table_properties(layout = "autofit") %>% flextable::save_as_docx(path = "Table1.docx")
save_xlsx1(table1$table_body %>% as.data.frame(), "Table1.out.xlsx"); print(table1)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图1. 🥕饮食与疾病的PheWAS富集分析⛪
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图1：PheWAS")
pacman::p_load(PheWAS, writexl)

dat1 <- dat
Xs <- diet.inc.s100
phewas.res <- plot_phewas(dat1, phecode = NA, Xs = Xs, varX = vars.basic)
phewas.res$plots <- Map(\(p, x)
	p + labs(title = var2lab(x)) + theme(
		plot.title = element_text(hjust = 0.5, face = "bold"),
		axis.title = element_text(face = "bold"),
		axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = "bold")
	), phewas.res$plots, Xs
)
names(phewas.res$plots) <- Xs; phewas.res$plots;

lab <- setNames(unname(diet.lst[gsub("^diet\\.|\\.(sum|pts|q5|s100|s100_q5)$", "", Xs)]), Xs)

res <- phewas.res$res %>%
	as_tibble() %>%
	mutate(
		Diet = recode(snp, !!!lab),
		category = tolower(category),
		description = stringr::str_squish(description),
		logp = -log10(p),
		FDR = p.adjust(p, "BH"),
		sig_bonf = bonferroni %in% TRUE
	) %>%
	filter(!is.na(Diet), is.finite(p), p > 0)

save_xlsx1(res %>% filter(logp >= 10), "Fig1.PheWAS.top.xlsx")

sum_score <- res %>% group_by(Diet) %>% summarise(
	n_bonf = sum(sig_bonf), n_fdr = sum(FDR < 0.05), mean_logp = mean(logp), median_logp = median(logp),
	n_neg_bonf = sum(sig_bonf & beta < 0), n_pos_bonf = sum(sig_bonf & beta > 0), .groups = "drop"
) %>% arrange(desc(n_bonf))

sum_cat <- res %>% filter(Diet %in% c("DASH", "MAHA")) %>%
	group_by(Diet, category) %>%
	summarise(n_bonf = sum(sig_bonf), prop_bonf = mean(sig_bonf), mean_logp = mean(logp), .groups = "drop")

dm <- res %>% filter(Diet %in% c("DASH", "MAHA")) %>%
	select(Diet, phenotype, description, category, beta, p, logp, sig_bonf) %>%
	pivot_wider(names_from = Diet, values_from = c(beta, p, logp, sig_bonf), names_sep = ".") %>%
	mutate(
		same_direction = sign(beta.DASH) == sign(beta.MAHA),
		delta_logp = logp.DASH - logp.MAHA,
		delta_beta = beta.DASH - beta.MAHA,
		sig_pattern = case_when(
			sig_bonf.DASH & sig_bonf.MAHA ~ "Both",
			sig_bonf.DASH & !sig_bonf.MAHA ~ "DASH_only",
			!sig_bonf.DASH & sig_bonf.MAHA ~ "MAHA_only", TRUE ~ "Neither"
		)
	)

dash_only <- dm %>% filter(sig_pattern == "DASH_only") %>% arrange(desc(logp.DASH))
maha_only <- dm %>% filter(sig_pattern == "MAHA_only") %>% arrange(desc(logp.MAHA))

cat_diff <- dm %>%
	group_by(category) %>%
	summarise(
		n_both = sum(sig_pattern == "Both"),
		n_dash_only = sum(sig_pattern == "DASH_only"),
		n_maha_only = sum(sig_pattern == "MAHA_only"),
		mean_delta_logp = mean(delta_logp),
		prop_same_direction = mean(same_direction),
		.groups = "drop"
	) %>%
	arrange(desc(mean_delta_logp))

t2dm_tbl <- res %>% filter(description == "Type 2 diabetes") %>% arrange(p) %>% select(Diet, beta, SE, p, logp, sig_bonf, n_total)

p_cat <- ggplot(sum_cat, aes(x = forcats::fct_reorder(category, n_bonf), y = n_bonf, fill = Diet)) +
	geom_col(position = "dodge") +
	coord_flip(clip = "off") +
	scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
	labs(x = NULL, y = "Bonferroni-significant phenotypes", fill = NULL) +
	theme_bw(base_size = 14) +
	theme(
		legend.position = c(0.72, 0.52), legend.direction = "vertical", legend.justification = c(0, 0.5),
		legend.background = element_rect(fill = scales::alpha("white", 0.85), colour = "grey80"),
		legend.key.size = unit(0.5, "cm"),
		axis.title = element_text(face = "bold"),
		axis.text = element_text(face = "bold")
	)

lab_dat <- dm %>% filter(pmin(p.DASH, p.MAHA) <= 1e-10) %>% slice_max(abs(delta_logp), n = 12) %>% mutate(dx = ifelse(logp.DASH >= logp.MAHA, 0.7, -0.7), dy = ifelse(logp.DASH >= logp.MAHA, -0.5, 0.5), hjust = ifelse(dx > 0, 0, 1))
p_cmp_p <- ggplot(dm, aes(logp.DASH, logp.MAHA)) +
	geom_point(alpha = 0.45) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
	geom_text(data = lab_dat, aes(x = logp.DASH + dx, y = logp.MAHA + dy, label = description, hjust = hjust), size = 3, fontface = "bold", check_overlap = TRUE) +
	labs(x = "DASH: -log10(P)", y = "MAHA: -log10(P)") +
	theme_bw(base_size = 14) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))

lab_dat_b <- dm %>% filter(pmin(p.DASH, p.MAHA) <= 1e-10) %>% slice_max(abs(delta_beta), n = 12) %>% mutate(dx = ifelse(beta.DASH >= beta.MAHA, 0.0006, -0.0006), dy = ifelse(beta.DASH >= beta.MAHA, -0.0006, 0.0006), hjust = ifelse(dx > 0, 0, 1))
p_cmp_beta <- ggplot(dm, aes(beta.DASH, beta.MAHA)) +
	geom_point(alpha = 0.45) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
	geom_text(data = lab_dat_b, aes(x = beta.DASH + dx, y = beta.MAHA + dy, label = description, hjust = hjust), size = 3, fontface = "bold", check_overlap = TRUE) +
	labs(x = "DASH: beta", y = "MAHA: beta") +
	theme_bw(base_size = 14) +
	theme(axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"))


Fig1 <- p_cat | (phewas.res$plots[["diet.maha.s100"]] / phewas.res$plots[["diet.dash.s100"]])
save_plot(Fig1, "Fig1.png", width = 14, height = 10, dpi = 320)
save_xlsx1(
	list(
		phewas_full = res,
		score_summary = sum_score,
		category_summary = sum_cat,
		dash_only = dash_only,
		maha_only = maha_only,
		category_difference = cat_diff,
		t2dm = t2dm_tbl
	),
	"Fig1.out.xlsx"
)
print(sum_score)
print(sum_cat, n = 50)
print(head(dash_only, 50))
print(maha_only, n = 20)
print(cat_diff, n = 30)
print(t2dm_tbl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图2. diet与19种疾病的关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图2：diet与19种疾病的关联")
dat1 <- dat
for (x in diet.inc.pts) dat1[[x]] <- as.numeric(scale(dat1[[x]]))

assoc.res <- assoc_reg(dat1, diet.inc.pts, covs, names(dx.lst), type = "t2e") %>%
	mutate(lab.Y = coalesce(dx.lst[Outcome], Outcome), Exposure0 = gsub("^diet\\.|\\.(pts|q5|qt)$", "", Exposure), Diet = diet.lst[Exposure0])

write_xlsx(
	assoc.res %>% transmute(
		Outcome = lab.Y, Exposure = Diet, HR = sprintf("%.3f", estimate),
		CI = sprintf("%.3f–%.3f", conf.low, conf.high), P = format(p.value, scientific = TRUE, digits = 3),
		HR_raw = estimate, CI_low = conf.low, CI_high = conf.high, P_raw = p.value, N_total, N_event
	), "Fig2.data.xlsx"
)

Y.inc1 <- Y.inc[1:3]; Y.col1 <- Y.cols[1:3]
Y.inc2 <- Y.inc[4:6]; Y.col2 <- Y.cols[4:6]

Fig2 <- ((plot_forest2(assoc.res, Y.inc1, Y.col1, xlim = c(0.84, 1.02), show_legend = TRUE) |
	plot_forest2(assoc.res, Y.inc2, Y.col2, xlim = c(0.78, 1.08), show_legend = TRUE)) +
	plot_layout(guides = "collect")) & theme(legend.position = "bottom")

save_plot(Fig2, "Fig2.png", width = 8, height = 8)
save_xlsx1(assoc.res %>% as.data.frame(), "Fig2.out.xlsx")
print(assoc.res)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图3. 增量模型比较
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图3：增量模型比较")
diet.pairs <- c(dash = "DASH", medi24 = "MEDI", mind = "MIND")
cols.pattern <- c(DASH = "#16B9C0", MEDI = "#7CAE00", MIND = "#B97AF7")
cols.info <- c("Traditional-only information" = "#16B9C0", "MAHA-only information" = "#F26D60")
outcome_levels <- rev(dx.lst[Y.inc]); pattern_levels <- unname(diet.pairs)

fit_disc_one <- function(dat, Y, trad_nm, trad_lab, covs) {
	trad_var <- paste0("diet.", trad_nm, ".hml")
	d0 <- dat %>%
		transmute(time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]],
			across(all_of(covs)), trad = .data[[trad_var]], maha = diet.maha.hml) %>%
		filter(trad %in% c("low", "high"), maha %in% c("low", "high")) %>% drop_na()
	g_hi <- paste0(trad_lab, " high + MAHA high"); g_hl <- paste0(trad_lab, " high + MAHA low")
	g_lh <- paste0(trad_lab, " low + MAHA high");  g_ll <- paste0(trad_lab, " low + MAHA low")
	d0 <- d0 %>% mutate(
		group = factor(case_when(trad == "high" & maha == "high" ~ g_hi, trad == "high" & maha == "low" ~ g_hl,
			trad == "low" & maha == "high" ~ g_lh, TRUE ~ g_ll), levels = c(g_hi, g_hl, g_lh, g_ll)),
		discord = factor(case_when(trad == "high" & maha == "low" ~ g_hl, trad == "low" & maha == "high" ~ g_lh,
			TRUE ~ NA_character_), levels = c(g_hl, g_lh))
	)
	res4 <- coxph(Surv(time, event) ~ group + ., data = d0 %>% select(time, event, all_of(covs), group)) %>%
		broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("^group", term)) %>%
		transmute(Pattern = trad_lab, Y, group = sub("^group", "", term), estimate, conf.low, conf.high, p.value)
	res4 <- bind_rows(tibble(Pattern = trad_lab, Y, group = g_hi, estimate = 1, conf.low = 1, conf.high = 1, p.value = NA_real_), res4)
	res2 <- coxph(Surv(time, event) ~ discord + ., data = d0 %>% filter(!is.na(discord)) %>% select(time, event, all_of(covs), discord)) %>%
		broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>% filter(grepl("^discord", term)) %>%
		transmute(Pattern = trad_lab, Y, contrast = paste0(g_lh, " vs ", g_hl), estimate, conf.low, conf.high, p.value)
	list(res4 = res4, res2 = res2)
}

fit_aic_one <- function(dat, Y, trad_nm, trad_lab, covs) {
	d0 <- dat %>% transmute(time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]],
		across(all_of(covs)), maha = as.numeric(scale(diet.maha.pts)),
		trad = as.numeric(scale(.data[[paste0("diet.", trad_nm, ".pts")]]))) %>% drop_na()
	get_aic <- function(v) AIC(coxph(Surv(time, event) ~ ., data = d0 %>% select(time, event, all_of(covs), all_of(v))))
	a0 <- get_aic(character(0)); am <- get_aic("maha"); at <- get_aic("trad"); ab <- get_aic(c("maha", "trad"))
	tibble(Pattern = trad_lab, Y, Base = a0, MAHA = am, Traditional = at, Both = ab,
		add_MAHA_given_trad = at - ab, add_trad_given_MAHA = am - ab,
		net_traditional_gain = (am - ab) - (at - ab))
}

res.lst <- purrr::imap(diet.pairs, ~ purrr::map(Y.inc, fit_disc_one, dat = dat, trad_nm = .y, trad_lab = .x, covs = covs))
tab.all <- purrr::map_dfr(res.lst, ~ purrr::map_dfr(.x, "res4")) %>% mutate(Outcome = dx.lst[Y]) %>%
	select(Pattern, Outcome, group, estimate, conf.low, conf.high, p.value)
tab.discord <- purrr::map_dfr(res.lst, ~ purrr::map_dfr(.x, "res2")) %>%
	mutate(Outcome = factor(dx.lst[Y], levels = outcome_levels), Pattern = factor(Pattern, levels = pattern_levels)) %>%
	select(Pattern, Outcome, contrast, estimate, conf.low, conf.high, p.value)
tab.ext <- purrr::imap_dfr(diet.pairs, ~ purrr::map_dfr(Y.inc, fit_aic_one, dat = dat, trad_nm = .y, trad_lab = .x, covs = covs)) %>%
	mutate(Outcome = dx.lst[Y], winner = case_when(abs(add_trad_given_MAHA - add_MAHA_given_trad) < 2 ~ "Similar",
		add_trad_given_MAHA > add_MAHA_given_trad ~ "Traditional stronger", TRUE ~ "MAHA stronger"))
tab.ext2 <- tab.ext %>% mutate(
	Outcome = factor(Outcome, levels = outcome_levels), Pattern = factor(Pattern, levels = pattern_levels),
	maha_pos = pmax(add_MAHA_given_trad, 0), trad_pos = pmax(add_trad_given_MAHA, 0), total_unique = maha_pos + trad_pos,
	pct_trad = ifelse(total_unique > 0, 100 * trad_pos / total_unique, 50), net_trad = trad_pos - maha_pos
)
tab.sum <- tab.ext2 %>% group_by(Pattern) %>% summarise(
	mean_total_unique = mean(total_unique, na.rm = TRUE), mean_pct_trad = mean(pct_trad, na.rm = TRUE),
	n_trad_better = sum(net_trad > 2), n_MAHA_better = sum(net_trad < -2),
	n_similar = sum(abs(net_trad) <= 2 & total_unique >= 5), n_minimal = sum(total_unique < 5), .groups = "drop"
)

save_xlsx1(list(`4group_results` = tab.all, discordant_contrast = tab.discord, incremental_AIC = tab.ext,
	incremental_dominance = tab.ext2, incremental_summary = tab.sum), "Fig3.out.xlsx")

shift_y <- c(DASH = -0.22, MEDI = 0, MIND = 0.22)
tab.a <- tab.discord %>% mutate(y0 = as.numeric(Outcome), y = y0 + unname(shift_y[as.character(Pattern)]))
xlim_A <- range(c(1, tab.a$conf.low, tab.a$conf.high), na.rm = TRUE); xpad_A <- 0.06 * diff(xlim_A)
pA <- ggplot(tab.a, aes(color = Pattern)) +
	geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.8) +
	geom_segment(aes(x = conf.low, xend = conf.high, y = y, yend = y), linewidth = 1.0) +
	geom_point(aes(x = estimate, y = y), size = 3.6) +
	scale_color_manual(values = cols.pattern) +
	scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.45, 0.45))) +
	coord_cartesian(xlim = c(xlim_A[1] - xpad_A, xlim_A[2] + xpad_A), clip = "off") +
	labs(x = "Hazard ratio", y = NULL, color = NULL) + theme_classic(base_size = 16) +
	theme(legend.position = "bottom", legend.text = element_text(size = 12, face = "bold"),
		axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"))

max_total <- max(tab.ext2$total_unique, na.rm = TRUE); if (!is.finite(max_total) || max_total <= 0) max_total <- 1
tab.b0 <- tab.ext2 %>% transmute(Pattern, Outcome, y = as.numeric(Outcome), trad = ifelse(round(trad_pos) == 0, 0, trad_pos), maha = ifelse(round(maha_pos) == 0, 0, maha_pos))
tab.b.bg <- expand_grid(Pattern = factor(pattern_levels, levels = pattern_levels), Outcome = factor(outcome_levels, levels = outcome_levels)) %>%
	mutate(y = as.numeric(Outcome), xmin = 0, xmax = max_total)
tab.b.rect <- bind_rows(
	tab.b0 %>% transmute(Pattern, Outcome, y, Source = "Traditional-only information", xmin = 0, xmax = trad, label = ifelse(round(trad) >= 10, as.character(round(trad)), "")),
	tab.b0 %>% transmute(Pattern, Outcome, y, Source = "MAHA-only information", xmin = trad, xmax = trad + maha, label = ifelse(round(maha) >= 10, as.character(round(maha)), ""))
) %>% filter(xmax > xmin)
tab.b.txt <- tab.b.rect %>% mutate(x = (xmin + xmax) / 2) %>% filter(label != "")
pB <- ggplot() +
	geom_rect(data = tab.b.bg, aes(xmin = xmin, xmax = xmax, ymin = y - 0.34, ymax = y + 0.34), fill = "grey94", color = NA) +
	geom_rect(data = tab.b.rect, aes(xmin = xmin, xmax = xmax, ymin = y - 0.34, ymax = y + 0.34, fill = Source), color = NA) +
	geom_text(data = tab.b.txt, aes(x = x, y = y, label = label), size = 3.4, fontface = "bold") +
	facet_grid(. ~ Pattern, switch = "x") + scale_fill_manual(values = cols.info, breaks = names(cols.info)) +
	scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
	scale_y_continuous(breaks = seq_along(outcome_levels), labels = outcome_levels, expand = expansion(add = c(0.5, 0.5))) +
	labs(x = "Unique information (ΔAIC units)", y = NULL, fill = NULL) + theme_classic(base_size = 15) +
	theme(legend.position = "top", legend.text = element_text(size = 11, face = "bold"), strip.placement = "outside",
		strip.background = element_blank(), strip.text.x = element_text(face = "bold", size = 13),
		axis.text.y = element_text(face = "bold"), axis.title.x = element_text(face = "bold"),
		panel.spacing.x = unit(0.8, "lines"), plot.margin = margin(2, 5.5, 5.5, 5.5))

Fig3 <- cowplot::plot_grid(pA, pB, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1.05))
save_plot(Fig3, "Fig3.png", width = 11.6, height = 10.2, dpi = 320)
print(tab.all, n = 30)
print(tab.discord, n = 30)
print(tab.ext, n = 30)
print(tab.sum)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图4. DASH 🤺 MAHA 📊🧱
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图4：DASH vs MAHA 的10年风险")
dat1 <- dat %>% mutate(diet.dash.3c = factor(diet.dash.3c, levels = c("low", "middle", "high")), diet.maha.3c = factor(diet.maha.3c, levels = c("low", "middle", "high"))) 

Y <- "t2dm"
plot_cuminc_joint <- function(dat1, Y_time, Y_event, Xfacet, Xline, covs, t0 = 10, Xfacetlab = NULL, Xlinelab = NULL, title = NULL, leg_position = "right", leg_direction = "vertical") {
	d <- dat1[!is.na(dat1[[Xfacet]]) & !is.na(dat1[[Xline]]), ]
	lvf <- levels(factor(d[[Xfacet]])); lvl <- levels(factor(d[[Xline]]))
	d$xf <- factor(d[[Xfacet]], levels = lvf); d$xl <- factor(d[[Xline]], levels = lvl)
	covs2 <- covs[covs %in% names(d)]
	d <- d[complete.cases(d[, c("xf", "xl", Y_time, Y_event, covs2), drop = FALSE]), ]
	form <- as.formula(paste0("survival::Surv(", Y_time, ", ", Y_event, ") ~ xf + xl", if (length(covs2)) paste0(" + ", paste(covs2, collapse = " + ")) else ""))
	fit <- survival::coxph(form, data = d)
	nd <- expand.grid(xf = factor(lvf, levels = lvf), xl = factor(lvl, levels = lvl))
	if (length(covs2)) for (v in covs2) nd[[v]] <- if (is.numeric(d[[v]])) mean(d[[v]], na.rm = TRUE) else names(which.max(table(d[[v]])))
	tt <- seq(0, t0, by = 0.1)
	res <- bind_rows(lapply(seq_len(nrow(nd)), function(i) { s <- summary(survival::survfit(fit, newdata = nd[i, , drop = FALSE]), times = tt, extend = TRUE); tibble(xf = nd$xf[i], xl = nd$xl[i], time = s$time, y = 1 - s$surv, cl = 1 - s$upper, cu = 1 - s$lower) }))
	ggplot(res, aes(time, y, color = xl, fill = xl)) +
		geom_ribbon(aes(ymin = cl, ymax = cu), alpha = 0.15, linewidth = 0) +
		geom_line(linewidth = 1.1) +
		facet_wrap(~xf, nrow = 1) +
		scale_x_continuous(limits = c(0, t0), breaks = c(0, 5, 10), labels = c("0", "5", "10")) +
		scale_y_continuous(labels = scales::label_percent(accuracy = 1), expand = expansion(mult = c(0, 0.08))) +
		scale_color_brewer(palette = "Set2") + scale_fill_brewer(palette = "Set2") + theme_minimal() +
		labs(title = title, x = "Follow-up time (years)", y = "Cumulative incidence", color = if (is.null(Xlinelab)) Xline else Xlinelab, fill = if (is.null(Xlinelab)) Xline else Xlinelab) +
		theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(face = "bold"), axis.text = element_text(size = 11, face = "bold"), legend.text = element_text(size = 10, face = "bold"), legend.title = element_text(size = 11, face = "bold"), strip.text = element_text(size = 10, face = "bold"), legend.position = leg_position, legend.direction = leg_direction)
}

pa <- plot_risk(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), X1 = "diet.dash.3c", X2 = "diet.maha.3c", covs = vars.basic, method = "10years", t0 = 10, X1lab = "DASH", X2lab = "MAHA", group_bgcolor = TRUE, title = NULL, leg_position = "top", leg_direction = "horizontal", tab = TRUE)  
pb <- plot_risk(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), X1 = "diet.maha.3c", X2 = "diet.dash.3c", covs = vars.basic, method = "10years", t0 = 10, X1lab = "MAHA", X2lab = "DASH", group_bgcolor = TRUE, title = NULL, leg_position = "top", leg_direction = "horizontal", tab = TRUE)  

pc <- plot_cuminc_joint(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), Xfacet = "diet.maha.3c", Xline = "diet.dash.3c", covs = vars.basic, t0 = 10, Xfacetlab = "MAHA", Xlinelab = "DASH", title = NULL, leg_position = "right")  
pd <- plot_cuminc_joint(dat1, Y_time = paste0(Y, ".t2e"), Y_event = paste0(Y, ".Yt2e"), Xfacet = "diet.dash.3c", Xline = "diet.maha.3c", covs = vars.basic, t0 = 10, Xfacetlab = "DASH", Xlinelab = "MAHA", title = NULL, leg_position = "right")  

Fig4 <- ((pa | pb) / (pc | pd)) +
	plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold")))
save_plot(Fig4, "Fig4.png", width = 10.8, height = 11.0, dpi = 320)
fig4_tab <- dat1 %>% count(diet.dash.3c, diet.maha.3c, name = "N") %>% mutate(Outcome = dx.lst[Y])
save_xlsx1(fig4_tab, "Fig4.out.xlsx")
print(fig4_tab)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图5. 🗺 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图5：地理分布")
# remotes::install_github("ropensci/rnaturalearthhires")
pacman::p_load(dplyr, sf, ggplot2, viridis, rnaturalearth, patchwork, rlang)

uk_map <- ne_countries(scale = "medium", country = "United Kingdom", returnclass = "sf") %>% st_transform(27700)

plot_center_diet <- function(dat, x_col, y_col, score_col, center.inc = NULL, center.topN = NA, ncol = 1, title = TRUE, legend = TRUE, col_title = NULL) {
	x_str <- as_name(enquo(x_col)); y_str <- as_name(enquo(y_col))
	if (!is.null(center.inc)) {
		m <- match(tolower(center.inc), tolower(center.lst))
		target_codes <- names(center.lst)[m[!is.na(m)]]
	} else {
		target_codes <- dat %>% filter(!is.na({{x_col}}), !is.na({{y_col}}), !is.na({{score_col}})) %>% count(center) %>% slice_max(n, n = center.topN) %>% pull(center) %>% as.character()
	}
	plot_list <- lapply(target_codes, function(code) {
		dat_c <- dat %>%
			filter(as.character(center) == code, !is.na({{x_col}}), !is.na({{y_col}}), !is.na({{score_col}})) %>%
			mutate(m_e = median({{x_col}}), m_n = median({{y_col}}), d = sqrt(({{x_col}} - m_e)^2 + ({{y_col}} - m_n)^2)) %>%
			filter(d <= quantile(d, 0.98))
		agg <- dat_c %>% group_by({{x_col}}, {{y_col}}) %>% summarise(score = mean({{score_col}}), .groups = "drop")
		box <- st_bbox(st_as_sf(agg, coords = c(x_str, y_str), crs = 27700))
		p <- ggplot() +
			geom_sf(data = uk_map, fill = "grey92", color = "white", linewidth = 0.3) +
			geom_tile(data = agg, aes(x = {{x_col}}, y = {{y_col}}, fill = score), width = 1000, height = 1000) +
			scale_fill_viridis_c(option = "turbo", name = NULL) +
			coord_sf(datum = NA, expand = FALSE, xlim = c(box["xmin"] - 1000, box["xmax"] + 1000), ylim = c(box["ymin"] - 1000, box["ymax"] + 1000)) +
			theme_void() +
			theme(plot.margin = margin(2, 0, 2, 0))
		if (title) p <- p + ylab(paste0(ifelse(is.na(center.lst[code]), code, center.lst[code]), " (N=", format(nrow(dat_c), big.mark = ","), ")")) + theme(axis.title.y = element_text(angle = 270, face = "bold", size = 13, margin = margin(r = 5)))
		if (legend) p <- p + theme(legend.position = "right", legend.key.height = unit(0.9, "cm"), legend.key.width = unit(0.45, "cm"), legend.text = element_text(size = 11)) else p <- p + theme(legend.position = "none")
		p
	})
	if (!is.null(col_title)) plot_list[[1]] <- plot_list[[1]] + ggtitle(col_title) + theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, margin = margin(b = 5)))
	wrap_plots(plot_list, ncol = ncol)
}

FigA <- plot_center_diet(dat, home_east, home_north, diet.maha.s100, center.inc = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), ncol = 1, title = TRUE, legend = FALSE, col_title = "MAHA")
FigB <- plot_center_diet(dat, home_east, home_north, diet.dash.s100, center.inc = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), ncol = 1, title = FALSE, legend = FALSE, col_title = "DASH")
FigC <- plot_center_diet(dat, home_east, home_north, diet.mind.s100, center.inc = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), ncol = 1, title = FALSE, legend = FALSE, col_title = "MIND")
FigD <- plot_center_diet(dat, home_east, home_north, diet.medi24.s100, center.inc = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), ncol = 1, title = FALSE, legend = FALSE, col_title = "MEDI-WebQ")
FigE <- plot_center_diet(dat, home_east, home_north, diet.medito.s100, center.inc = c("Croydon", "Hounslow", "Bristol", "Edinburgh"), ncol = 1, title = FALSE, legend = TRUE, col_title = "MEDI-Touch")
Fig5 <- FigA | FigB | FigC | FigD | FigE; print(Fig5)
save_plot(Fig5, "Fig5.png", width = 14, height = 14, dpi = 320)

center.lst.inc <- c("c11020" = "Croydon", "c11018" = "Hounslow", "c11011" = "Bristol", "c11005" = "Edinburgh")
rename_mapping <- setNames(paste0("diet.", diet.inc, ".s100"), diet.inc)

dat_stat <- dat %>%
	filter(as.character(center) %in% names(center.lst.inc)) %>%
	rename(!!!rename_mapping) %>%
	filter(if_all(all_of(diet.inc), ~ !is.na(.))) %>%
	mutate(center_name = factor(center.lst.inc[as.character(center)], levels = center.lst.inc))

fig5_mean <- dat_stat %>% group_by(center_name) %>% summarise(N = n(), across(all_of(diet.inc), ~ sprintf("%.1f ± %.1f", mean(.x), sd(.x)), .names = "{.col}_Score"), .groups = "drop") %>% as.data.frame()
fig5_cmp <- dat_stat %>% group_by(center_name) %>% summarise(across(all_of(diet.inc[-1]), list(Diff = ~ sprintf("%+.2f", mean(.data[[diet.inc[1]]] - .x)), Cor = ~ sprintf("%.3f", cor(.data[[diet.inc[1]]], .x))), .names = "{.col}_vs_ref_{.fn}"), .groups = "drop") %>% as.data.frame()

save_xlsx1(list(center_score_summary = fig5_mean, center_comparison = fig5_cmp), "Fig5.out.xlsx")
print(fig5_mean)
print(fig5_cmp)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S1. 不同diet评分之间的相关性🍉
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图S1：不同diet评分之间的相关性")
pacman::p_load(pheatmap, irr, ggdist)
dat1 <- dat

s100_vars <- diet.s100.cols[diet.s100.cols %in% names(dat1)]
q5_vars <- diet.q5.cols[diet.q5.cols %in% names(dat1)]

plot_dat <- dat1 %>%
	select(all_of(s100_vars)) %>%
	pivot_longer(everything(), names_to = "var", values_to = "score") %>%
	filter(is.finite(score)) %>%
	mutate(Diet_score = factor(var2lab(var), levels = rev(diet.order)))

plot_sum <- plot_dat %>%
	group_by(Diet_score) %>%
	summarise(
		N = n(), Mean = mean(score), SD = sd(score), Median = median(score),
		P25 = quantile(score, 0.25) %>% as.numeric(), P75 = quantile(score, 0.75) %>% as.numeric(),
		Min = min(score), Max = max(score), .groups = "drop"
	) %>%
	mutate(across(c(Mean, SD, Median, P25, P75, Min, Max), ~ round(.x, 1))) %>%
	as.data.frame()

set.seed(1)
plot_dat_jit <- split(plot_dat, plot_dat$Diet_score) %>% lapply(\(d) d[sample(nrow(d), min(5000, nrow(d))), ]) %>% bind_rows()

pacman::p_load(pheatmap, irr, ggdist, cowplot, grid)

p_top <- ggplot(plot_dat, aes(x = score, y = Diet_score, fill = Diet_score)) +
	ggdist::stat_halfeye(adjust = 0.7, width = 0.6, .width = 0, justification = -0.2, point_colour = NA, alpha = 0.7) +
	geom_jitter(data = plot_dat_jit, height = 0.10, width = 0, alpha = 0.03, size = 0.25, shape = 16, color = "black") +
	geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.85, median.linewidth = 1.0) + # ❗comment: fatten 已弃用，改为 median.linewidth
	scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0.01, 0)) +
	labs(x = "Harmonized dietary score (0–100)", y = NULL, title = "Harmonized dietary pattern scores") +
	theme_classic(base_size = 13) +
	theme(legend.position = "none", plot.title = element_text(face = "bold"), axis.text.y = element_text(face = "bold"))

cor_mat <- cor(dat1 %>% select(all_of(s100_vars)) %>% mutate(across(everything(), as.numeric)), use = "pairwise.complete.obs", method = "spearman")
dimnames(cor_mat) <- rep(list(var2lab(colnames(cor_mat))), 2)

p_cor <- pheatmap(
	cor_mat,
	main = "Spearman correlation",
	display_numbers = matrix(sprintf("%.2f", cor_mat), nrow(cor_mat)),
	number_color = "black", number_fontface = "bold",
	fontsize_number = 12, fontsize_row = 11, fontsize_col = 11,
	angle_col = 90, cellwidth = 32, cellheight = 32,
	clustering_method = "complete", border_color = NA,
	silent = TRUE
)

get_wkappa <- function(x, y) {
	d0 <- data.frame(
		x = factor(x, levels = c(0, 25, 50, 75, 100), ordered = TRUE),
		y = factor(y, levels = c(0, 25, 50, 75, 100), ordered = TRUE)
	) %>% drop_na()
	if (nrow(d0) < 50) return(NA_real_)
	irr::kappa2(d0, weight = "squared")$value
}

agree_stat <- function(x, y) {
	d0 <- data.frame(x = as.numeric(x), y = as.numeric(y)) %>% drop_na()
	if (nrow(d0) < 50) return(c(n = nrow(d0), exact = NA, within1 = NA))
	c(n = nrow(d0), exact = mean(d0$x == d0$y), within1 = mean(abs(d0$x - d0$y) <= 25))
}

wkappa_mat <- outer(q5_vars, q5_vars, Vectorize(\(a, b) if (a == b) 1 else get_wkappa(dat1[[a]], dat1[[b]])))
dimnames(wkappa_mat) <- rep(list(var2lab(q5_vars)), 2)

p_kappa <- pheatmap(
	wkappa_mat,
	main = "Weighted kappa",
	display_numbers = matrix(sprintf("%.2f", wkappa_mat), nrow(wkappa_mat)),
	number_color = "black", number_fontface = "bold",
	fontsize_number = 12, fontsize_row = 11, fontsize_col = 11,
	angle_col = 90, cellwidth = 32, cellheight = 32,
	clustering_method = "complete", border_color = NA,
	silent = TRUE
)

FigS1 <- cowplot::plot_grid(
	cowplot::ggdraw() + cowplot::draw_grob(ggplotGrob(
		p_top + theme(plot.margin = margin(5.5, 5.5, 20, 5.5))
	)),
	ggplot() + theme_void(),
	cowplot::plot_grid(
		cowplot::ggdraw() + cowplot::draw_grob(p_cor$gtable),
		cowplot::ggdraw() + cowplot::draw_grob(p_kappa$gtable),
		ncol = 2, rel_widths = c(1, 1)
	),
	ncol = 1, rel_heights = c(1.12, 0.10, 1)
)
save_plot(FigS1, "FigS1.png", width = 12.5, height = 12.5, dpi = 500, bg = "white")

ref_var_q5 <- "diet.maha.q5" # 📍
agree_tbl <- do.call(rbind, lapply(setdiff(q5_vars, ref_var_q5), \(v) {
	ag <- agree_stat(dat1[[ref_var_q5]], dat1[[v]])
	data.frame(
		var = v,
		label = var2lab(v),
		n = ag["n"],
		spearman_s100 = suppressWarnings(cor(dat1[[sub("\\.q5$", ".s100", ref_var_q5)]], dat1[[sub("\\.q5$", ".s100", v)]], use = "complete.obs", method = "spearman")),
		wkappa = get_wkappa(dat1[[ref_var_q5]], dat1[[v]]),
		exact_agreement = ag["exact"],
		within_1_level = ag["within1"]
	)
})) %>% mutate(label = factor(label, levels = diet.order[diet.order != var2lab(ref_var_q5)])) %>% arrange(label)
save_xlsx1(
	list(
		score_distribution = plot_sum,
		spearman_matrix = as.data.frame(cor_mat, check.names = FALSE),
		weighted_kappa_matrix = as.data.frame(wkappa_mat, check.names = FALSE),
		agreement_vs_ref = agree_tbl
	),
	"FigS1.out.xlsx"
)
print(agree_tbl)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S2. DASH vs MAHA 的 PheWAS ⛪ 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图S2：DASH vs MAHA 的 PheWAS")
pacman::p_load(tidyverse, ggrepel, patchwork)

dm2 <- dm %>% transmute(phenotype, description, beta.MAHA = as.numeric(beta.MAHA), beta.DASH = as.numeric(beta.DASH), logp.MAHA = as.numeric(logp.MAHA), logp.DASH = as.numeric(logp.DASH), delta_logp = as.numeric(delta_logp), delta_beta = as.numeric(delta_beta), sig_pattern = recode(as.character(sig_pattern), DASH_only = "DASH only", MAHA_only = "MAHA only")) %>% filter(if_all(c(beta.MAHA, beta.DASH, logp.MAHA, logp.DASH), is.finite)) %>% mutate(sig_pattern = ifelse(sig_pattern %in% c("Both", "DASH only", "MAHA only"), sig_pattern, "Neither"), logp.MAHA.cap = pmin(logp.MAHA, 15), logp.DASH.cap = pmin(logp.DASH, 15))
bonf <- -log10(0.05 / nrow(dm2)); cols <- c("Both" = "#3B82F6", "DASH only" = "#F97316", "MAHA only" = "#10B981")
lab_p <- bind_rows(dm2 %>% filter(sig_pattern == "DASH only") %>% slice_max(delta_logp, n = 10), dm2 %>% filter(sig_pattern == "MAHA only") %>% slice_min(delta_logp, n = 4), dm2 %>% filter(sig_pattern == "Both") %>% slice_max(abs(delta_logp), n = 4)) %>% distinct(phenotype, .keep_all = TRUE)
lab_b <- bind_rows(dm2 %>% filter(sig_pattern == "DASH only") %>% slice_max(abs(delta_beta), n = 10), dm2 %>% filter(sig_pattern == "MAHA only") %>% slice_max(abs(delta_beta), n = 4), dm2 %>% filter(sig_pattern == "Both") %>% slice_max(abs(delta_beta), n = 4)) %>% distinct(phenotype, .keep_all = TRUE)
blim <- max(quantile(abs(c(dm2$beta.DASH, dm2$beta.MAHA)), 0.995, na.rm = TRUE), max(abs(c(lab_b$beta.DASH, lab_b$beta.MAHA)), na.rm = TRUE)) * 1.08
thm <- theme_bw(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(face = "bold"), axis.text = element_text(face = "bold"), legend.position = "top", legend.title = element_blank(), legend.text = element_text(face = "bold", size = 14))

p1 <- ggplot() +
	geom_point(data = dm2 %>% filter(sig_pattern == "Neither"), aes(logp.DASH.cap, logp.MAHA.cap), color = "grey70", alpha = 0.35, size = 1) +
	geom_point(data = dm2 %>% filter(sig_pattern != "Neither"), aes(logp.DASH.cap, logp.MAHA.cap, color = sig_pattern, shape = sig_pattern), alpha = 0.9, size = 2.3) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + geom_vline(xintercept = bonf, linetype = 3, color = "grey60") + geom_hline(yintercept = bonf, linetype = 3, color = "grey60") +
	geom_text_repel(data = lab_p, aes(logp.DASH.cap, logp.MAHA.cap, label = description, color = sig_pattern), size = 3.2, fontface = "bold", box.padding = 0.25, point.padding = 0.15, segment.alpha = 0.5, max.overlaps = Inf, show.legend = FALSE) +
	scale_color_manual(values = cols) + scale_shape_manual(values = c("Both" = 16, "DASH only" = 17, "MAHA only" = 15)) +
	coord_cartesian(xlim = c(0, 15), ylim = c(0, 15), expand = FALSE) + scale_x_continuous(breaks = seq(0, 15, 3)) + scale_y_continuous(breaks = seq(0, 15, 3)) +
	labs(x = "DASH: -log10(P)", y = "MAHA: -log10(P)", title = "Association strength") + thm

p2 <- ggplot() +
	geom_point(data = dm2 %>% filter(sig_pattern == "Neither"), aes(beta.DASH, beta.MAHA), color = "grey70", alpha = 0.35, size = 1) +
	geom_point(data = dm2 %>% filter(sig_pattern != "Neither"), aes(beta.DASH, beta.MAHA, color = sig_pattern, shape = sig_pattern), alpha = 0.9, size = 2.3) +
	geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") + geom_hline(yintercept = 0, linetype = 3, color = "grey60") + geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
	geom_text_repel(data = lab_b, aes(beta.DASH, beta.MAHA, label = description, color = sig_pattern), size = 3.2, fontface = "bold", box.padding = 0.25, point.padding = 0.15, segment.alpha = 0.5, max.overlaps = Inf, show.legend = FALSE) +
	scale_color_manual(values = cols) + scale_shape_manual(values = c("Both" = 16, "DASH only" = 17, "MAHA only" = 15)) +
	coord_cartesian(xlim = c(-blim, blim), ylim = c(-blim, blim), expand = FALSE) +
	labs(x = "DASH: beta", y = "MAHA: beta", title = "Effect size") + thm

FigS2 <- (p1 / p2) + plot_layout(guides = "collect") & theme(legend.position = "top"); FigS2
save_plot(FigS2, "FigS2.png", width = 16, height = 16, dpi = 500, bg = "white")
save_xlsx1(list(phewas_compare = dm2, label_p = lab_p, label_beta = lab_b), "FigS2.out.xlsx")
print(dm2, n = 30)
print(lab_p)
print(lab_b)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S3. 联合打分✝年风险 🔗
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图S3：联合打分10年风险")
run_joint_risk_one <- function(dat, Y, covs, t0 = 10) {
	d2 <- dat %>% transmute(time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(all_of(covs)), dash = factor(diet.dash.3c, levels = c("low", "middle", "high")), maha = factor(diet.maha.3c, levels = c("low", "middle", "high"))) %>% drop_na() 
	fit <- coxph(Surv(time, event) ~ dash * maha + ., data = d2 %>% select(time, event, dash, maha, all_of(covs)))
	nd <- crossing(dash = factor(c("low", "middle", "high"), levels = c("low", "middle", "high")), maha = factor(c("low", "middle", "high"), levels = c("low", "middle", "high"))); for (v in covs) nd[[v]] <- make_typical_value(d2[[v]]); stopifnot(!anyNA(nd))
	bind_rows(lapply(seq_len(nrow(nd)), function(i) { ss <- summary(survfit(fit, newdata = nd[i, , drop = FALSE]), times = t0, extend = TRUE); data.frame(Y = Y, dash = as.character(nd$dash[i]), maha = as.character(nd$maha[i]), risk10 = round(100 * (1 - ss$surv), 2)) }))
}

tab2.risk <- bind_rows(lapply(Y.inc, function(Y) run_joint_risk_one(dat, Y, covs)))
tab2.risk.out <- tab2.risk %>% mutate(Outcome = dx.lst[Y]) %>% select(Outcome, dash, maha, risk10)
save_xlsx1(tab2.risk.out, "FigS3.out.xlsx")

plot_joint_heat <- function(d, ttl) {
	d <- d %>% mutate(dash = factor(dash, levels = c("low", "middle", "high")), maha = factor(maha, levels = c("low", "middle", "high")))
	ggplot(d, aes(x = maha, y = dash, fill = risk10)) +
		geom_tile(color = "white", linewidth = 0.8) +
		geom_text(aes(label = sprintf("%.2f%%", risk10)), size = 4.0) +
		scale_x_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
		scale_y_discrete(labels = c(low = "Low", middle = "Middle", high = "High")) +
		scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
		labs(title = ttl, x = "MAHA adherence", y = "DASH adherence", fill = "10-year risk") +
		theme_classic(base_size = 14) + theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

plist.risk <- lapply(Y.inc, function(Y) plot_joint_heat(tab2.risk %>% filter(Y == !!Y), dx.lst[Y]))
FigS3 <- wrap_plots(plist.risk, ncol = 2)
save_plot(FigS3, "FigS3.png", width = 8, height = 8)
print(tab2.risk.out)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🚩 图S4. 从MAHA和DASH的每个组成部分
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
step_header("图S4：MAHA和DASH的组成部分")
if (!exists("field")) field <- readxl::read_excel(paste0(dir0, "/files/foods.xlsx"), sheet = "Sheet2")

diet.dict <- read.delim(
	ifelse(file.exists(paste0(indir, "/common/diet.lst")), paste0(indir, "/common/diet.lst"), "diet.lst"),
	sep = "\t", header = FALSE, fill = TRUE, quote = "", comment.char = ""
) %>%
	as_tibble() %>%
	transmute(item = trimws(V1), item_label = gsub("_", " ", trimws(V2), fixed = TRUE)) %>%
	distinct(item, .keep_all = TRUE)

collapse1 <- function(dat, v, per = 0.8) { cc <- grep(paste0("^", v, "_i[0-4]$"), names(dat), value = TRUE); if (length(cc) == 0) rep(NA_real_, nrow(dat)) else rowMeans1(dat[, cc, drop = FALSE], per = per) }

fit_disc <- function(v, dat.item) {
	d0 <- dat.item %>% transmute(y_disc, x = zstd(.data[[v]]), across(any_of(covs))) %>% filter(!is.na(y_disc), !is.na(x)) %>% drop_na(any_of(covs))
	if (nrow(d0) < 200 || is.na(sd(d0$x, na.rm = TRUE)) || sd(d0$x, na.rm = TRUE) == 0) return(NULL)
	glm(as.formula(paste0("y_disc ~ x + ", paste(covs, collapse = " + "))), family = binomial(), data = d0) %>%
		broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
		filter(term == "x") %>%
		transmute(item = v, OR = estimate, conf.low, conf.high, p.value, beta = log(OR), n = nrow(d0))
}

fit_cox <- function(v, Y, dat.item0) {
	d0 <- dat %>%
		transmute(eid = as.character(eid), time = .data[[paste0(Y, ".t2e")]], event = .data[[paste0(Y, ".Yt2e")]], across(any_of(covs))) %>%
		left_join(dat.item0 %>% select(eid, all_of(v)), by = "eid") %>%
		transmute(time, event, x = zstd(.data[[v]]), across(any_of(covs))) %>%
		filter(!is.na(time), !is.na(event), !is.na(x)) %>%
		drop_na(any_of(covs))
	if (nrow(d0) < 200 || sum(d0$event == 1, na.rm = TRUE) < 20 || is.na(sd(d0$x, na.rm = TRUE)) || sd(d0$x, na.rm = TRUE) == 0) return(NULL)
	coxph(as.formula(paste0("Surv(time, event) ~ x + ", paste(covs, collapse = " + "))), data = d0, ties = "efron") %>%
		broom::tidy(exponentiate = TRUE, conf.int = TRUE) %>%
		filter(term == "x") %>%
		transmute(item = v, Y, Outcome = dx.lst[Y], HR = estimate, conf.low, conf.high, p.value, beta = log(HR))
}

dairy.vars <- c("p26150", "p26096", "p102820", "p102830", "p102840", "p102860", "p102880", "p102890", "p102900", "p102910")
veg.vars <- c("p104090", "p104240", "p104370", "p104140", "p104160", "p104180", "p104300", "p104310", "p104220", "p104230", "p104260", "p104340", "p104350", "p104150", "p104170", "p104290", "p104330", "p104130", "p104190", "p104270", "p104360", "p104060", "p104070", "p104200", "p104210", "p104250", "p104320", "p104380")
fruit.vars <- c("p104450", "p104560", "p104470", "p104480", "p104490", "p104530", "p104540", "p104430", "p104420", "p104550", "p104410", "p104580", "p104500", "p104460", "p104510", "p104570", "p104520", "p104440", "p104590")
fat.vars <- c("p26110", "oil", "p26062", "p26063", "p102440", "p102450", "p103160", "p104100")
wholegrain.vars <- c("p100770", "p100810", "p100800", "p100840", "p100850", "bread.wholemeal", "baguette.wholemeal", "bap.wholemeal", "roll.wholemeal", "bread.seeded", "baguette.seeded", "bap.seeded", "roll.seeded", "p102720", "p102740", "p102780")
refined.vars <- c("p100820", "p101230", "p101240", "p101250", "p101260", "p101270", "p102710", "p102730", "p102750", "p102760", "p102770", "bread.white", "baguette.white", "bap.white", "roll.white", "p100830", "p100860")
ssb.vars <- c("p100170", "p100180", "p100550", "p100160", "p100540", "p100530", "p100230", "p100190", "p100200", "p100210", "p100220")
sweet.vars <- c("p26064", "p102120", "p102140", "p102150", "p102170", "p102180", "p102190", "p102200", "p102210", "p102220", "p102230", "p102260", "p102270", "p102280", "p102290", "p102300", "p102310", "p102320", "p102330", "p102340", "p102350", "p102360", "p102370", "p102380", "p102060", "p102010", "p102020", "p102050", "p102070", "p101970", "p101980", "p101990", "p102030")
hp.vars <- c("p102460", "p102470", "p102480", "p102500", "p102040", "p102000", "p104020", "p103050", "p103170", "p103180", "p103010", "p103070", "p103080")
alcohol.vars <- c("p26067", "p26138", "p26151", "p26152", "p26153")
meat.vars <- c("p103010", "p103050", "p103070", "p103080")

raw.items <- unique(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars))
raw.items <- raw.items[grepl("^p\\d+$", raw.items)]
derived.items <- unique(c(setdiff(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars), raw.items), "w.p26005", "w.p26052"))

diet.raw <- readRDS(paste0(indir, "/Rdata/diet0.rds")) %>% mutate(eid = as.character(eid))
dat.food <- readRDS(paste0(indir, "/Rdata/diet.fudan.rds")) %>% mutate(eid = as.character(eid))
need.raw <- c("eid", grep(paste0("^(", paste(raw.items, collapse = "|"), ")_i[0-4]$"), names(diet.raw), value = TRUE))

dat.raw <- diet.raw[, need.raw, drop = FALSE] %>% select(eid)
for (v in raw.items) dat.raw[[v]] <- collapse1(diet.raw, v)

dat.item0 <- dat.raw %>% left_join(dat.food %>% select(eid, any_of(derived.items), energy.kJ), by = "eid")

dash.items <- unique(c(field %>% filter(DASH > 0) %>% pull(name), "w.p26052"))
maha.items <- unique(c(dairy.vars, veg.vars, fruit.vars, fat.vars, wholegrain.vars, refined.vars, ssb.vars, sweet.vars, hp.vars, alcohol.vars, "w.p26005", "w.p26052"))

extra.dict <- tribble(
	~item, ~item_label,
	"oil", "Healthy oil", "bread.white", "White bread", "bread.wholemeal", "Wholemeal bread",
	"bread.seeded", "Seeded bread", "baguette.white", "White baguette", "baguette.wholemeal", "Wholemeal baguette",
	"bap.white", "White bap", "bap.wholemeal", "Wholemeal bap", "roll.white", "White roll",
	"roll.wholemeal", "Wholemeal roll", "w.p26005", "Protein intake", "w.p26052", "Sodium"
)

item.map <- tibble(item = union(raw.items, derived.items)) %>%
	filter(item %in% names(dat.item0)) %>%
	mutate(
		in_dash = item %in% dash.items,
		in_maha = item %in% maha.items,
		origin = case_when(in_dash & in_maha ~ "Both", in_dash ~ "DASH only", in_maha ~ "MAHA only", TRUE ~ "Other"),
		origin2 = case_when(
			item %in% alcohol.vars ~ "Alcohol",
			item %in% dairy.vars ~ "Dairy",
			item %in% meat.vars ~ "Meat",
			item %in% ssb.vars ~ "Drinks",
			item %in% refined.vars ~ "Refined grains",
			item %in% c(sweet.vars, setdiff(hp.vars, meat.vars)) ~ "Sweets/processed foods",
			item %in% fruit.vars ~ "Fruit",
			item %in% veg.vars ~ "Vegetables",
			item %in% wholegrain.vars ~ "Whole grains",
			item %in% fat.vars ~ "Fats/oils",
			item == "w.p26005" ~ "Protein",
			item == "w.p26052" ~ "Sodium",
			TRUE ~ "Other"
		)
	) %>%
	left_join(diet.dict, by = "item") %>%
	left_join(extra.dict, by = "item", suffix = c("", ".x")) %>%
	mutate(item_label = coalesce(item_label.x, item_label, item)) %>%
	select(-item_label.x)

dat.item <- dat %>%
	transmute(eid = as.character(eid), maha = diet.maha.sum, dash = diet.dash.sum, across(any_of(covs))) %>%
	mutate(
		maha.hml = mk_hml(maha, group_pct_th),
		dash.hml = mk_hml(dash, group_pct_th),
		disc.grp = case_when(maha.hml == "high" & dash.hml == "low" ~ "MH_DL", maha.hml == "low" & dash.hml == "high" ~ "ML_DH", TRUE ~ NA_character_)
	) %>%
	filter(!is.na(disc.grp)) %>%
	left_join(dat.item0, by = "eid") %>%
	mutate(y_disc = as.integer(disc.grp == "MH_DL"))

res.item.disc <- purrr::map_dfr(item.map$item, fit_disc, dat.item = dat.item) %>% left_join(item.map, by = "item") %>% arrange(desc(abs(beta)))
top.items <- res.item.disc %>% slice_max(order_by = abs(beta), n = 20) %>% pull(item)
res.item.cox <- purrr::map_dfr(top.items, \(v) purrr::map_dfr(Y.inc, \(Y) fit_cox(v, Y, dat.item0))) %>% left_join(item.map, by = "item")
origin2.levels <- c("Alcohol", "Dairy", "Drinks", "Fats/oils", "Fruit", "Vegetables", "Meat", "Protein", "Refined grains", "Whole grains", "Sweets/processed foods", "Sodium", "Other")

item.class <- item.map %>%
	filter(item %in% top.items) %>%
	select(item, origin2) %>%
	mutate(origin2 = factor(origin2, levels = origin2.levels))

ord.tbl <- res.item.disc %>%
	filter(item %in% top.items) %>%
	select(item, item_label, beta.disc = beta) %>%
	distinct() %>%
	left_join(item.class, by = "item") %>%
	arrange(origin2, desc(beta.disc)) %>%
	mutate(item_label = factor(item_label, levels = rev(unique(item_label))))

res.item.cox2 <- res.item.cox %>%
	filter(item %in% top.items) %>%
	select(item, Outcome, HR, p.value, beta, item_label) %>%
	left_join(ord.tbl %>% select(item, item_label.ord = item_label, origin2), by = "item") %>%
	mutate(
		item_label = item_label.ord,
		origin2 = factor(as.character(origin2), levels = origin2.levels),
		Outcome = factor(Outcome, levels = dx.lst[Y.inc]),
		lab = ifelse(p.value < 0.05, sprintf("%.2f", HR), "")
	)

p4 <- ggplot(res.item.cox2, aes(Outcome, item_label, fill = beta)) +
	geom_tile(color = "white", linewidth = 0.5) +
	geom_text(aes(label = lab), size = 3.0) +
	scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
	facet_grid(origin2 ~ ., scales = "free_y", space = "free_y", switch = "y", drop = TRUE) +
	labs(x = NULL, y = NULL, fill = "log(HR)") +
	theme_bw(base_size = 13) +
	theme(
		strip.placement = "outside",
		strip.background = element_rect(fill = "grey95", color = NA),
		strip.text.y.left = element_text(angle = 0, face = "bold"),
		axis.text.x = element_text(angle = 35, hjust = 1),
		panel.spacing.y = unit(0.15, "lines")
	)

save_plot(p4, "FigS4.png", width = 10.5, height = 8.4, dpi = 320)

save_xlsx1(
	list(
		item_map = item.map,
		item_discordance = res.item.disc,
		top_items = data.frame(item = top.items),
		ordering_table = ord.tbl %>% as.data.frame(),
		top_item_outcome = res.item.cox,
		top_item_outcome_heatmap = res.item.cox2 %>% as.data.frame()
	),
	"FigS4.out.xlsx"
)
print(head(res.item.disc, 50))
print(res.item.cox, n = 50)
