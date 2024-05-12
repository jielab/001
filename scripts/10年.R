setwd("C:/Users/Chong You/Desktop/pku/HuangJie/walkingpace+VTE")
library(survival); library(ggplot2)

dat <- readRDS("ukb.final.rds")
cox_model <- coxph(Surv(follow_years, Y_yes) ~ X+Z, data = dat)
#	surv_pred <- survfit(cox_model, newdata = new_data)
#	surv_summary <- summary(surv_pred, times = 10)
	base_haz <- basehaz(cox_model, centered = FALSE)
	cum_base_haz <- approx(base_haz$time, base_haz$hazard, xout = 10, method = "linear")$y
	new_data <-  expand.grid(X = levels(dat$X), Z = levels(dat$Z))
	new_data_num <- model.matrix(~ X + Z, data = new_data)[,-1]
	Ht <- exp(new_data_num %*% cox_model$coefficients) * cum_base_haz
	St <- exp(-Ht)
	var_beta_x <- new_data_num%*%diag(cox_model$var) 
	se_Ht <- sqrt(var_beta_x) * Ht
	se_St <- St * se_Ht
	z_value <- qnorm(0.975)  # 正态分布的97.5%分位数
	lower_CI <- St * exp(-z_value * se_St)
	upper_CI <- St * exp(z_value * se_St)
res <- data.frame(X = new_data$X, Z = new_data$Z, SurvivalProb = St, LowerCI = lower_CI, UpperCI = upper_CI)
	ggplot(res, aes(x = Z, y = SurvivalProb, fill = X)) +
	geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.7) +
	geom_errorbar(aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = position_dodge(width = 0.7)) +
	geom_text(aes(label = sprintf("%.2f", UpperCI), y = UpperCI), vjust = -0.5, position = position_dodge(width = 0.7), size = 2) +
	geom_text(aes(label = sprintf("%.2f", LowerCI), y = LowerCI), vjust = 1.5, position = position_dodge(width = 0.7), size = 2) +
	labs(x = "Genetic Factor", y = "Survival Probability", title = "Survival Probabilities by Walking Pace and Genetic Factor") +
	scale_fill_brewer(palette = "Set1", name = "Walking Pace") + theme_minimal() + coord_cartesian(ylim = c(0.915, 0.99))
