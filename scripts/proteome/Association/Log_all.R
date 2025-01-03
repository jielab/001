library(data.table)
library(dplyr)
library(parallel)
library(stats)

results_summary <- function(tgt_out_df) {
  oratio_out_lst <- sapply(1:nrow(tgt_out_df), function(i) {
    oratio <- sprintf("%.2f", tgt_out_df$oratio[i])
    lbd <- sprintf("%.2f", tgt_out_df$or_lbd[i])
    ubd <- sprintf("%.2f", tgt_out_df$or_ubd[i])
    paste0(oratio, " [", lbd, "-", ubd, "]")
  })
  
  p_out_lst <- sapply(1:nrow(tgt_out_df), function(i) {
    if (tgt_out_df$pval_bfi[i] < 0.001) {
      "***"
    } else if (tgt_out_df$pval_bfi[i] < 0.01) {
      "**"
    } else if (tgt_out_df$pval_bfi[i] < 0.05) {
      "*"
    } else {
      ""
    }
  })
  
  list(oratio_out_lst = oratio_out_lst, p_out_lst = p_out_lst)
}

process <- function(pro_f, cov_f_lst, tmp_tgt_df, pro_cov_df) {
  tmp_pro_df <- pro_cov_df[, c("eid", pro_f, paste0(pro_f, "_SampAge"))]
  tmp_df <- merge(tmp_tgt_df, tmp_pro_df, by = "eid", all.x = TRUE)
  colnames(tmp_df) <- gsub(pro_f, "x_pro", colnames(tmp_df))
  tmp_df$x_pro_sa[is.na(tmp_df$x_pro_sa)] <- median(tmp_df$x_pro_sa, na.rm = TRUE)
  
  tryCatch({
    Y <- tmp_df$target_y
    X <- tmp_df[, c(cov_f_lst, "x_pro", "x_pro_sa"), with = FALSE]
    log_mod <- glm(Y ~ ., data = X, family = binomial(link = "logit"))
    oratio <- exp(coef(log_mod)["x_pro"])
    ci_mod <- confint(log_mod, parm = "x_pro", level = 0.95)
    lbd <- exp(ci_mod[1])
    ubd <- exp(ci_mod[2])
    pval <- summary(log_mod)$coefficients["x_pro", "Pr(>|z|)"]
    nb_all <- nrow(tmp_df)
    nb_case <- sum(tmp_df$target_y)
    prop_case <- round(nb_case / nb_all * 100, 3)
    c(pro_f, nb_all, nb_case, prop_case, oratio, lbd, ubd, pval)
  }, error = function(e) {
    c(pro_f, nb_all, nb_case, prop_case, NA, NA, NA, NA)
  })
}

sort_nicely <- function(l) {
  l[order(as.numeric(gsub("[^0-9]", "", l)))]
}

# Paths and Data Loading
dpath <- "/home1/jiayou/Documents/Projects/ProDisAtlas/"
target_file_lst <- sort_nicely(list.files(paste0(dpath, "Data/Target/Targets2Analysis/"), full.names = TRUE))
target_info_df <- fread(paste0(dpath, "Data/Target/TargetVsProtein.csv"))
pro_cov_df <- fread(paste0(dpath, "Data/ProteinData/ProteinData_n_Cov.csv"))

pro_f_lst <- colnames(pro_cov_df)[14:2933]
pro_cov_df$Race <- ifelse(pro_cov_df$Race == 1, 1, 0)

cov_f_lst_in_sex <- c("Age", "Sex", "Race", "TDI", "BMI", "smk", "fastingtime", "season")
cov_f_lst_non_sex <- c("Age", "Race", "TDI", "BMI", "smk", "fastingtime", "season")
cov_df <- pro_cov_df[, c("eid", cov_f_lst_in_sex)]

bad_tgt <- list()
cl <- makeCluster(detectCores() - 1)
clusterExport(cl, list("process", "results_summary", "cov_f_lst_in_sex", "cov_f_lst_non_sex", "pro_cov_df"))

for (tgt_file in target_file_lst) {
  tgt_name <- basename(tgt_file)
  tgt_info <- target_info_df[target_info_df$NAME == tgt_name, ]
  cov_f_lst <- ifelse(tgt_info$SEX %in% c(1, 2), cov_f_lst_non_sex, cov_f_lst_in_sex)
  
  tmp_tgt_df <- fread(tgt_file)[, .(eid, target_y, BL2Target_yrs)]
  tmp_tgt_df <- tmp_tgt_df[tmp_tgt_df$BL2Target_yrs <= 0 | tmp_tgt_df$target_y == 0]
  tmp_tgt_df <- merge(tmp_tgt_df, cov_df, by = "eid", all.x = TRUE)
  
  if (sum(tmp_tgt_df$target_y) > 50) {
    res <- parLapply(cl, pro_f_lst, process, cov_f_lst, tmp_tgt_df, pro_cov_df)
    results <- rbindlist(res)
    colnames(results) <- c("Pro_code", "nb_individuals", "nb_case", "prop_case(%)", "oratio", "or_lbd", "or_ubd", "pval_raw")
    _, p_f_bfi <- p.adjust(results$pval_raw, method = "bonferroni")
    results$pval_bfi <- p_f_bfi
    results$pval_bfi[results$pval_bfi >= 1] <- 1
    summary <- results_summary(results)
    results$or_output <- summary$oratio_out_lst
    results$pval_significant <- summary$p_out_lst
    fwrite(results, paste0(dpath, "Results/Association/CrossSectionalAnalysis/All/", tgt_name, ".csv"))
  } else {
    bad_tgt <- append(bad_tgt, list(list(tgt_name, sum(tmp_tgt_df$target_y))))
  }
}

stopCluster(cl)
bad_tgt_df <- as.data.frame(do.call(rbind, bad_tgt))
fwrite(bad_tgt_df, paste0(dpath, "Results/Association/CrossSectionalAnalysis/bad_targets_cs_all.csv"))
