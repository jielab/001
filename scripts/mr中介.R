#包的安装与运行

install.packages("remotes")
remotes::install_github("cran/RMediation",force = TRUE)


library(RMediation)

library(tidyr)
library(readr)
library(dplyr)
library(kableExtra)

library(ggplot2)
library(cowplot)


#中介：中介：TV(ukb-b-5192)到吸烟(ieu-b-25)到LC(ieu-a-966)

#第一步：TV到LC，总效应C

library(TwoSampleMR)

TV<- extract_instruments(outcomes = 'ukb-b-5192')

LC <- extract_outcome_data(snps = TV$SNP, outcomes = 'ieu-a-966')

dat <- harmonise_data(TV, LC)

res <- mr(dat)

exposure_total_beta <- res %>%filter(method == "Inverse variance weighted") %>%pull(b)
exposure_total_se <- res %>%filter(method == "Inverse variance weighted")%>% pull(se)



#差异法

#总效应=直接效应+中介效应
#直接效应=TV到LC的直接效应 (将TV和吸烟纳入暴露做多变量)

mvmr<-mv_extract_exposures(c("ukb-b-5192","ieu-b-25"), clump_r2=0.001, clump_kb=10000, harmonise_strictness=2, access_token = ieugwasr::check_access_token(), find_proxies=TRUE, force_server=FALSE, pval_threshold=5e-8, pop="EUR") #提取TV与吸烟的工具变量

mvmr_outcome_dat<-extract_outcome_data(mvmr$SNP, "ieu-a-966") #提取LC中的结局信息

mvmr_dat<-mv_harmonise_data(mvmr, mvmr_outcome_dat, harmonise_strictness=2) #合并

mvmr_res<-mv_multiple(mvmr_dat) #MVMR结果

#TV到LC的直接效应

direct_beta<-mvmr_res[["result"]][["b"]][2] #提取出beta值

direct_se<-mvmr_res[["result"]][["se"]][2] #提取出se

difference_method_PoE <- function(total_beta, total_se, direct_beta, direct_se, verbose = F){
  # 计算中介效应
  
  # 计算中介效应的beta（中介效应=总效应-直接效应）
  
indirect_beta = total_beta -  direct_beta
  
  #indirect_beta = round(indirect_beta,2)  将变量 indirect_beta 的值四舍五入到小数点后两位
  if (verbose) {print(paste("Indirect effect = ", 
                            round(total_beta, 2)," - ", round(direct_beta, 2), 
                            " = ", round(indirect_beta,2)))}
  
  
  # 计算中介效应的se 
  ### 误差传染法
  # SE of INDIRECT effect (difference) = sqrt(SE TOTAL^2 + SE DIRECT^2) 
  indirect_se = round(sqrt(total_se^2 + direct_se^2), 4)
  if (verbose) {print(paste("SE of indirect effect = sqrt(",
                            round(total_se, 2),"^2 + ", round(direct_se,2), 
                            "^2) = ", indirect_se))}
  
  
  # 把数据整理成df，数据整理成符合"tidy data"规范的数据框（data frame），即每个变量形成一列，每个观测形成一行，每个单元格中只有一个值
  df <-data.frame(b= indirect_beta,
                  se = indirect_se)
  
  # 计算置信区间 CIs 和 OR
  df$lo_ci    <- df$b - 1.96 * df$se
  df$up_ci    <- df$b + 1.96 * df$se
  df$or        <-  exp(df$b)
  df$or_lci95 <- exp(df$lo_ci)
  df$or_uci95 <- exp(df$up_ci)
  
  df<-round(df,3)
  return(df)
}

indirect_effect<-difference_method_PoE(exposure_total_beta,exposure_total_se,direct_beta,direct_se) #将数据带入计算 
indirect_effect  #查看结果



#中介效应比
indirect_effect[,1]/exposure_total_beta

#乘积法

#两种办法算乘积法中介效应A*B
## INDIRECT = TOTAL (exposure -> mediator) x TOTAL (mediator -> outcome)
## INDIRECT = TOTAL (exposure -> mediator) x DIRECT (of mediator , mvmr) 

#TV暴露到中介吸烟的效应,A 两样本
TV<- extract_instruments(outcomes = 'ukb-b-5192')

outcome_dat_M <- extract_outcome_data(snps = TV$SNP, outcomes = 'ieu-b-25')

dat_M <- harmonise_data(TV, outcome_dat_M)

res_M <- mr(dat_M)


EM_beta<- res_M %>%filter(method == "Inverse variance weighted") %>%pull(b)

EM_se<- res_M %>%filter(method == "Inverse variance weighted") %>%pull(se)


#中介吸烟到结局LC的效应,B 两样本
mediation<- extract_instruments(outcomes = 'ieu-b-25')
#去掉A的IVs
outcome_dat_M0 <- extract_outcome_data(snps = mediation$SNP, outcomes = 'ieu-a-966')

dat_M0 <- harmonise_data(mediation, outcome_dat_M0)

res_M0 <- mr(dat_M0)

MO_beta_total<- res_M0 %>%filter(method == "Inverse variance weighted") %>%pull(b)

MO_se_total<- res_M0 %>%filter(method == "Inverse variance weighted") %>%pull(se)



#中介到结局的直接效应,B

MO_beta<-mvmr_res[["result"]][["b"]][1]

MO_se<-mvmr_res[["result"]][["se"]][1]



#两种办法算乘积法的se
#  1) Delta法
#  2) 误差传染法


#  1) Delta
product_method_Delta <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
  
  
  # method 1
  # INDIRECT = TOTAL (exposure -> mediator) x TOTAL (mediator -> outcome)
  # method 2
  # INDIRECT = TOTAL (exposure -> mediator) x DIRECT (of mediator , mvmr) 
  
  
  # 计算中介效应beta
  EO <- EM_beta * MO_beta
  
  if (verbose) {print(paste("Indirect effect = ", round(EM_beta, 2)," x ", round(MO_beta,2), " = ", round(EO, 3)))}
  
  
  #  RMediation package计算置信区间
  CIs = medci(EM_beta, MO_beta, EM_se, MO_se, type="dop")
  
  # 把数据 放进数据集
  df <-data.frame(b = EO,
                  se = CIs$SE,
                  lo_ci = CIs[["95% CI"]][1],
                  up_ci= CIs[["95% CI"]][2])
  # 计算 OR
  df$or        <-  exp(df$b)
  df$or_lci95 <- exp(df$lo_ci)
  df$or_uci95 <- exp(df$up_ci)
  
  df<-round(df,3)
  return(df)
}

#第一种INDIRECT = TOTAL (exposure -> mediator) x TOTAL (mediator -> outcome)
product_method_Delta(EM_beta, EM_se, MO_beta_total, MO_se_total)

#第二种 INDIRECT = TOTAL (exposure -> mediator) x DIRECT (of mediator , mvmr)
product_method_Delta(EM_beta, EM_se, MO_beta, MO_se)

#  2) 误差传染法
product_method_PoE <- function(EM_beta, EM_se, MO_beta, MO_se, verbose=F){
  
  # 计算中介效应
  EO_beta <- EM_beta * MO_beta
  
  if (verbose) {print(paste("Indirect effect = ", round(EM_beta, 2)," x ", round(MO_beta,2), " = ", round(EO, 3)))}
  
  
  # 计算中介效应se
  ### 误差传染法
  # SE of INDIRECT effect (difference) = sqrt(SE EM^2 + SE MO^2) 
  EO_se = round(sqrt(EM_se^2 + MO_se^2), 4)
  if (verbose) {print(paste("SE of indirect effect = sqrt(",
                            round(EM_se, 2),"^2 + ", round(MO_se,2), 
                            "^2) = ", indirect_se))}
  
  
  # put data into a tidy df
  df <-data.frame(b= EO_beta,
                  se = EO_se)
  
  # calculate CIs and OR
  df$lo_ci    <- df$b - 1.96 * df$se
  df$up_ci    <- df$b + 1.96 * df$se
  df$or        <-  exp(df$b)
  df$or_lci95 <- exp(df$lo_ci)
  df$or_uci95 <- exp(df$up_ci)
  
  df<-round(df,3)
  return(df)
}

#第一种 中介效应 =两样本(暴露 -> 中介) x 两样本 (中介 -> 暴露)
product_method_PoE(EM_beta, EM_se, MO_beta_total, MO_se_total)

#第二种 中介效应 = 两样本 (暴露 -> 中介) x 多变量 (中介→暴露)
product_method_PoE(EM_beta, EM_se, MO_beta, MO_se)

