setwd("D:/101/female")
pacman::p_load(data.table, readxl, lubridate, tidyverse, plotly, dplyr, survival, survminer, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival分析入门 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data(cancer, package="survival")
lung$sex <- as.factor(lung$sex)
lung %>% group_by(status, sex) %>% summarize(average = mean(time)) #tally()
surv.obj <- Surv(time=lung$time, event=lung$status) 
fit.surv <- survfit(surv.obj ~ sex, data=lung)
	ggsurvplot(fit.surv, ylim=c(0,1), pval=TRUE, size=2, conf.int=FALSE, ggtheme=theme_classic(), risk.table=FALSE, data=lung)
fit.cox <- coxph(surv.obj ~ ph.ecog+ph.karno+pat.karno+meal.cal+wt.loss +age+sex, data=lung)
	summary(fit.cox); ggforest(fit.cox)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 实战：读入数据并生成几个大家都用的变量 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") %>% 
	mutate (
	a = ifelse(abo=="A", "A", "non-A"),
	o = ifelse(abo=="O", "O", "non-O"),
	se = ifelse(fut2.rs601338_A==2, "non-se", "se"), # 两个A是“非分泌性”
	a_se = factor(paste(a, se, sep="."), levels=c("non-A.se", "non-A.non-se", "A.se", "A.non-se")), # 把最多的组放在前面，作为ref
	o_se = factor(paste(o, se, sep="."), levels=c("non-O.se", "non-O.non-se", "O.se", "O.non-se")), # 把最多的组放在前面，作为ref
	s = ifelse(sp1.S==0, "non-S", "S"),
	z = ifelse(sp1.Z==0, "non-Z", "Z"),
	shbg.rs6257_T = hardcall(shbg.rs6257_T),
	rh.RHD.rs590787_A = hardcall(rh.RHD.rs590787_A)
	) %>% 
	dplyr::select(-c(lct.microbe.rs3940549_A, lct.microbe.rs182549_T)) # 这两个SNP跟lct.MCM6.rs4988235_A高度相关：r>0.9


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 数据 Sanity check 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=TRUE)), facet=sex)
grep("abo", names(dat0), value=TRUE) # 找变量名字
attach(dat0) #这样，下面就可以不用写 dat0$，直接写变量名字就行
	hist(icdDate_covid, breaks="weeks", freq=TRUE); table(icd_covid)
	table(abo, abo.O1.rs8176719_T); table(abo, abo.AB.rs8176746_G) # ABO血型
	cor(abo.O1.rs8176719_T, abo.microbe.rs545971_C, use="complete.obs") # 0.93, (abo.microbe.rs550057_C) 0.80, (abo.microbe.rs8176645_T) 0.87
	prop.table(table(abo, fut2.rs601338_A),1) 
	table(sport.ACTN3.rs1815739_C); table(sport.ACE.rs4343_G); hist(sport.RBFOX1.rs7191721_G) 
	table(icd_lungcancer); hist(icdDate_lungcancer, breaks="months")
	group_by(dat, Y_yes, o_se) %>% summarise(count=n(), mean=mean(X, na.rm=TRUE))
	aggregate(X ~ Y_yes*o_se, dat, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	ggplot(dat, aes(x=Y_yes, y=X, fill=o_se)) + geom_boxplot() + stat_summary(fun.y=mean, color="darkred", position=position_dodge(0.75), geom="point", shape=18, size=3)
	bp <- boxplot(X ~ Y_yes*se, dat, las=2, col=rainbow(6), font=2); bp$stats
detach(dat0)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=TRUE)) %>% na.omit() %>% dplyr::sample_n(10000)
	car::scatterplotMatrix(dat, spread=FALSE, smoother.args=list(lty=0.1))
	psych::pairs.panels(dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 北大公卫数据 Sanity check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PKU_le8.xlsx sex +dashpts +PA_pts +smoke_pts +sleep_pts +bmi_pts +nonhdl_pts +hba1c_pts +BP_pts +LE8score +CVH_cat
# rs1815739 in ACTN3, rs7191721 in RBFOX1
dir <- "D:/data/ukb/phe/"
pku <- read.csv(paste0(dir,"PKU_lungcancer.csv"))
	surv.obj <- Surv(time=pku$lung_cancer_time, event=pku$lung_cancer_inc)
	fit.cox <- coxph(surv.obj ~ walkingpace + rs1815739_C + walkingpace*rs1815739_C +age+sex+centre+Townsend_i+alcohol+Smoking+pa_3c+HDS2_i +BMI_i+grip_i+hpt_baseline+diabetes_baseline_first_report+cvd_baseline+famhis_lungcancer+genebatch+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data=pku)
	summary(fit.cox)  
pku <- dat0 %>% merge(pku, by.x="eid", by.y="n_eid") %>% rename(age=age.x, sex=sex.x) 
	nrow(pku); grep("\\.y", names(pku), value=TRUE)
	table(pku$walking_pace, pku$walkingpace) # walkingpace来自北大的数据
	summary(pku$walking_freq); summary(pku$walking_time)
	plot(pku$walking_freq, pku$walking_time)
	summary(pku$Townsend_i); summary(pku$deprivation)
	plot(pku$Townsend_i, pku$deprivation)
dat1 <- pku %>% # 继续比较survival分析用到的变量
	mutate(
		Y_date = icdDate_lungcancer,
		Y_yes = ifelse( is.na(icdDate_lungcancer), 0,1),
		time2lungcancer = as.numeric(Y_date) - as.numeric(date_attend),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_year = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_year >0)
	table(dat1$Y_yes, dat1$lung_cancer_inc)
	summary(dat1$icdDate_lungcancer); summary(dat1$time2lungcancer); summary(dat1$lung_cancer_time)
	plot(dat1$follow_year, dat1$lung_cancer_time) # 对不上
		gg <- ggplot(data=dat1, aes(x=follow_year, y=lung_cancer_time)) + geom_point() # + stat_summary(fun="count") ??????????????
		ggplotly(gg)
	plot(dat1$time2lungcancer, dat1$lung_cancer_time) # 看着对，但是数据很少
		hist(dat1$time2lungcancer); hist(dat1$lung_cancer_time)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y-Z 批量分析示例
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% drop_na(abo, se, apoe) %>%
	filter(ethnic_cat=="White" & sex==1) %>% # ethnicity_gen==1
	mutate(
		walk = inormal(walking_time * walking_freq * walking_pace),
		across(grep("walking", names(dat0), value=T), ~factor(.x))
	) 
Xs <- c("walking_pace", "bb_ALP", "bb_AST", "bb_BUN", "bb_CRE", "bb_CRP", "bb_CYS", "bb_LPA", "bb_oestradiol", "bb_SHBG", "bb_TES") # grep("^bb_", names(dat), value=TRUE)
Ys <- c("icdDate_lungcancer", "icdDate_t2dm", "icdDate_chd", "icdDate_stroke", "icdDate_asthma", "icdDate_copd") # grep("^icdDate", names(dat), value=TRUE)
Zs <- grep("^o$|^se$|^rh|shbg|^apoe$|\\.rs", names(dat), value=TRUE) #  
outfile="res.txt"; sink(outfile)
for (Y in Ys) {
	print(Y)
	dat1 <- dat %>%
	mutate(
		Y_date = dat[[Y]],
		Y_yes = ifelse( is.na(dat[[Y]]), 0,1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	for (X in Xs) {
		dat1$X = dat1[[X]]
		for (Z in Zs) {
			dat1$Z <- dat1[[Z]]
			#fit.glm <- glm(Y_yes ~ X +Z +X*Z +age+sex+bmi+education+deprivation +PC1+PC2, data=dat1, family="binomial")
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +shbg.rs6259_G +age+sex+bmi +PC1+PC2, data=dat1)
			res.cox <- coef(summary(fit.cox))
			p_int <- signif(tail(res.cox,1)[,5] ,2)
			if (!is.na(p_int) & p_int < 0.05) {
				print(paste(X, Y, Z, p_int)) 
			}
		}
	}
}
sink()
# 如果上面代码用了write.table，可用下面代码整合分析结果
outfile2="surv.p.tsv"; file.create(outfile2)
dat <- read.table(outfile, sep='|', header=FALSE, as.is=TRUE) 
	names(dat) <- c('X', 'Y', 'cnt', 'b', 'se', 'p')
	dat$p=signif(dat$p,2)
	dat.b <- dat %>% dplyr::select(X, Y, b) %>% acast(X ~ Y, value.var='b'); dat.b=round(dat.b,1)
	dat.p <- dat %>% dplyr::select(X, Y, p) %>% acast(X ~ Y, value.var='p')
	write.table(dat.p, file=outfile2, sep='\t', row.names=TRUE, col.names=TRUE, append=FALSE, quote=FALSE)
	#plt <- ggcorrplot(dat.b, lab=TRUE, p.mat=dat.p, sig.level=1e-4, insig ='blank') + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某一*显著*结果的精细分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X="bb_TES"
Y="icdDate_t2dm"
Z="se"
dat1 <- dat %>% 
	mutate(
		X = dat[[X]],
		X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low","middle","high")),
		Y_date = dat[[Y]],
		Y_yes = ifelse(is.na(dat[[Y]]), 0, 1),
		Y_status = ifelse(!is.na(Y_date), "Disease", ifelse(!is.na(date_lost), "Lost", ifelse(!is.na(death_date), "Decease", "OK"))),
		follow_end_day = data.table::fifelse(!is.na(Y_date), Y_date, data.table::fifelse(!is.na(death_date), death_date, as.Date("2021-10-01"))), # fifelse preserves the type and class of the inputs.
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
		follow_years_int = ceiling(follow_years),
		Z = dat[[Z]]
	) %>% filter( follow_years >0 ) 
plot(dat1$date_attend, dat1$follow_end_day)
table(dat1$Y_yes); table(dat1$Death_yes); table(dat1$Y_status, dat1$follow_years_int)
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	surv.obj[1:20]
	fit.surv <- survfit(surv.obj ~1); str(fit.surv)
	fit.surv %>% ggsurvfit::ggsurvfit() + add_confidence_interval() + add_risktable() #Kaplan-Meier curve
	ggsurvplot(survfit(surv.obj ~1, data=dat1), ylim=c(0.96,1),risk.table=TRUE)
	ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), palette=c("green","green4", "gold","gold4", "tomato","tomato4"), ylim=c(0.8,1), conf.int=FALSE, pval=FALSE, pval.method=TRUE, test.for.trend=FALSE, surv.median.line="hv", risk.table=TRUE, cumevents=FALSE)
	dat1 %>% ggscatter(x=X, y=Y, color=Z, xlab=X, ylab=Y, add="reg.line", ellipse=TRUE, conf.int=FALSE, mean.point=TRUE)
fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+bmi+PC1+PC2, data=dat1); coef(summary(fit.cox))
	res.cox <- coef(summary(fit.cox)) %>% as.data.frame(row.names=dimnames(.)[[1]]) %>% data.frame(X=row.names(.), row.names=NULL) 
	data.table::setDT(res.cox, keep.rownames="X")
	survminer::ggforest(fit.cox, main="", fontsize=1.2, data=dat1)
	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	res.cox as.data.frame(res.cox, row.names=dimnames(res.cox)[[1]])
	res.add <- rbind(res.cox) ## 把多个分析结果综合起来
	ggforestplot::forestplot(df=df, name=name, estimate=beta, se=se, pvalue=pvalue, psignif=0.002, colour=study, shape=study, xlab="", title="", logodds=TRUE) + 
		ggplot2::scale_shape_manual(values = c(23L, 21L, 21L, 21L, 21L), labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS"))
	
