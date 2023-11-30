setwd("D:/101/female")
pacman::p_load(data.table, readxl, lubridate, tidyverse, plotly, dplyr, survival, survminer, ggsurvfit, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读入数据，生成几个新变量，sanity check一下 
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
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=TRUE)), facet=sex)
attach(dat0) #这样，下面就可以不用写 dat0$，直接写变量名字就行
	hist(icdDate_covid, breaks="weeks", freq=TRUE); table(icd_covid)
	table(abo, abo.O1.rs8176719_T); table(abo, abo.AB.rs8176746_G) # ABO血型
	cor(abo.O1.rs8176719_T, abo.microbe.rs545971_C, use="complete.obs") # 0.93, (abo.microbe.rs550057_C) 0.80, (abo.microbe.rs8176645_T) 0.87
	prop.table(table(abo, fut2.rs601338_A),1) 
	group_by(dat0, abo, se) %>% summarise(count=n(), mean=mean(bb_TES, na.rm=TRUE))
	aggregate(bb_TES ~ abo*se, dat0, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	bp <- boxplot(bb_TES ~ abo*se, dat0, las=2, col=rainbow(6), font=2); bp$stats
	dat0 %>% drop_na(abo, se) %>% ggplot(aes(x=abo, y=bb_TES, fill=se)) + geom_boxplot() + stat_summary(fun.y=mean, color="darkred", position=position_dodge(0.75), geom="point", shape=18, size=3)
detach(dat0)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=TRUE)) %>% na.omit() %>% dplyr::sample_n(10000)
	car::scatterplotMatrix(dat, spread=FALSE, smoother.args=list(lty=0.1))
	psych::pairs.panels(dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 北大公卫数据 Sanity check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PKU_le8.xlsx sex +dashpts +PA_pts +smoke_pts +sleep_pts +bmi_pts +nonhdl_pts +hba1c_pts +BP_pts +LE8score +CVH_cat
# 文章中的两个SNP：power.ACTN3.rs1815739_C in ACTN3, rs7191721_G in RBFOX1
dir <- "D:/data/ukb/phe/"
pku <- read.csv(paste0(dir,"PKU_lungcancer.csv"))
	surv.obj <- Surv(time=pku$lung_cancer_time, event=pku$lung_cancer_inc)
	fit.cox <- coxph(surv.obj ~ walkingpace + rs1815739_C + walkingpace*rs1815739_C +age+sex+centre+Townsend_i, data=pku); coef(summary(fit.cox))  
pku <- dat0 %>% merge(pku, by.x="eid", by.y="n_eid", all.y=T) %>% rename(age=age.x, sex=sex.x) 
	nrow(pku); grep("\\.y", names(pku), value=TRUE); table(pku$ethnic_cat)
	table(pku$walking_pace, pku$walkingpace, useNA="always") # walkingpace来自北大的数据
	plot(pku$Townsend_i, pku$deprivation); table(!is.na(pku$Townsend_i), !is.na(pku$deprivation))
dat1 <- pku %>% filter(sex==1) %>% # 继续比较survival分析用到的变量
	mutate(
		Y_date = icdDate_lungcancer,
		Y_yes = ifelse( is.na(icdDate_lungcancer), 0,1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_year = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) 
	table(dat1$Y_yes, dat1$lung_cancer_inc, useNA="always")
	plot(dat1$follow_year, dat1$lung_cancer_time) # 基本对上了
		table(!is.na(dat1$follow_year), !is.na(dat1$lung_cancer_time))
		gg <- ggplot(data=dat1, aes(x=follow_year, y=lung_cancer_time)) + geom_point()
		ggplotly(gg)
	surv.obj <- Surv(time=dat1$follow_year, event=dat1$Y_yes)
	fit.cox <- coxph(surv.obj ~ walking_pace + endurance.RBFOX1.rs7191721_G  + walking_pace*endurance.RBFOX1.rs7191721_G  +age+sex+centre+deprivation, data=dat1); coef(summary(fit.cox))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y-Z 批量分析示例
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% drop_na(abo, se, apoe) %>%
	filter(ethnic_cat=="White" & sex==0) %>% # ethnicity_gen==1
	mutate(
		walk = inormal(walking_time * walking_freq * walking_pace),
		across(grep("walking", names(dat0), value=T), ~factor(.x))
	) 
Xs <- c("bb_ALP", "bb_CRE", "bb_CRP", "bb_CYS", "bb_LPA", "bb_oestradiol", "bb_SHBG", "bb_TES") # grep("^bb_", names(dat), value=TRUE)
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
			dat1$X_qt <- cut(dat1$X, breaks=quantile(dat1$X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5))
			dat1$X_qt <- factor(ifelse(dat1$X_qt=="q1", "low", ifelse(dat1$X_qt=="q5", "high", "middle")), levels=c("low","middle","high"))
			#fit.glm <- glm(Y_yes ~ X +Z +X*Z +age+sex+bmi+education+deprivation +PC1+PC2, data=dat1, family="binomial")
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+bmi +PC1+PC2, data=dat1)
			res.cox <- coef(summary(fit.cox))
			p_int <- signif(tail(res.cox,1)[,5] ,2)
			if (!is.na(p_int) & p_int < 0.05) {
				print(paste(X, Y, Z, p_int)) 
				png(paste(X,Y,Z,"png",sep="."))
				if (Z %like% "\\.rs" & length(unique(dat1$Z)) >3) dat1$Z <- hardcall(dat1$Z) # 要不然没法画出Z的3个类型
				print(ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), ylim=c(0.8,1), risk.table=FALSE))	
				dev.off()
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
Z="o.se"
dat1 <- dat %>% 
	mutate(
		X = dat[[X]],
		X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low","middle","high")),
		Y_date = dat[[Y]],
		Y_yes = ifelse(is.na(dat[[Y]]), 0, 1),
		Y_status = ifelse(!is.na(Y_date), "Disease", ifelse(!is.na(date_lost), "Lost", ifelse(!is.na(death_date), "Decease", "OK"))),
		follow_end_day = data.table::fifelse(!is.na(Y_date), Y_date, data.table::fifelse(!is.na(death_date), death_date, as.Date("2021-12-31"))), # fifelse preserves the type and class of the inputs.
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
		follow_years_int = ceiling(follow_years),
		Z = dat[[Z]]
	) %>% filter( follow_years >0 ) 
plot(dat1$date_attend, dat1$follow_end_day)
table(dat1$Y_status, dat1$follow_years_int)
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	surv.obj[1:20]
	fit.surv <- survfit(surv.obj ~1); str(fit.surv)
	fit.surv %>% ggsurvfit::ggsurvfit() + add_confidence_interval() + add_risktable() #Kaplan-Meier curve
	ggsurvplot(survfit(surv.obj ~1, data=dat1), ylim=c(0.9,1),risk.table=TRUE)
	#palette=c("green","green4", "gold","gold4", "tomato","tomato4"), 
	ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), ylim=c(0.8,1), conf.int=FALSE, pval=FALSE, pval.method=TRUE, test.for.trend=FALSE, surv.median.line="hv", risk.table=TRUE, cumevents=FALSE)
	#dat1 %>% ggscatter(x=X, y=Y, color=Z, xlab=X, ylab=Y, add="reg.line", ellipse=TRUE, conf.int=FALSE, mean.point=TRUE)
fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+bmi+PC1+PC2, data=dat1); coef(summary(fit.cox))
	res.cox <- coef(summary(fit.cox)) #%>% as.data.frame(row.names=dimnames(.)[[1]]) %>% data.frame(X=row.names(.), row.names=NULL) 
	data.table::setDT(res.cox, keep.rownames="X")
	survminer::ggforest(fit.cox, main="", fontsize=1.2, data=dat1)
	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	res.cox as.data.frame(res.cox, row.names=dimnames(res.cox)[[1]])
	res.add <- rbind(res.cox) ## 把多个分析结果综合起来
	ggforestplot::forestplot(df=df, name=name, estimate=beta, se=se, pvalue=pvalue, psignif=0.002, colour=study, shape=study, xlab="", title="", logodds=TRUE) + 
		ggplot2::scale_shape_manual(values = c(23L, 21L, 21L, 21L, 21L), labels = c("Meta-analysis", "NFBC-1997", "DILGOM", "FINRISK-1997", "YFS"))
	
