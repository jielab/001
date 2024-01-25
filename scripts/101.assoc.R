setwd("D:/tmp")
pacman::p_load(data.table, readxl, lubridate, tidyverse, plotly, dplyr, survival, survminer, ggsurvfit, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读入数据，生成几个ABO相关新变量，sanity check一下 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") %>%
	rename (chunk=height_sitting) %>% 
	mutate (
	a = ifelse(abo=="A", "A", "non-A"), o = ifelse(abo=="O", "O", "non-O"),
	se = ifelse(fut2.rs601338_A==2, "non-se", "se"), # [PMID: 30345375]
	a_se = factor(paste(a, se, sep="."), levels=c("non-A.se", "non-A.non-se", "A.se", "A.non-se")), # 把最多的组作为ref
	o_se = factor(paste(o, se, sep="."), levels=c("non-O.se", "non-O.non-se", "O.se", "O.non-se")),
	s = ifelse(sp1.S==0, "non-S", "S"), z = ifelse(sp1.Z==0, "non-Z", "Z"),
	leg = height - chunk, leg_ratio = leg / chunk, chunk_ratio = chunk / leg,
	whr = waist / hip, whr.2.0 = waist.2.0 / hip.2.0, whr_diff = whr.2.0 - whr,
	bmi_diff = bmi.2.0 - bmi
	) 
summary(dat0$chunk)
	naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=TRUE)), facet=sex)
	hist(dat0$icdDate_covid, breaks="weeks", freq=TRUE)
	prop.table(table(dat0$abo, dat0$fut2.rs601338_A),1) 
	group_by(dat0, abo, se) %>% summarise(count=n(), mean=mean(bb_TES, na.rm=TRUE))
	aggregate(bb_TES ~ abo*se, dat0, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	bp <- boxplot(bb_TES ~ abo*se, dat0, las=2, col=rainbow(6), font=2); bp$stats
	dat0 %>% drop_na(abo, se) %>% ggplot(aes(x=abo, y=bb_TES, fill=se)) + geom_boxplot() + stat_summary(fun.y=mean, color="darkred", position=position_dodge(0.75), geom="point", shape=18, size=3)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=TRUE)) %>% na.omit() %>% dplyr::sample_n(10000)
	#car::scatterplotMatrix(dat, spread=FALSE, smoother.args=list(lty=0.1))
	#psych::pairs.panels(dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y 或 X-Y-Z交互作用 批量分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% drop_na(age, sex) %>% filter(ethnic_cat=="White") %>% 
	mutate( across(grep("walking", names(dat0), value=T), ~factor(.x)))
Xs <- grep("^height$|^chunk$|^leg|score_sum$", names(dat), value=TRUE) 
Ys <- grep("^icdDate", names(dat), value=TRUE)
Zs <- grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
outfile="101.assoc.tsv"; file.create(outfile)
#sink("101.assoc.log")
for (Y in Ys) {
	writeLines(paste('\n\n-->Run:', Y))
	dat1 <- dat %>%
	mutate(
		Y_date = dat[[Y]],
		Y_yes = ifelse( is.na(dat[[Y]]), 0,1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	for (X in Xs) {
		if (X==Y) next
		print(paste("RUN", X, Y))
		dat1$X = inormal(dat1[[X]])
		fit.cox <- coxph(surv.obj ~ X +age+sex +PC1+PC2, data=dat1)
		res.cox <- coef(summary(fit.cox))
		b=round(res.cox[1,1],4); se=round(res.cox[1,3],4); p=signif(res.cox[1,5],2)
		write.table(paste(Y, X, b, se, p), file=outfile, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
		dat1$X_qt <- cut(dat1$X, breaks=quantile(dat1$X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5))
		dat1$X_qt <- factor(ifelse(dat1$X_qt=="q1", "low", ifelse(dat1$X_qt=="q5", "high", "middle")), levels=c("low","middle","high"))
		next
		for (Z in Zs) {
			print(paste("RUN", X, Y, Z))
			dat1$Z <- dat1[[Z]]
			if (Z %like% "\\.rs" & length(unique(dat1$Z)) >3) dat1$Z <- hardcall(dat1$Z) # 要不然没法画出Z的3个类型
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+sex +PC1+PC2, data=dat1)
			res.cox <- coef(summary(fit.cox)); print(res.cox)
			png(paste(X, Y, Z, "png", sep="."))
			print(ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), ylim=c(0.5,1), risk.table=FALSE))	
			dev.off()
		}
	}
}
#sink()
# 可用下面代码 transpose 汇总数据，画 heatmap 图
dat <- read.table('D:/tmp/101.assoc.tsv', header=F)
dat$bnp <- paste(dat$V3, "(", dat$V5, ")", sep="")
bnp <- dat %>% reshape2::acast(V2 ~ V1, value.var='bnp')
write.table(bnp, file="101.new.tsv", sep='\t', row.names=TRUE, col.names=TRUE, append=FALSE, quote=FALSE)
#plt <- ggcorrplot(dat.b, lab=TRUE, p.mat=dat.p, sig.level=1e-4, insig ='blank') + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某一*显著*结果的精细分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X="whr_diff"
Y="icdDate_t2dm" 
dat$X = dat[[X]]
dat1 <- dat %>% drop_na(age, sex, whr_diff, bmi_diff) 
dat1 <- dat1 %>% 
	mutate(
		X = std(X),
		Z = std(y.t2d.score_sum),
		X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		Z_qt = cut(Z, breaks=quantile(Z, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
		X_qt = factor(ifelse(X_qt=="q1", "loss", ifelse(X_qt=="q5", "gain", "same")), levels=c("loss", "same", "gain")),
		Z_qt = factor(ifelse(Z_qt=="q1", "low", ifelse(Z_qt=="q5", "high", "middle")), levels=c("low", "high", "middle")),
		Y_date = dat1[[Y]],
		Y_yes = ifelse(is.na(dat1[[Y]]), 0, 1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 ) 
	aggregate(Y_yes ~ X_qt, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
	aggregate(Y_yes ~ Z_qt, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
	summary(glm(Y_yes ~ X_qt +Z_qt + X_qt*Z_qt + age+sex+bmi+whr +smoke_status+alcohol_status + PC1+PC2, data=dat1))
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	ggsurvplot(survfit(surv.obj ~X_qt + Z_qt, data=dat1), ylim=c(0.8,1), risk.table=FALSE)
	fit.cox <- coxph(surv.obj ~ X +Z + X*alcohol_status + age+sex+bmi+whr +smoke_status+alcohol_status + PC1+PC2, data=dat1)
	print(coef(summary(fit.cox)))
	survminer::ggforest(fit.cox, main="", fontsize=1.2, data=dat1) # 不能显示interaction值
	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	





