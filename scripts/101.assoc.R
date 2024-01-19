setwd("D:/tmp")
pacman::p_load(data.table, readxl, lubridate, tidyverse, plotly, dplyr, survival, survminer, ggsurvfit, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读入数据，生成几个ABO相关新变量，sanity check一下 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") %>%
	dplyr::select(-c(lct.microbe.rs3940549_A, lct.microbe.rs182549_T)) %>% # 这两个SNP跟lct.MCM6.rs4988235_A高度相关：r>0.9
	mutate (
	a = ifelse(abo=="A", "A", "non-A"),
	o = ifelse(abo=="O", "O", "non-O"),
	se = ifelse(fut2.rs601338_A==2, "non-se", "se"), # 两个A是“非分泌性”
	a_se = factor(paste(a, se, sep="."), levels=c("non-A.se", "non-A.non-se", "A.se", "A.non-se")), # 把最多的组放在前面，作为ref
	o_se = factor(paste(o, se, sep="."), levels=c("non-O.se", "non-O.non-se", "O.se", "O.non-se")), # 把最多的组放在前面，作为ref
	s = ifelse(sp1.S==0, "non-S", "S"),
	z = ifelse(sp1.Z==0, "non-Z", "Z"),
	leg = height - height_sitting,
	leg_ratio = leg / height_sitting
	) %>% rename (chunk=height_sitting)
	summary(dat0$chunk)
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
# X-Y 或 X-Y-Z交互作用 批量分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% drop_na(age, sex) %>%
	filter(ethnic_cat=="White") # %>% mutate( across(grep("score_sum", names(dat0), value=T), ~std(.x))) # 批量变成 ~factor(.x)
	cor(dat$bmi, dat$bmi.EUR77.score_sum, use="complete.obs") # bmi.EUR941, bmi.EUR2446, height.EUR697, height.EUR3290
	coef(summary(lm(bmi ~ bmi.EUR77.score_sum, data=dat)))
Xs <- grep("^height$|^leg$|^chunk$", names(dat), value=TRUE) # |score_sum$
Ys <- grep("^icdDate", names(dat), value=TRUE)
Zs <- grep("^o$|^se$", names(dat), value=TRUE) # |^rh|shbg|^apoe$|\\.rs
outfile="101.assoc.tsv"; file.create(outfile)
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
		writeLines(paste('\n\n-->Run:', X, Y))
		dat1$X = inormal(dat1[[X]])
		fit.cox <- coxph(surv.obj ~ X +age+sex +PC1+PC2, data=dat1)
		res.cox <- coef(summary(fit.cox)); print(res.cox[1,])
		b=round(res.cox[1,1],4); se=round(res.cox[1,3],4); p=signif(res.cox[1,5],2)
		write.table(paste(Y, X, b, se, p), file=outfile, sep='\t', row.names=FALSE, col.names=FALSE, append=TRUE, quote=FALSE)
		dat1$X_qt <- cut(dat1$X, breaks=quantile(dat1$X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5))
		dat1$X_qt <- factor(ifelse(dat1$X_qt=="q1", "low", ifelse(dat1$X_qt=="q5", "high", "middle")), levels=c("low","middle","high"))
		png(paste(X,Y,"png",sep="."))
		print(ggsurvplot(survfit(surv.obj ~ X_qt, data=dat1), ylim=c(0.5,1), risk.table=FALSE))	
		dev.off()
		next # 有了这个，loop 后面的代码就代码就不运行了
		for (Z in Zs) {
			dat1$Z <- dat1[[Z]]
			if (Z %like% "\\.rs" & length(unique(dat1$Z)) >3) dat1$Z <- hardcall(dat1$Z) # 要不然没法画出Z的3个类型
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+bmi +PC1+PC2, data=dat1)
			res.cox <- coef(summary(fit.cox)); print(res.cox[1,])
			p_int <- signif(tail(res.cox,1)[,5] ,2)
			if (!is.na(p_int) & p_int < 0.05) {
				print(paste(X, Y, Z, p_int)) 
				png(paste(X,Y,Z,"png",sep="."))
				print(ggsurvplot(survfit(surv.obj ~X_qt + Z, data=dat1), ylim=c(0.8,1), risk.table=FALSE))	
				dev.off()
			}
		}
	}
}
# 可用下面代码 transpose 汇总数据，画 heatmap 图
dat <- read.table('D:/analysis/mr/pheno/res.p.txt', header=F); dat$V3 <- signif(dat$V3,2)
pval <- dat %>% reshape2::acast(V1 ~ V2, value.var='V3'); pval=signif(pval,2) 
write.table(pval, file=outfile, sep='\t', row.names=TRUE, col.names=TRUE, append=FALSE, quote=FALSE)
#plt <- ggcorrplot(dat.b, lab=TRUE, p.mat=dat.p, sig.level=1e-4, insig ='blank') + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 针对某一*显著*结果的精细分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X="walking_pace" # "bb_SHBG"
Y="icdDate_lungcancer" # "icdDate_t2dm"
dat <- dat0 %>% filter(ethnic_cat=="White") # 剩下 742,590 人
summary(dat[[Y]])
dat$walking_pace <- factor(dat$walking_pace, levels=c(2,1,3), labels=c("steady", "slow", "fast"))
dat1 <- dat %>%  
	mutate(
		X = dat[[X]],
	#	X_qt = cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
	#	X_qt = factor(ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")), levels=c("low","middle","high")),
		Y_date = dat[[Y]],
		Y_yes = ifelse(is.na(dat[[Y]]), 0, 1),
		follow_end_day = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(death_date), death_date, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 ) # 最后这一步，过滤掉了在date_attend之前就有该疾病的人（N=314）
	table(dat1$Y_yes, floor(dat1$follow_years))
	aggregate(Y_yes ~ walking_pace, dat1, FUN=function(x) { paste( length(x), sum(x), round(sum(x)/length(x),3)) } )
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
for (Z in c("power.ACTN3.rs1815739_C", "endurance.RBFOX1.rs7191721_G", "cvd.9p21.rs2383206_G", "bb.TP53.rs1042522_C")) {
	dat1$Z <- hardcall(dat1[[Z]])
	png(paste(X,Y,Z,"png",sep="."))
	print(ggsurvplot(survfit(surv.obj ~ X + Z, data=dat1), ylim=c(0.95,1), conf.int=FALSE, risk.table=FALSE)) # palette=c("green","green4", "gold","gold4", "tomato","tomato4"), 
	dev.off()
}
fit.cox <- coxph(surv.obj ~ X + Z +X*Z +age+bmi+PC1+PC2, data=dat1); coef(summary(fit.cox))
	survminer::ggforest(fit.cox, main="", fontsize=1.2, data=dat1) # 不能显示interaction值
	fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	# ggforestplot::forestplot
