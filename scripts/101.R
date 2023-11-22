setwd("D:/101")
pacman::p_load(data.table, readxl, lubridate, tidyverse, dplyr, survival, survminer, reshape2, psych)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读入数据并生成几个大家都用的变量 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") %>% 
	mutate (
	o = ifelse(abo=="O", "O", "non-O"),
	se = ifelse(fut2.rs601338_A==2, "non-se", "se"), # 两个A是“非分泌性”
	o_se = factor(paste(o, se, sep="."), levels=c("non-O.se", "non-O.non-se", "O.se", "O.non-se")), # 把最多的组放在前面，作为ref
	s = ifelse(sp1.S==0, "non-S", "S"),
	z = ifelse(sp1.Z==0, "non-Z", "Z")
	) %>% 
	dplyr::select(-c(lct.microbe.rs3940549_A, lct.microbe.rs182549_T)) # 这两个SNP跟lct.MCM6.rs4988235_A高度相关：r>0.9


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 数据 Sanity check 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
grep("abo", names(dat0), value=T) # 找变量名字
attach(dat0) #这样，下面就可以不用写 dat0$，直接写变量名字就行
	hist(icdDate_covid, breaks="weeks", freq=TRUE); table(icd_covid)
	table(abo, abo.O1.rs8176719_T); table(abo, abo.AB.rs8176746_G) # ABO血型
	cor(abo.O1.rs8176719_T, abo.microbe.rs545971_C, use="complete.obs") # 0.93, (abo.microbe.rs550057_C) 0.80, (abo.microbe.rs8176645_T) 0.87
	prop.table(table(abo, fut2.rs601338_A),1) 
	table(sport.ACTN3.rs1815739_C); table(sport.ACE.rs4343_G); hist(sport.RBFOX1.rs7191721_G) 
	table(icd_lungcancer); hist(icdDate_lungcancer, breaks="months")
	group_by(dat, Y_yes, o_se) %>% summarise(count=n(), mean=mean(X, na.rm=T))
	aggregate(X ~ Y_yes*o_se, dat, FUN=function(x) {round(c(length(x), mean(x), sd(x), quantile(x,probs=c(0,0.5,1))), 2)} )
	ggplot(dat, aes(x=Y_yes, y=X, fill=o_se)) + geom_boxplot() + stat_summary(fun.y=mean, color="darkred", position=position_dodge(0.75), geom="point", shape=18, size=3)
	bp <- boxplot(X ~ Y_yes*se, dat, las=2, col=rainbow(6), font=2); bp$stats
detach(dat0)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=T)) %>% na.omit() %>% dplyr::sample_n(10000)
	car::scatterplotMatrix(dat, spread=F, smoother.args=list(lty=0.1))
	psych::pairs.panels(dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 北大公卫数据 Sanity check
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PKU_le8.xlsx sex +dashpts +PA_pts +smoke_pts +sleep_pts +bmi_pts +nonhdl_pts +hba1c_pts +BP_pts +LE8score +CVH_cat
#两个运动SNP：rs1815739_C, rs7191721_G
dir <- "D:/data/ukb/phe/"
pku <- read.csv(paste0(dir,"PKU_lungcancer.csv")); hist(pku$walkingpace); table(pku$lung_cancer_baseline)
	surv.obj <- Surv(time=pku$lung_cancer_time, event=pku$lung_cancer_inc)
	fit.cox <- coxph(surv.obj ~ walkingpace*rs7191721_G + walkingpace+rs7191721_G +age+sex+
	centre+Townsend_i+alcohol+Smoking+pa_3c+HDS2_i +BMI_i+grip_i+hpt_baseline+diabetes_baseline_first_report+cvd_baseline+famhis_lungcancer+genebatch+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data=pku); summary(fit.cox)  
pku <- dat0 %>% merge(pku, by.x="eid", by.y="n_eid") %>% rename(age=age.x, sex=sex.x) 
	nrow(pku); grep("\\.y", names(pku), value=T)
	plot(pku$age, pku$age.y); plot(pku$icdDate_lungcancer, pku$lung_cancer_time) #比较PKU和我们的数据
dat1 <- pku %>% # 继续比较survival分析用到的变量
	mutate(
		Y_date = icdDate_lungcancer,
		Y_yes = ifelse( is.na(icdDate_lungcancer), 0,1),
		follow_end_day = ifelse(!is.na(Y_date), Y_date, ifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	)
	table(dat1$Y_yes, dat1$lung_cancer_inc)
	plot(dat1$follow_years, dat1$lung_cancer_time) 
	surv.obj <- Surv(time=dat1$lung_cancer_time, event=dat1$lung_cancer_inc)
	fit.cox <- coxph(surv.obj ~ walkingpace+rs7191721_G + walkingpace*rs7191721_G + age+sex+genePC1+genePC2, data=dat1); summary(fit.cox)  # rs1815739_C, rs7191721_G


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# X-Y-Z绕圈圈分析示例
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- pku %>% filter(ethnic_cat=="White" & ethnicity_gen==1 & sex==1 & !is.na(o) & !is.na(se)) # ethnic_cat 见 phe.R
Xs <- "bb_TES" # grep("^bb_|walkingpace", names(dat), value=T)
Ys <- "icdDate_covid" # grep("^icdDate_", names(dat), value=T)
Zs <- "se" #grep(".rs|^rs", names(dat), value=T)
outfile="surv.res"; file.create(outfile)
for (Y in Ys) {
	print(paste("Y变量:", Y))
	dat1 <- dat %>%
	mutate(
		Y_date = dat[[Y]],
		Y_yes = ifelse( is.na(dat[[Y]]), 0,1),
		follow_end_day = ifelse(!is.na(Y_date), Y_date, ifelse(!is.na(death_date), death_date, as.Date("2022-01-01"))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(attend_date)) / 365.25
	) %>% filter( follow_years >0 )
	print(table(dat1$Y_yes))
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	for (X in Xs) {
		print(paste("X变量:", X))
		dat1$X <- dat1[[X]]
		for (Z in Zs) {
			print(paste("Z变量:", Z))
			dat1$Z <- dat1[[Z]]
			fit.glm <- glm(Y_yes ~ X +Z +X*Z +age+sex +PC1+PC2+Townsend_i, data=dat1, family="binomial")
			fit.cox <- coxph(surv.obj ~ X + Z + X*Z +age+sex +PC1+PC2+Townsend_i, data=dat1); summary(fit.cox)  
			print(res.glm <- coef(summary(fit.glm)))
			print(res.glm <- coef(summary(fit.glm)))
			#png(file=paste(X,Y,"frt.png",sep="."), w=1200, h=1600)
			#print(ggforest(fit.cox, main=paste("X:", X, "| Y:", Y), fontsize=2.2, data=dat1)); dev.off()
			#res_str <- paste(X, Y, Z, nrow(dat1), res[1], res[3], res[5], sep='|')
			#write.table(res_str, outfile, append=T, quote=F, row.names=F, col.names=F)
		}
	}
}
# 如果上面代码用了write.table，可用下面代码整合分析结果
outfile2="surv.p.tsv"; file.create(outfile2)
dat <- read.table(outfile, sep='|', header=F, as.is=T) 
	names(dat) <- c('X', 'Y', 'cnt', 'b', 'se', 'p')
	dat$p=signif(dat$p,2)
	dat.b <- dat %>% dplyr::select(X, Y, b) %>% acast(X ~ Y, value.var='b'); dat.b=round(dat.b,1)
	dat.p <- dat %>% dplyr::select(X, Y, p) %>% acast(X ~ Y, value.var='p')
	write.table(dat.p, file=outfile2, sep='\t', row.names=T, col.names=T, append=F, quote=F)
#plt <- ggcorrplot(dat.b, lab=T, p.mat=dat.p, sig.level=1e-4, insig ='blank') + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
