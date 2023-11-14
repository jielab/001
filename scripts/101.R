setwd("D:/101")
pacman::p_load(data.table, readxl, lubridate, tidyverse, dplyr, survival, reshape2)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read and check data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
grep("abo", names(dat0), value=T) # 找变量名字
attach(dat0) #这样，下面就可以不用写 dat0$，直接写变量名字就行
table(abo, abo.O1.rs8176719_T); table(abo, abo.AB.rs8176746_G) # ABO血型
cor(abo.O1.rs8176719_T, abo.microbe.rs545971_C, use="complete.obs") # 0.93, (abo.microbe.rs550057_C) 0.80, (abo.microbe.rs8176645_T) 0.87
prop.table(table(abo, fut2.rs601338_A),1) # wild-type G allele encodes the "secretor" (Se) allele
hist(lct.MCM6.rs4988235_A, freq=F) # In Europe, T (或者说A) allele is common and lactase persistence, C/C（或者说G/G）is lactose intolerant.
cor(lct.MCM6.rs4988235_A, lct.microbe.rs3940549_A, use="complete.obs") # 0.94, (lct.microbe.rs182549_T) 0.99, (lct.microbe.rs3940549_A) 0.94
table(sp1); table(sp1.M); table(sp1.S); table(sp1.Z)
hist(age_mn) #初潮年龄
table(sport.ACTN3.rs1815739_C); table(sport.ACE.rs4343_G); hist(sport.RBFOX1.rs7191721_G) 
table(icd_lungcancer); hist(icdDate_lungcancer, breaks="months") # 癌症
table(icd_chd); hist(icdDate_chd, breaks="months") # 心血管
table(icd_covid); hist(icdDate_covid, breaks="months") # COVID-19
detach(dat0)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 读取北医的数据 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#PKU_le8.xlsx sex +dashpts +PA_pts +smoke_pts +sleep_pts +bmi_pts +nonhdl_pts +hba1c_pts +BP_pts +LE8score +CVH_cat
dir <- "D:/data/ukb/phe/"
dat <- read.csv(paste0(dir,"PKU_lungcancer.csv"))
	hist(dat$walkingpace); table(dat$lung_cancer_baseline)
surv.obj <- Surv(time=dat$lung_cancer_time, event=dat$lung_cancer_inc)
fit.cox <- coxph(surv.obj ~ walkingpace*rs7191721_G + walkingpace+rs7191721_G +age+sex+centre+Townsend_i+alcohol+Smoking+pa_3c+HDS2_i +BMI_i+grip_i+hpt_baseline+diabetes_baseline_first_report+cvd_baseline+famhis_lungcancer+genebatch+genePC1+genePC2+genePC3+genePC4+genePC5+genePC6+genePC7+genePC8+genePC9+genePC10, data=dat)
summary(fit.cox)  # rs1815739_C, rs7191721_G


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 万有引力之关联分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnicity_gen==1) # 欧洲白人
Xs <- grep(".rs", names(dat), value=T)
Ys <- grep("^icdDate|^fod_", names(dat), value=T)
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
		fit.cox <- coxph(surv.obj ~ X + age+sex+PC1+PC2, data=dat1) #print(summary(fit.cox))
		#png(file=paste(X,Y,"frt.png",sep="."), w=1200, h=1600)
		#print(ggforest(fit.cox, main=paste("X:", X, "| Y:", Y), fontsize=2.2, data=dat1)); dev.off()
		res <- coef(summary(fit.cox))[1,]
		res_str <- paste(X, Y, nrow(dat1), res[1], res[3], res[5], sep='|')
		write.table(res_str, "surv.res", append=T, quote=F, row.names=F, col.names=F)
	}
}
# 整合分析结果
dat <- read.table('D:/101/surv.res', sep='|', header=F, as.is=T) 
	names(dat) <- c('X', 'Y', 'cnt', 'b', 'se', 'p')
	dat$p=signif(dat$p,2)
	dat.b <- dat %>% dplyr::select(X, Y, b) %>% acast(X ~ Y, value.var='b'); dat.b=round(dat.b,1)
	dat.p <- dat %>% dplyr::select(X, Y, p) %>% acast(X ~ Y, value.var='p')
	write.table(dat.p, file="surv.p.tsv", sep='\t', row.names=T, col.names=T, append=F, quote=F)
#plt <- ggcorrplot(dat.b, lab=T, p.mat=dat.p, sig.level=1e-4, insig ='blank') 
#	plt + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
#corrplot(dat.b, is.corr=F, method='shade', bg='black', col=colorRampPalette(c('white','green','gold'))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig='pch', pch.cex=2, tl.srt=45, outline=T)
