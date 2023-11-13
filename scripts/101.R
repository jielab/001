setwd("D:/101")
pacman::p_load(data.table, lubridate, tidyverse, dplyr, survival, reshape2)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read and check data 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata")
naniar::gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=T)), facet=sex)
grep("abo", names(dat0), value=T) # 找变量名字
attach(dat0) #这样，下面就可以不用写 dat0$，直接写变量名字就行
table(abo, abo.O1.rs8176719_T); table(abo, abo.AB.rs8176746_G) # ABO血型
	cor(abo.microbe.rs545971_C, abo.microbe.rs8176645_T, use="complete.obs") # 0.85
	cor(abo.microbe.rs545971_C, abo.microbe.rs550057_C, use="complete.obs") # 0.83
	cor(abo.microbe.rs8176645_T, abo.microbe.rs550057_C, use="complete.obs") # 0.72
prop.table(table(abo, fut2.rs601338_A),1) # wild-type G allele encodes the "secretor" (Se) allele
hist(lct.MCM6.rs4988235_A, freq=F) # In Europe, T (或者说A) allele is common and lactase persistence, C/C（或者说G/G）is lactose intolerant.
	cor(lct.MCM6.rs4988235_A, lct.microbe.rs3940549_A, use="complete.obs") # 0.94 选这个就可以了
	cor(lct.MCM6.rs4988235_A, lct.microbe.rs182549_T, use="complete.obs") # 0.99
	cor(lct.microbe.rs182549_T, lct.microbe.rs3940549_A, use="complete.obs") # 0.94
table(sp1); table(sp1.M); table(sp1.S); table(sp1.Z)
hist(age_mn) #初潮年龄
table(sport.ACTN3.rs1815739); hist(sport.RBFOX1.rs7191721) 
table(icd_lungcancer); hist(icdDate_lungcancer, breaks="months") # 癌症
table(icd_chd); hist(icdDate_chd, breaks="months") # 心血管
table(icd_covid); hist(icdDate_covid, breaks="months") # COVID-19
detach(dat0)


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
		fit.cox <- coxph(surv.obj ~ X + age+sex, data=dat1) #print(summary(fit.cox))
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
	dat.b <- dat %>% dplyr::select(X, Y, b) %>% acast(X ~ Y, value.var='b'); dat.b[is.na(dat.b)] =0; dat.b=round(dat.b,1)
	dat.p <- dat %>% dplyr::select(X, Y, p) %>% acast(X ~ Y, value.var='p'); dat.p[is.na(dat.p)] =1
	write.table(dat.p, file="p.txt", sep='\t', row.names=T, col.names=T, append=F, quote=F)
#plt <- ggcorrplot(dat.b, lab=T, p.mat=dat.p, sig.level=1e-4, insig ='blank') 
#	plt + theme(axis.text=element_text(size=12, face='bold', color=c("black","blue")))
#corrplot(dat.b, is.corr=F, method='shade', bg='black', col=colorRampPalette(c('white','green','gold'))(100), tl.col='black', tl.cex=1.3, addCoef.col='black', number.cex=0.9, insig='pch', pch.cex=2, tl.srt=45, outline=T)
