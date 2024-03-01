pacman::p_load(dplyr, tidyverse, ggplot2)
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
setwd("D:/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PCA图及人种划分
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
races=c("EUR", "SAS", "AFR", "EAS") 
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") 
dat <- dat0 %>% subset(select=grepl("^eid|^age|^sex$|^ethnic|^garray|^batch|^PC|height$|height.AFR.score_sum|height.EUR.score_sum|height.EAS.score_sum|height.HIS.score_sum|height.SAS.score_sum", names(dat0))) %>% 
	filter(!is.na(PC1), !is.na(ethnic_cat), ethnic_cat !="Other", ethnic_cat !="Mixed") %>%
	mutate(race=ifelse((PC1> -25 & PC1<25 & PC2> -25 & PC2<20), "EUR",ifelse((PC1>25 & PC1<125 & PC2> -160 & PC2< -60), "SAS", ifelse((PC1>130 & PC1<180 & PC2> -290 & PC2< -250), "EAS", ifelse((PC1>250 & PC1<425 & PC2>25 & PC2<90), "AFR", NA)))))
dat1 <- dat %>% group_by(ethnic_cat) %>% slice_sample(n=1000)
ggplot(dat1, aes(PC1, PC2, color=ethnic_cat)) + geom_point(size=1) + scale_color_manual(labels = c("白人", "南亚人", "黑人", "中国人"), values = c("blue", "purple", "black", "red")) + theme_bw() + guides(color=guide_legend("人种")) + 
	geom_rect(aes(xmin=-25, xmax=25, ymin=-25, ymax=20), linetype=5, fill="transparent", color="blue", size=1) +
	geom_rect(aes(xmin=25, xmax=125, ymin=-160, ymax=-60), linetype=5, fill="transparent", color="purple", size=1) +
	geom_rect(aes(xmin=130, xmax=180, ymin=-290, ymax=-250), linetype=5, fill="transparent", color="red", size=1) + 
	geom_rect(aes(xmin=250, xmax=425, ymin=25, ymax=90), linetype=5, fill="transparent", color="black", size=1)
dat <- dat %>% drop_na(race)
summary(lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, data=dat))
for (r in races){
	dat_sub <- subset(dat, race==r); print(paste(r,nrow(dat_sub)))
	dat_sub$prs.sameRace <- std(dat_sub[[paste0("height.",r,".score_sum")]])
	assign(paste0("PC1.center.",r), median(dat_sub$PC1))
	assign(paste0("PC2.center.",r), median(dat_sub$PC2))
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train <- lm(height ~ prs.sameRace +age+sex, data=dat_sub.train)
	lm.train.ge <- lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, data=dat_sub.train)
	y.hat <- predict(lm.train, dat_sub.valid); print((cor(y.hat, dat_sub.valid$height, use="complete.obs"))^2)
	y.hat.ge <- predict(lm.train.ge, dat_sub.valid); print((cor(y.hat.ge, dat_sub.valid$height, use="complete.obs"))^2)
	}
for (r in races){
	eval(parse(text=paste0( 'dat$d2.', r, ' <- 1 / sqrt((dat$PC1 - PC1.center.', r, ')^2 + (dat$PC2 - PC2.center.', r, ')^2)' )))
	}
	dat$dtot <- dat$d2.AFR + dat$d2.EAS + dat$d2.EUR + dat$d2.SAS
for (r in races){
	eval(parse(text=paste0( 'dat$height.', r, '.score_adj <- dat$d2.', r, ' / dat$dtot * dat$height.', r, '.score_sum' )))
	}
summary(lm(height ~ height.AFR.score_adj + height.EAS.score_adj + height.EUR.score_adj + height.SAS.score_adj + age+sex, data=dat))
for (r in races){
	dat_sub <- subset(dat, race==r)
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train.adj <- lm(height ~ height.AFR.score_adj + height.EAS.score_adj + height.EUR.score_adj + height.SAS.score_adj +age+sex+PC1+PC2, data=dat_sub.train)
	y.hat.adj <- predict(lm.train.adj, dat_sub.valid); print((cor(y.hat.adj, dat_sub.valid$height, use="complete.obs"))^2)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check on basic stuff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=="White") %>% drop_na(height, height.EUR.score_sum)
	Xs = c('height.EUR697', 'height.EUR3290', 'height.AFR', 'height.EAS', 'height.EUR', 'height.ALL')
	for (X in Xs) {
		dat$X <- inormal(dat[[paste0(X,'.score_sum')]])
		print(X); print(cor(dat$height, dat$X, use="complete.obs"))^2;  
	}
dat$trait <- dat$bmi; dat$prs <- dat$bmi.EUR2446.score_sum
myhist <- hist(dat$prs, breaks=100) 
	avg <- by(dat$trait, cut(dat$prs, breaks=myhist$breaks), function(x) mean(x,na.rm=T))
	par(new=T); plot(myhist$mids, avg, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
	axis(side=4); mtext(side=4, line=3, 'measured')	
dat <- dat0 %>% filter(ethnic_cat=="White") %>% 
	rename(bmi.prs=bmi.EUR2446.score_sum, height.prs=height.EUR.score_sum, cad.prs=cad.mult.score_sum) %>% 
	dplyr::select(eid, age, sex, bmi, bmi.prs, height.prs, cad.prs, date_attend, icdDate_cad) 
	saveRDS(dat, "ukb.rds")

