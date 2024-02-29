pacman::p_load(dplyr, tidyverse, ggplot2)
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
setwd("D:/")

set.seed(12345)
dat <- readRDS("ukb.dat")
summary(lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex+PC1+PC2, data=dat))
for (r in c("EUR", "SAS", "AFR", "EAS")){
	dat_sub <- subset(dat, race==r); print(paste(r,nrow(dat_sub)))
	dat_sub$prs.sameRace <- std(dat_sub[[paste0("height.",r,".score_sum")]])
	dat_sub$d2center = 1 / sqrt( (dat_sub$PC1 - median(dat_sub$PC1))^2 + (dat_sub$PC2 - median(dat_sub$PC2))^2 ) 
	assign(paste0("dat.",r), dat_sub) # 对于4个race, 分别生成一个新的dat
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train <- lm(height ~ prs.sameRace +age+sex+PC1+PC2, data=dat_sub.train)
	lm.train.ge <- lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex+PC1+PC2, data=dat_sub.train)
	y.hat <- predict(lm.train, dat_sub.valid); print((cor(y.hat, dat_sub.valid$height, use="complete.obs"))^2)
	y.hat.ge <- predict(lm.train.ge, dat_sub.valid); print((cor(y.hat.ge, dat_sub.valid$height, use="complete.obs"))^2)
}
dat_new <- rbind(dat.EUR, dat.SAS, dat.AFR, dat.EAS)
dat_new$tt <- sum(dat_new$d2center) #??
for (r in c("EUR", "SAS", "AFR", "EAS")){
	dat_sub <- subset(dat_new, race==r)
	for (r2 in c("EUR", "SAS", "AFR", "EAS")){
		dat_sub[[paste0("height.",r2,".score_adj")]] <- std(dat_sub[[paste0("height.",r2,".score_sum")]]) / dat_sub$tt
	}
	assign(paste0("dat.",r), dat_sub) # 对于4个race, 分别生成一个新的dat
}
dat_new <- rbind(dat.EUR, dat.SAS, dat.AFR, dat.EAS)
summary(lm(height ~ height.AFR.score_adj + height.EAS.score_adj + height.EUR.score_adj + height.SAS.score_adj + age+sex+PC1+PC2, data=dat_new))
for (r in c("EUR", "SAS", "AFR", "EAS")){
	dat_sub <- subset(dat_new, race==r)
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train.adj <- lm(height ~ height.AFR.score_adj + height.EAS.score_adj + height.EUR.score_adj + height.SAS.score_adj +age+sex+PC1+PC2, data=dat_sub.train)
	y.hat.adj <- predict(lm.train.adj, dat_sub.valid); print((cor(y.hat.adj, dat_sub.valid$height, use="complete.obs"))^2)
}