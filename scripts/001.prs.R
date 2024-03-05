pacman::p_load(dplyr, tidyverse, ggplot2)
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)

set.seed(12345)
races=c("EUR", "SAS", "AFR", "EAS") 
dat <- readRDS("ukb.rds")
summary(lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, data=dat))

for (r in races){
	dat$PC1 <- std(dat$PC1); dat$PC2 <- std(dat$PC2) # 标准化，避免一个风头完全掩盖另外一个
	dat_sub <- subset(dat, race==r); print(paste(r,nrow(dat_sub)))
	dat_sub$prs.sameRace <- std(dat_sub[[paste0("height.",r,".score_sum")]])
	assign(paste0("PC1.center.",r), median(dat_sub$PC1))
	assign(paste0("PC2.center.",r), median(dat_sub$PC2))
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train <- lm(height ~ prs.sameRace +age+sex, data=dat_sub.train)
	lm.train.ge <- lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, data=dat_sub.train)
	y.hat <- predict(lm.train, dat_sub.valid)
	y.hat.ge <- predict(lm.train.ge, dat_sub.valid)
	print((cor(y.hat, dat_sub.valid$height, use="complete.obs"))^2)
	print((cor(y.hat.ge, dat_sub.valid$height, use="complete.obs"))^2)
}
for (r in races){ 
	eval(parse(text=paste0( 'dat$d2.', r, ' <- sqrt((dat$PC1 - PC1.center.', r, ')^2 + (dat$PC2 - PC2.center.', r, ')^2)' )))
}