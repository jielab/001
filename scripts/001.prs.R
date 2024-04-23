pacman::p_load(dplyr, tidyverse, ggplot2, cluster, factoextra, pROC)
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
setwd("D:/")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRS-GRID
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
races=c("EUR", "SAS", "AFR", "EAS") 
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") 
dat <- dat0 %>% subset(select=grepl("^eid|^age|^sex$|^ethnic|^PC|^umap|t2dm|height", names(dat0))) %>% 
	filter(!is.na(PC1), !is.na(umap1), ethnic_cat %in% c("White","Black","Asian","Chinese")) %>%
	mutate(
		race=ifelse((PC1> -25 & PC1<25 & PC2> -25 & PC2<20), "EUR",ifelse((PC1>25 & PC1<125 & PC2> -160 & PC2< -60), "SAS", ifelse((PC1>130 & PC1<180 & PC2> -290 & PC2< -250), "EAS", ifelse((PC1>250 & PC1<425 & PC2>25 & PC2<90), "AFR", NA)))),
		#Y_yes=ifelse(is.na(icdDate_t2dm), 0, 1)
	)
	for (r in races){
		dat[[paste0('height.', r, '.score_sum')]] = std(dat[[paste0('height.', r, '.score_sum')]])
	}
ggplot(dat, aes(PC1, PC2, color=ethnic_cat)) + geom_point(size=1) + scale_color_manual(labels = c("EUR", "SAS", "AFR", "EAS"), values = c("blue", "purple", "black", "red")) + theme_bw() + guides(color=guide_legend("race")) + 
	geom_rect(aes(xmin=-25, xmax=25, ymin=-25, ymax=20), linetype=5, fill="transparent", color="blue", size=1) +
	geom_rect(aes(xmin=25, xmax=125, ymin=-160, ymax=-60), linetype=5, fill="transparent", color="purple", size=1) +
	geom_rect(aes(xmin=130, xmax=180, ymin=-290, ymax=-250), linetype=5, fill="transparent", color="red", size=1) + 
	geom_rect(aes(xmin=250, xmax=425, ymin=25, ymax=90), linetype=5, fill="transparent", color="black", size=1)
	#km = kmeans(subset(dat, select=c("PC1","PC2")), centers=4, nstart=50); fviz_cluster(km, data=dat)	
dat$PC1 <- residuals(lm(PC1 ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, na.action=na.exclude, data=dat ))
dat$PC2 <- residuals(lm(PC2 ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, na.action=na.exclude, data=dat ))
dat <- dat %>% drop_na(race, height, height.EUR.score_sum) %>% group_by(race) %>% slice_sample(n=1000) %>% dplyr::select(age, sex, race, ethnic_cat, PC1, PC2, PC3, PC4, height, height.AFR.score_sum, height.EAS.score_sum, height.EUR.score_sum, height.SAS.score_sum)
for (r in races){
	dat_sub <- subset(dat, race==r); print(paste(r,nrow(dat_sub)))
	dat_sub$prs.sameRace <- std(dat_sub[[paste0("height.",r,".score_sum")]])
	assign(paste0("PC1.center.",r), median(dat_sub$PC1))
	assign(paste0("PC2.center.",r), median(dat_sub$PC2))
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train <- lm(height ~ height.AFR.score_sum + height.EAS.score_sum + height.EUR.score_sum + height.SAS.score_sum +age+sex, data=dat_sub.train)
	y.hat <- predict(lm.train, dat_sub.valid)
	print((cor(y.hat, dat_sub.valid$height, use="complete.obs"))^2)
	# dfroc <- roc(dat_sub.valid$height, y.hat.ge)
	# print(plot(dfroc,col="red", legacy.axes=TRUE, print.auc=TRUE, print.thres=TRUE, grid=c(0.2,0.2),grid.col=c("blue","yellow")))
	# print(auc(dfroc))
	}
for (r in races){ 
	eval(parse(text=paste0( 'dat$d2.', r, ' <- sqrt((dat$PC1 - PC1.center.', r, ')^2 + (dat$PC2 - PC2.center.', r, ')^2)' )))
	eval(parse(text=paste0( 'dat$d2inv.', r, ' <- 1 / dat$d2.', r )))
	}
	dat$dtot <- dat$d2inv.AFR + dat$d2inv.EAS + dat$d2inv.EUR + dat$d2inv.SAS
for (r in races){
	eval(parse(text=paste0( 'dat$height.', r, '.score_adj <- dat$d2inv.', r, ' / dat$dtot * dat$height.', r, '.score_sum' )))
	}
for (r in races){
	print(paste("PROCESS", r))
	dat_sub <- subset(dat, race==r)
	dat_sub$prs.sameRace <- std(dat_sub[[paste0("height.",r,".score_sum")]])
	ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*0.5)))
	dat_sub.train <- dat_sub[ii,]
	dat_sub.valid <- dat_sub[-ii,]
	lm.train.adj <- lm(height ~ height.AFR.score_adj + height.EAS.score_adj + height.EUR.score_adj + height.SAS.score_adj +age+sex+PC1+PC2, data=dat_sub.train)
	y.hat.adj <- predict(lm.train.adj, dat_sub.valid)
	print((cor(y.hat.adj, dat_sub.valid$height, use="complete.obs"))^2)
	# dfroc <- roc(dat_sub.valid$height, y.hat.adj)
	# print(plot(dfroc,col="red", legacy.axes=TRUE, print.auc=TRUE, print.thres=TRUE, grid=c(0.2,0.2),grid.col=c("blue","yellow")))
	# print(auc(dfroc))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check on basic stuff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") 
dat <- dat0 %>% filter(ethnic_cat=="White") %>% drop_na(height, height.EUR.score_sum)
	Xs = c('height.EUR697', 'height.EUR3290', 'height.AFR', 'height.EAS', 'height.EUR')
	for (X in Xs) {
		dat$X <- inormal(dat[[paste0(X,'.score_sum')]])
		print(X); print(cor(dat$height, dat$X, use="complete.obs"))^2;  
	}
dat$trait <- dat$bmi; dat$prs <- dat$bmi.EUR2446.score_sum
myhist <- hist(dat$prs, breaks=20) 
	avg <- by(dat$trait, cut(dat$prs, breaks=myhist$breaks), function(x) mean(x,na.rm=T))
	par(new=T); plot(myhist$mids, avg, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
	axis(side=4); mtext(side=4, line=3, 'measured')	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRS对比结果展示图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:/Users/jiehu/Desktop")
pacman::p_load(readxl, tidyverse, ggplot2)
dat <- read_excel("PRS.xlsx",sheet="Sheet2")
barplot(t(as.matrix(dat[, 2:3])), 
        beside = TRUE,
        names.arg = dat$人种,
		main="两种方法对height的预测",
        legend.text = TRUE,
		ylim = c(0,0.15),
		col = c("green", "blue"),
        ylab = expression(italic("Nagelkerke R2")),
        xlab = "人种")


