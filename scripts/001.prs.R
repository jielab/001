setwd("D:/")
pacman::p_load(data.table, dplyr, tidyr, ggplot2, survival, pROC)
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
inormal <- function(x) qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PRS-GRID
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(12345)
trait='t2dm'; model='cox' # lm 或 glm 或 cox
races=c("EUR", "AFR", "EAS", "SAS") 
train_n=0.5
pcs=4 # 40个 PC
folds=10 # 10-fold cross validation

dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") 
dat <- dat0 %>% subset(select=grepl("^eid|^age|^sex$|^ethnic|^PC|^umap|t2dm|height|date_attend|date_lost|date_death", names(dat0))) %>% 
	drop_na(PC1, umap1) %>% filter(ethnic_cat %in% c("White", "Black", "Chinese", "Asian")) %>%
	mutate(race=ifelse((PC1> -25 & PC1<25 & PC2> -25 & PC2<20), "EUR",ifelse((PC1>25 & PC1<125 & PC2> -160 & PC2< -60), "SAS", ifelse((PC1>130 & PC1<180 & PC2> -290 & PC2< -250), "EAS", ifelse((PC1>250 & PC1<425 & PC2>25 & PC2<90), "AFR", NA)))))
# ggplot(dat, aes(PC1, PC2, color=ethnic_cat)) + geom_point(size=1) + scale_color_manual(labels = c("EUR", "SAS", "AFR", "EAS"), values = c("blue", "purple", "black", "red")) + theme_bw() + guides(color=guide_legend("race")) + 
	# geom_rect(aes(xmin=-25, xmax=25, ymin=-25, ymax=20), linetype=5, fill="transparent", color="blue", size=1) +geom_rect(aes(xmin=25, xmax=125, ymin=-160, ymax=-60), linetype=5, fill="transparent", color="purple", size=1) +
	# geom_rect(aes(xmin=130, xmax=180, ymin=-290, ymax=-250), linetype=5, fill="transparent", color="red", size=1) +geom_rect(aes(xmin=250, xmax=425, ymin=25, ymax=90), linetype=5, fill="transparent", color="black", size=1)
	# km = kmeans(subset(dat, select=c("PC1","PC2")), centers=4, nstart=50); fviz_cluster(km, data=dat)	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 处理PC & 计算d
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (r in races){
	dat[[paste0(r, '.prs')]] = std(dat[[paste0(trait, '.', r, '.score_sum')]]) # 可尝试去掉std()
}
for (i in 1:pcs) {
  PCi <- paste0("PC", i)
  dat[[PCi]] <- std(residuals(lm(as.formula(paste0(PCi, " ~ AFR.prs + EAS.prs + EUR.prs + SAS.prs + age + sex")), na.action=na.exclude, data=dat)))
}
# dat <- dat %>% drop_na(race, trait, EUR.prs) %>% group_by(race) %>% slice_sample(n=1000) 
for (r in races){
	dat_sub <- subset(dat, race==r)
	for (i in 1:pcs) {
		assign(paste0("PC", i, ".center.", r), median(dat_sub[[paste0("PC", i)]],na.rm=TRUE)) # 🐕 必须加上na.rm
	}
}
for (r in races){
	dat$d1 <- 0
	for (i in 1:pcs) {
		eval(parse(text = paste0('dat$d1 <- dat$d1 + (dat$PC',i, ' - PC',i,'.center.',r,')^2'))) 
	}
	eval(parse(text=paste0( 'dat$d2.', r, ' <- sqrt(dat$d1)' ))) 
	eval(parse(text=paste0( 'dat$d2inv.', r, ' <- 1 / dat$d2.', r )))
}
dat$dtot <- dat$d2inv.AFR + dat$d2inv.EAS + dat$d2inv.EUR + dat$d2inv.SAS
for (r in races){
	eval(parse(text=paste0( 'dat$', r, '.prs_adj <- dat$d2inv.', r, ' / dat$dtot * dat$', r, '.prs' )))
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 进行拟合 regression
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
if (model=="lm") {dat$trait=dat[[trait]]
} else if (model=="glm") {dat$trait=ifelse(is.na(dat[[paste0('icdDate_',trait)]]), 0, 1)
} else if (model=="cox") {
	dat <- dat %>% mutate(
		Y_date = dat[[paste0('icdDate_',trait)]], Y_yes = ifelse( is.na(Y_date), 0,1),
		follow_end_day = data.table::fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
		follow_years = (as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	) %>% filter( follow_years >0 )
}
for (r in races){
	print(paste("PROCESS", r))
	fit_C = fit_G = list()
	fit_folds_C = fit_folds_G = numeric(folds)
	dat_sub <- subset(dat, race==r)
	for (i in 1:folds) { 
		ii <- sort(sample(1:nrow(dat_sub), round(nrow(dat_sub)*train_n)) )
		dat_sub.train <- dat_sub[ii,]
		dat_sub.valid <- dat_sub[-ii,]
		if (model=="lm") { 
			model_C.train <- lm(trait ~ AFR.prs + EAS.prs + EUR.prs + SAS.prs +age+sex, data=dat_sub.train)
			model_G.train <- lm(trait ~ AFR.prs_adj + EAS.prs_adj + EUR.prs_adj + SAS.prs_adj +age+sex, data=dat_sub.train)
			y_C.hat <- predict(model_C.train, dat_sub.valid)
			y_G.hat <- predict(model_G.train, dat_sub.valid)
			fit_folds_C[i] <- (cor(y_C.hat, dat_sub.valid$trait, use="complete.obs"))^2
			fit_folds_G[i] <- (cor(y_G.hat, dat_sub.valid$trait, use="complete.obs"))^2
		} else if (model=="glm") {
			model_C.train <- glm(trait ~ AFR.prs + EAS.prs + EUR.prs + SAS.prs +age+sex, data=dat_sub.train)
			model_G.train <- glm(trait ~ AFR.prs_adj + EAS.prs_adj + EUR.prs_adj + SAS.prs_adj +age+sex, data=dat_sub.train)
			y_C.hat <- predict(model_C.train, dat_sub.valid)
			y_G.hat <- predict(model_G.train, dat_sub.valid)
			fit_folds_C[i] <- auc(roc(dat_sub.valid$trait, y_C.hat))
			fit_folds_G[i] <- auc(roc(dat_sub.valid$trait, y_G.hat))
		} else if (model=="cox") {
			surv.obj <- Surv(time=dat$follow_years, event=dat$Y_yes)
			train_cox(surv.obj ~ AFR.prs + EAS.prs + EUR.prs + SAS.prs +age+sex, data=dat_sub.train, newdata=JIE, trControl=NULL, parallel=FALSE, mc.cores=2, seed=12345)
			xxx <- get.risk.coxph(xxx, 365.25*3)
		}
	}
	print(fit_C[[r]] <- fit_folds_C)
	print(fit_G[[r]] <- fit_folds_G)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 画图
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
fit_C <- as.data.frame(fit_C); fit_C$Method <- "PRS-CSx"
fit_G <- as.data.frame(fit_G); fit_G$Method <- "PRS-GRID"
combined_data <- rbind(fit_C, fit_G)
long_data <- pivot_longer(combined_data,  cols = c("EUR", "SAS", "AFR", "EAS"),   names_to = "Race", values_to = "value")
summary_stats <- long_data %>% group_by(Race, Method) %>% summarize(mean = mean(value), sd = sd(value))
F1 <- ggplot(summary_stats, aes(x = fct_relevel(Race, "EUR", "AFR", "SAS", "EAS"), y = mean, fill = fct_relevel(Method, "PRS-CSx", "PRS-GRID",))) +
	geom_bar(stat = "identity", position = position_dodge(), linewidth = 0.8, color = "white", alpha = 0.8) +
	geom_text(aes(label=round(summary_stats$mean,3)),size=3.5,vjust=7,position =position_dodge(0.9))+
	geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(0.9), linewidth = 0.8, color = "grey50") +
	labs(title = "Mean AUC with SD given by 10-fold cross validation", x = "Race", y = "Mean Value", fill = "Method") +
	theme_bw() +
	theme(plot.tag = element_text(size = 15), axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), axis.title.y= element_text(size = 15), strip.text = element_text(size = 12))+scale_fill_manual(values = c("grey70","#4DBBD5FF"))
tiff(file = "Figure2.tiff", compression="lzw",width = 6000, height = 4000, res = 600)
dev.off()
#dat <- read_excel("PRS.xlsx",sheet="Sheet2")
#barplot(t(as.matrix(dat[, 2:3])), beside = TRUE, names.arg = dat$人种, main="两种方法对height的预测", legend.text = TRUE,
#	ylim = c(0,0.15), col = c("green", "blue"), ylab = expression(italic("Nagelkerke R2")), xlab = "人种")
	

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sanity check on basic stuff
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS("D:/data/ukb/Rdata/all.Rdata") 
dat <- dat0 %>% filter(ethnic_cat=="White") %>% drop_na(height, height.EUR.prs)
	Xs = c('height.EUR697', 'height.EUR3290', 'height.AFR', 'height.EAS', 'height.EUR')
	for (X in Xs) {
		dat$X <- inormal(dat[[paste0(X,'.prs')]])
		print(X); print(cor(dat$height, dat$X, use="complete.obs"))^2;  
	}
dat$trait <- dat$bmi; dat$prs <- dat$bmi.EUR2446.prs
myhist <- hist(dat$prs, breaks=20) 
	avg <- by(dat$trait, cut(dat$prs, breaks=myhist$breaks), function(x) mean(x,na.rm=T))
	par(new=T); plot(myhist$mids, avg, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2)
	axis(side=4); mtext(side=4, line=3, 'measured')	