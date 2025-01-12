pacman::p_load(data.table, tidyverse, lubridate, survival, vivid, rcssci)

dir0='D:/'
source(paste0(dir0, '/scripts/f/f.R'))

dat0 <- readRDS(file=paste0(dir0, "/data/ukb/phe/Rdata/all.plus.rds")); sum(!is.na(dat0$date_death))
dat <- dat0 %>% filter(ethnic_cat=="White", sex==1) %>% mutate(
	bald12=ifelse(bald==1,0, ifelse(bald==2,1,NA)),
	bald13=ifelse(bald==1,0, ifelse(bald==3,1,NA)),
	bald14=ifelse(bald==1,0, ifelse(bald==4,1,NA)),
	sexf_cat=factor(ifelse(age_fsex<16, "early", ifelse(age_fsex<22, "average", "late")),levels=c("early", "average", "late")),
	sexp_cat=factor(ifelse(sex_partner<=1, "low", ifelse(sex_partner<=10, "average", "high")),levels=c("low", "average", "high"))
)
table(dat$sexf_cat); table(dat$sexp_cat)

covs <- "age bmi smoke_status alcohol_status deprivation sexf_cat sexp_cat bb_TES bb_SHBG bb_VITD bb_IGF1" %>% strsplit(" ") %>% unlist()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table 1: 基本信息及Pheno关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(compareGroups); descrTable(formula=as.formula(paste(c("bald ~ ", paste(covs, collapse="+")), collapse="")), data=dat, digits=1)
library(gtsummary); dat %>% filter(sex==1) %>% select(all_of(c("bald", covs))) %>% tbl_summary(by=bald, missing="no") %>% add_p()
Y="bald12"; fit <- glm(formula=as.formula(paste(c(Y, " ~ ", paste(covs, collapse="+")), collapse="")), data=dat, family="binomial"); summary(fit)
library(forestmodel); forest_model(fit)
# fit <- lm(age_fsex ~ factor(bald) +..., data=dat); 
# fit <- nnet::multinom(bald ~ covs, data=dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome🥚的PCA和cluster
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)
prots <- grep("prot_", names(dat), value=T)
dat <- dat[!apply(dat[prots], 1, function(x) all(is.na(x))), ] # 去掉没有任何蛋白质数据的人
dat1 <- dat %>% select(all_of(prots)) # 只保留protein的所有列
	print(mis.top <- dat1 %>% summarise(across(everything(), ~ mean(is.na(.)))) %>% gather(variable, missing_rate) %>% arrange(desc(missing_rate)) %>% head(., 20))
	to_exclude <- mis.top[mis.top$missing_rate>0.25, "variable"]
	dat1 <- dat1 %>% select(-all_of(to_exclude)) # 去掉缺失数据太多的protein
	prots=setdiff(prots, to_exclude) # 更新 protein 名单
	names(dat1) <- gsub("prot_", "", names(dat1))
dat1 <- as.data.frame(lapply(dat1, function(x) { x[is.na(x)] <- mean(x, na.rm=TRUE); return(x)})) # 填补protein缺失数据
	dat1 <- t(dat1); pca_result <- prcomp(dat1, scale.=TRUE) # t()之后就能计算🥚PCA，否则计算🙍的PCA
	screeplot(pca_result, main="Scree Plot of PCA", col="blue", type="lines")
	biplot(pca_result, main="PCA Biplot")
	kmeans_result <- kmeans(dat1, centers=3)
	dat1$cluster <- factor(kmeans_result$cluster)
	dat1.pca <- data.frame(pca_result$x)
	dat1.pca$cluster <- factor(kmeans_result$cluster)
	ggplot(dat1.pca, aes(x=PC1, y=PC2, color=cluster)) + geom_point(alpha=0.6) +labs(title="", x="PC1", y="PC2") + theme_minimal()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome🥚的关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y="bald12"
prots <- grep("prot_", names(dat), value=T) # 🏮 prot_ met_ bb_ bc_
res <- data.frame(prot=character(), BETA=numeric(), P=numeric(), stringsAsFactors=FALSE)
for (X in prots) {
	if (sum(!is.na(dat[[X]])) < 10000) {
		res <- rbind(res, data.frame(prot=X, BETA=NA, P=NA))
	} else {
		model <- glm(as.formula(paste(Y, "~", X, "+", paste(covs, collapse=" + "))), data=dat, family="binomial")
		res <- rbind(res, data.frame(prot=X, BETA=summary(model)$coef[2,1], P=summary(model)$coef[2,4]))
	}
}
valcano(res, "prot", "BETA", "P", 0.05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome🥚的 cisMr 🏮
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Y="bald14"
res <- read.table("D:/analysis/assoc.sum/bald.cisMr.log", header=TRUE) %>% filter(Y== !!Y) # 🍬
valcano(res, "X", "BETA", "P", 0.05)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 一般的 survival 🏃‍分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ys <- 'death'
X <- "bb_LDL"
dat1 <- dat %>% mutate(
	Y_date=dat[[paste0('icdDate_',Y)]],
	Y_yes=ifelse( is.na(Y_date), 0,1),
	follow_end_day=fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
	follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25,
	X=dat[[X]],
	X_qt=cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
	X_qt=ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")),
) %>% filter(follow_years>0)
	median(dat1$follow_years)); sum(dat1$follow_years) # person-years
	round(prop.table(table(dat1$alcohol_status, dat1$Y_yes), 1), 3)
rcssci_cox(data=dat1, time="follow_years", y="Y_yes", x="X", covs=c("age","sex"), prob=0.1, filepath= 'D:/tmp') # 🐂
surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	fit.cox <- coxph(surv.obj ~ X + age+sex+bmi+smoke_status+alcohol_status +PC1+PC2+PC3+PC4, data=dat1); print(coef(summary(fit.cox)))	
	# fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	km.obj <- survfit(surv.obj ~ X_qt, data=dat1) # KM 是不能带协变量的，M会被作为分层变量
	# survdiff(surv.obj ~ X_qt, data=dat1) # log-rank test
	# plot(km.obj, fun=function(x) 1-x)
	survminer::ggsurvplot(km.obj, ylim=c(0,0.2), fun="event", pval=TRUE, risk.table=TRUE, ncensor.plot=TRUE, ggtheme=theme_bw(), palette=c("green","gray","orange"))
	survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1)
	library(forplo) 🏮
		exdf <- data.frame(cbind(OR=c(1.21,0.90,1.02, 1.54,1.32,0.79,1.38,0.85,1.11, 1.58,1.80,2.27),LCI=c(0.82,0.61,0.66,1.08,0.91,0.48,1.15,0.39,0.91,0.99,1.48,0.92),UCI=c(1.79,1.34,1.57, 2.19,1.92,1.32,1.64,1.87,1.34, 2.54,2.19,5.59), groups=c(1,1,1, 2,2,2,2,2,2, 3,3,3)))
		rownames(exdf) <- c('Barry, 2005', 'Frances, 2000', 'Rowley, 1995', 'Biro, 2000', 'Crowe, 2010', 'Harvey, 1996', 'Johns, 2004', 'Parr, 2002', 'Zhang, 2011', 'Flint, 1989', 'Mac Vicar, 1993', 'Turnbull, 1996')
		knitr::kable(exdf)
		forplo(exdf[,1:3], groups=exdf$groups, grouplabs=c('Low risk of bias', 'Some concerns', 'High risk of bias'),
			left.align=TRUE, col=2, ci.edge=FALSE, diamond=c(4,11,15), diamond.col='#b51b35')
## 10年风险计算
	#	XnM=factor(paste(M, X, sep="|"), levels=paste(rep(levels(M),each=3), rep(levels(X),3), sep="|")),
	fit.cox <- coxph(surv.obj ~ X+M, data=dat1) 
	dat1$risk <- nricens::get.risk.coxph(mdl, 10) # 每人的risk 🔔
	dat2 <- expand.grid(X=levels(dat1$X), M=levels(dat1$M))
	surv.pred <- survfit(fit.cox, newdata=dat2)
	surv.10y <- summary(surv.pred, times=10)
	dat2 <- dat2 %>% mutate(risk=1-t(surv.10y$surv), ci_lower=t(1-surv.10y$upper), ci_upper=t(1-surv.10y$lower))
	ggplot(dat2, aes(x=M, y =risk, fill=X)) + geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7) +
		geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
		geom_text(aes(label=sprintf("%.2f", risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
		labs(x="X labels", y="10-year risk (%)", title="") + scale_fill_manual(values=c("green", "gray", "orange"), name="Legend name") + theme_minimal() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome与常见病的热力图 🌋
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ggpairs(tips, mapping=aes(color=sex), columns=c("total_bill", "time", "tip"))
library(ComplexHeatmap); library(GGally) # BiocManager::install(version="3.20")
	Ys <- c("asthma", "breast_cancer", "cad", "cancer", "circulatory", "copd", "covid", "lung_cancer", "metabolic", "mi", "pancreatic_cancer", "pe", "prostate_cancer", "ra", "stroke", "t1dm", "t2dm", "varicose", "vte")
	icdDate_Ys <- paste0("icdDate_", Ys)
	time2_Ys <- paste0("time2_", Ys)
	Xs <- grep("bb_", names(dat), value=TRUE) 
	dat1 <- subset(dat, select=Xs[1:10])[1:10000,]
		ggpairs(dat1)+theme_bw()
		data(tips, package="reshape")
		ggpairs(tips[, c(1, 3, 4, 2)], upper=list(continuous="density", combo="box_no_facet"), lower=list(continuous="points", combo="dot_no_facet"))
	dat1 <- dat %>% mutate(across(.cols=intersect(names(dat)[matches("icdDate_", ignore.case=TRUE)], icdDate_Ys), .fns=~ as.numeric(. - date_attend), .names="time2_{col}")) %>%
		rename_with(~ gsub("icdDate_", "", .), .cols=starts_with("time2_icdDate_"))
		dat1 <- dat1 %>% select(all_of(grep("^(bb_|time2_)", names(dat1), value=TRUE)))
	cor_matrix <- matrix(NA, nrow=length(Xs), ncol=length(Ys))
	rownames(cor_matrix) <- Xs; colnames(cor_matrix) <- Ys
	for (i in 1:length(Xs)) {
		for (j in 1:length(Ys)) {
			cor_matrix[i,j] <- cor(dat1[[Xs[i]]], dat1[[time2_Ys[j]]], method="spearman", use="pairwise.complete.obs")
		}
	}
col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
Heatmap( cor_matrix, col=col_fun, column_title="diseases onset", row_title="blood biochemistry", heatmap_legend_param=list(title="r2", at=c(-1, 0, 1), labels=c("-1", "0", "1")))

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 桑基图 ⛵
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(networkD3)
	nodes <- data.frame(name=c(
		"FA10", "APOB", "APOA1", "APOA2", "ITB3",
		"Retinoid metabolism", "PPARA activation", "Caspase activation", "Chylomicron assembly", "ROS Detoxification",
		"Lipid metabolism", "Cellular stress response", "Inflammation", "Metabolism of vitamins", "Post-translational modification",
		"Cardiovascular Disease", "Metabolic Syndrome", "Immune Dysfunction", "Neurodegenerative Disease", "General Disease"
	))
	links <- data.frame(
		source=c( 0, 0, 1, 2, 2, 3, 4, 4, 4,5, 5, 6, 7, 7, 8, 9, 9,10, 11, 12, 12, 13, 14, 14, 10, 11),
		target=c(5, 6, 6, 7, 8, 7, 8, 9, 9,10, 11, 11, 12, 12, 13, 13, 14, 15, 16, 16, 17, 18, 18, 19, 15, 19),
		value=c(10, 8, 6, 5, 7, 4, 6, 3, 8, 5, 6, 7, 8, 4, 5, 6, 4, 10, 9, 8, 7, 6, 5, 4, 3, 2) # Flow values: importance of connections
	)
	sankeyNetwork(Links=links, Nodes=nodes, Source="source", Target="target", Value="value", NodeID="name", units="Importance", fontSize=12, nodeWidth=30, sinksRight=TRUE)
library("ggsankeyfier")
library(ggplot2)
data("ecosystem_services")
ggplot(ecosystem_services_pivot1, aes(x=stage, y=RCSES, group=node, connector=connector, edge_id=edge_id, fill=node)) + 
	geom_sankeynode(v_space="auto") + geom_sankeyedge(v_space="auto")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 众多因素 vivid 交互作用分析 💃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	dat1$Y_yes <- as.factor(dat1$Y_yes); dat1$follow_years <- NULL
	fit.rforest <- randomForest(Y_yes ~ ., na.action=na.omit, data=dat1)
	#fit.rforest <- rfsrc(surv.obj ~ ., data=dat1)
	fit.vivi <- vivi(fit=fit.rforest, response="Y_yes", data=dat1); print(fit.vivi, digits=1)
	pdf(paste(Y,X,'heatmap.pdf',sep='.')); print(viviHeatmap(mat=fit.vivi)); dev.off()
	pdf(paste(Y,X,'network.pdf',sep='.')); print(viviNetwork(mat=fit.vivi)); dev.off()
