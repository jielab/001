pacman::p_load(data.table, tidyverse, lubridate, survival, rcssci)

dir0='D:'
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- readRDS(file=paste0(dir0, "/data/ukb/phe/Rdata/all.plus.rds")); dat0$bb_shbg.tp53.rs1042522_C <- NULL
	covs_else <- "age sex bmi smoke_status alcohol_status PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()
	covs_bald <- "age bmi smoke_status alcohol_status bb_TES bb_SHBG bb_VITD bb_IGF1 PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()
	hist(dat0$date_attend, breaks="months", freq=TRUE); hist(dat0$icdDate_sle, breaks="months", freq=TRUE)
	dat0 %>% drop_na(date_attend, icdDate_sle) %>% nrow()
	subset(dat0, prot.yes==1 & icdDate_sle.2 < date_attend) %>% nrow() # 🏮


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table 1: 基本信息
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=="White") %>% mutate(
	bald12=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==2, 1, NA))),
	bald13=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==3, 1, NA))),
	bald14=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==4, 1, NA))),
	bald134=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald %in% 3:4, 1, NA))),
	bald1234=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald %in% 2:4, 1, NA))),
	sexf_cat=factor(ifelse(age_fsex <16, "early", ifelse(age_fsex<22, "average", "late")),levels=c("early", "average", "late")),
	sexp_cat=factor(ifelse(sex_partner <=1, "low", ifelse(sex_partner<=10, "average", "high")),levels=c("low", "average", "high"))
)
table(dat$sex, dat$bald12); table(dat$sexp_cat)

dat1 <- dat %>% filter(sex==1) %>% drop_na(bald, prot_eda2r) 
	dat1 %>% group_by(bald) %>% summarise(count=n(), mean=mean(prot_eda2r, na.rm=TRUE))
	library(ggpubr); ggboxplot(dat1, x="bald", y="prot_eda2r", add="mean_se", color="bald") + # 🏮🥚
		stat_compare_means(method="anova", label="p.format", label.x=1.5, label.y=min(dat1$prot_eda2r) - (max(dat1$prot_eda2r) * 0.1)) + 
		stat_compare_means(aes(label=paste0("p=", round(after_stat(p), 3))), method="t.test", ref.group="1", comparisons=list(c("1", "2"), c("1", "3"), c("1", "4"))) +
		stat_summary(fun="mean", geom="text", aes(label=round(after_stat(y), 2)), vjust=-1.5, color="black") +theme_minimal()

dat1 <- subset(dat, sex==1)
	library(compareGroups); descrTable(formula=as.formula(paste(c("bald ~ ", paste(covs_bald, collapse="+")), collapse="")), data=dat1, digits=1)
	library(gtsummary); dat1 %>% select(all_of(c("bald", covs_bald))) %>% tbl_summary(by=bald, missing="no") %>% add_p()
	Y="bald12"; fit <- glm(formula=as.formula(paste(c(Y, " ~ ", paste(covs_bald, collapse="+")), collapse="")), data=dat1, family="binomial"); summary(fit)
	library(forestmodel); forest_model(fit)
	# fit <- lm(age_fsex ~ factor(bald) +..., data=dat); 
	# fit <- nnet::multinom(bald ~ covs_bald, data=dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 批量运行🏃‍Proteome🥚的表型关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ys <- c("bald12", "bald13", "bald14", "bald134", "bald1234") #  
for (Y in Ys) {
	res <- data.frame(Y=character(), X=character(), BETA=numeric(), SE=numeric(), P=numeric(), stringsAsFactors=FALSE)
	print(Y); if (grepl("bald", Y)) covs=covs_bald else covs=covs_else
	for (X in grep("bb_|bc_|hla_|prot_|met_", names(dat), value=T)) {
	
	print(X)
	if (sum(!is.na(dat[[X]])) < 10000) next
	if (paste0("icdCt_", Y) %in% names(dat)) { # 存在发病日期，cox分析 🔫
		dat1 <- dat %>% mutate(
			Y_date=dat[[paste0('icdDate_',Y)]], Y_yes=ifelse(is.na(Y_date), 0,1),
			X=dat[[X]],
			#X_qt=cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0("q",1:5)),
			#X_qt=ifelse(X_qt=="q1", "low", ifelse(X_qt=="q5", "high", "middle")),
			follow_end_day=fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date("2021-12-31")))),
			follow_years=(as.numeric(follow_end_day) - as.numeric(date_attend)) / 365.25
		) %>% filter(follow_years>0)
		print(sum(dat1$follow_years)) # person-years
		surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
		fit <- coxph(as.formula(paste("surv.obj ~ ", X, "+", paste(covs, collapse="+"))), data=dat1)
		if (!is.na(fit$coef[[X]])) {res <- rbind(res, data.frame(Y=Y, X=X, BETA=rb(summary(fit)$coef[1,1]), SE=rb(summary(fit)$coef[1,3]), P=rp(summary(fit)$coef[1,5])))}	
	} else if (na.omit(all(dat[[Y]]) %in% c(0, 1)))  { # 0/1变量，logistic分析 🔫
		fit <- glm(as.formula(paste(Y, "~", X, "+", paste(covs, collapse="+"))), data=dat, family="binomial")
		if (!is.na(fit$coef[[X]])) {res <- rbind(res, data.frame(Y=Y, X=X, BETA=rb(summary(fit)$coef[2,1]), SE=rb(summary(fit)$coef[2,2]), P=rp(summary(fit)$coef[2,4])))}	
	}
	
	}
	write.table(res, paste0(Y,".assoc.txt"), append=FALSE, quote=FALSE, row.names=FALSE)
	gc()
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Survival🏊‍相关分析结果展示🌳
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
X="bb_VITD"; Y="cad" # 以这个为例，接上面的代码
	rcssci_cox(data=dat1, time="follow_years", y="Y_yes", x="X", covs=c("age","sex"), prob=0.1, filepath= 'D:/tmp') # 🐂
	km.obj <- survfit(surv.obj ~ smoke_status, data=dat1) # KM 是不能带协变量的，M会被作为分层变量
	survdiff(surv.obj ~ smoke_status, data=dat1) # log-rank test
	plot(km.obj, fun=function(x) 1-x)
	survminer::ggsurvplot(km.obj, ylim=c(0,0.2), fun="event", pval=TRUE, risk.table=TRUE, ncensor.plot=TRUE, ggtheme=theme_bw(), palette=c("green","gray","orange"))
	survminer::ggforest(fit, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1)
library(forplo) #🏮
exdf <- data.frame(cbind(OR=c(1.21,0.90,1.02, 1.54,1.32,0.79,1.38,0.85,1.11, 1.58,1.80,2.27),LCI=c(0.82,0.61,0.66,1.08,0.91,0.48,1.15,0.39,0.91,0.99,1.48,0.92),UCI=c(1.79,1.34,1.57, 2.19,1.92,1.32,1.64,1.87,1.34, 2.54,2.19,5.59), groups=c(1,1,1, 2,2,2,2,2,2, 3,3,3)))
	rownames(exdf) <- c('Barry, 2005', 'Frances, 2000', 'Rowley, 1995', 'Biro, 2000', 'Crowe, 2010', 'Harvey, 1996', 'Johns, 2004', 'Parr, 2002', 'Zhang, 2011', 'Flint, 1989', 'Mac Vicar, 1993', 'Turnbull, 1996')
	forplo(exdf[,1:3], groups=exdf$groups, grouplabs=c('Low risk of bias', 'Some concerns', 'High risk of bias'), left.align=TRUE, col=2, ci.edge=FALSE, diamond=c(4,11,15), diamond.col='#b51b35')

## 10年风险计算
	#	XnM=factor(paste(M, X, sep="|"), levels=paste(rep(levels(M),each=3), rep(levels(X),3), sep="|")),
	fit <- coxph(surv.obj ~ X+M, data=dat1) 
	dat1$risk <- nricens::get.risk.coxph(mdl, 10) # 每人的risk 🔔
	dat2 <- expand.grid(X=levels(dat1$X), M=levels(dat1$M))
	surv.pred <- survfit(fit, newdata=dat2)
	surv.10y <- summary(surv.pred, times=10)
	dat2 <- dat2 %>% mutate(risk=1-t(surv.10y$surv), ci_lower=t(1-surv.10y$upper), ci_upper=t(1-surv.10y$lower))
	ggplot(dat2, aes(x=M, y =risk, fill=X)) + geom_bar(stat="identity", position=position_dodge(width=0.7), width=0.7) +
		geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
		geom_text(aes(label=sprintf("%.2f", risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
		labs(x="X labels", y="10-year risk (%)", title="") + scale_fill_manual(values=c("green", "gray", "orange"), name="Legend name") + theme_minimal() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 表型关联和cisMr结果比较
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phes <- c("bald12", "bald13", "bald14", "bald134", "bald1234")
anas <- c("phe", "clp08", "cml08", "g1", "g2")
dat.phe <- read.table(paste0(dir0, '/analysis/assoc.sum/xian/all.assoc.txt'), header=TRUE) %>%
	mutate(X=ifelse(X %like% "prot_", gsub("PROT_", "prot_", toupper(X)), X), Analysis="phe") 
	Y="bald12"; valcano(Y, dat.phe[dat.phe$Y==Y & dat.phe$X %like% "prot_", ], "X", "BETA", "P", 0.05, "", "") # 👀
dat.mr <- read.table(paste0(dir0, '/analysis/assoc.sum/xian/all.mr.log'), header=TRUE)[,c(1:4,6:8)] %>% rename(P=P.ivw) %>% 
	mutate(X_new=ifelse(Analysis %like% "way2", Y, X), Y_new=ifelse(Analysis %like% "way2", X, Y), X=X_new, Y=Y_new) %>% dplyr::select(-X_new, -Y_new)
dat.cis <- read.table(paste0(dir0, '/analysis/assoc.sum/xian/all.cisMr.log'), header=TRUE)[,-5] 
dat <- rbind(dat.mr, dat.cis) %>% mutate(X=paste0("prot_", X), Y=gsub("pancre", "pancre_cancer", Y)) %>%
	mutate(Analysis=paste(Analysis, sprintf("%02d", p_t), sep='.')) %>% dplyr::select(-p_t) %>% 
	mutate(Analysis=recode(Analysis, "gwIV.way1.08"="g1", "gwIV.way2.08"="g2", "cis.clump.06"="clp06", "cis.clump.08"="clp08", "Mr.08"="cml08", "Mr.06"="cml06"))
dat <- dat[, colnames(dat.phe)] # 🏮
dat <- rbind(dat, dat.phe) 
dat <- dat %>% 
	pivot_wider(names_from=c(Y, Analysis), values_from=c(BETA, SE, P), names_glue="{Y}.{.value}.{Analysis}"); names(dat)
	dat <- dat %>% select(-matches("06$|SE")) 
	cols <- apply(expand.grid(c("BETA","P"), anas, phes)[,c(3,1,2)], 1, function(x) paste(x, collapse = "."))
	cols <- c("X", cols); dat <- dat[, cols] # 🏮
dat.coloc <- read.table(paste0(dir0, '/analysis/assoc.sum/xian/all.coloc.log'), header=TRUE)[,c(1,2,8)] %>%
	mutate(X=paste0("prot_",X)) %>% 
	pivot_wider(names_from=Y, values_from=H4, names_glue="{Y}.{.value}"); names(dat.coloc)
dat <- merge(dat, dat.coloc, by="X", all=TRUE) 
	valcano("bald12", dat, "X", "bald12.BETA.g1", "bald12.P.g1", 0.05, "", "") # 👀
	for (i in phes) { for (j in anas) {
		assign(paste0("plot.", i,".",j), valcano(i, dat, "X", paste0(i,".BETA.",j), paste0(i,".P.",j), 0.05, '', '')) # 最后两个空白表示 label_x, label_y
	}}
	plots <- mget( outer(phes, anas, function(i, j) paste0("plot.", i, ".", j)) |> as.vector() )
	gridExtra::grid.arrange(grobs=plots, ncol=length(phes), nrow=length(anas))
	bbplot("bald12 Clump^cML BETA", dat, "X", "bald12.BETA.clp08", "bald12.BETA.cml08", "yes", "bald12.P.clp08", "bald12.P.cml08", "no", 0.05)	
	bbplot("bald12 Clump^cML P", dat, "X", "bald12.BETA.clp08", "bald12.BETA.cml08", "no", "bald12.P.clp08", "bald12.P.cml08", "yes", 0.05)	
	bbplot("bald12 Phe^Mr BETA", dat, "X", "bald12.BETA.phe", "bald12.BETA.cml08", "yes", "bald12.P.phe", "bald12.P.cml08", "no", 0.05)	
	bbplot("bald12 Phe^Mr BETA", dat, "X", "bald12.BETA.phe", "bald12.BETA.cml08", "no", "bald12.P.phe", "bald12.P.cml08", "yes", 0.05)	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prioritize 挑选🏮🏄‍
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppp <- read.table(paste0(dir0, '/data/gwas/ppp/map_3k_v1.tsv'), sep="\t", quote="", header=TRUE, fill=TRUE)[,c(1,3,5)]
drug <- read.table(paste0(dir0, '/files/prot-drug.txt'), sep="\t", header=TRUE, fill=TRUE)
drug <- merge(ppp, drug, by.x="UniProt", by.y="UNIPROT_ACCESSION") %>% 
	mutate(drug=ifelse(APPROVED_DRUG_TARG_CONF !="", "approved", ifelse(CLINICAL_DRUG_TARG_CONF !="", "clinical", ifelse(DRUG_POTENTIAL_TARGET !="", "potential", NA)))) %>%
	select(ID, olink_target_fullname, drug) %>% mutate(ID=paste0("prot_", ID)) %>% arrange("drug")
	names(drug) <- c("X", "olink_target", "drug")
	drug[duplicated(drug$X) | duplicated(drug$X, fromLast=TRUE), ]
	drug <- drug[!duplicated(drug$X), ]
dat <- merge(dat, drug, by="X", all.x=TRUE)
write.table(dat, "analysis.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome与常见病的热力图 🌋
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=="White")
library(GGally)
	ggpairs(tips, mapping=aes(color=sex), columns=c("total_bill", "time", "tip"))
	data(tips, package="reshape"); ggpairs(tips[, c(1, 3, 4, 2)], upper=list(continuous="density", combo="box_no_facet"), lower=list(continuous="points", combo="dot_no_facet"))
library(ComplexHeatmap) # BiocManager::install(version="3.20")
	Ys <- c("asthma", "breast_cancer", "cad", "cancer", "circulatory", "copd", "covid", "lung_cancer", "metabolic", "mi", "pancreatic_cancer", "pe", "prostate_cancer", "ra", "stroke", "t1dm", "t2dm", "varicose", "vte")
		icdDate_Ys <- paste0("icdDate_", Ys)
		time2_Ys <- paste0("time2_", Ys)
	Xs <- grep("bb_", names(dat), value=TRUE) 
		dat1 <- subset(dat, select=Xs[1:10])[1:10000,]; ggpairs(dat1)+theme_bw()
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
	library(circlize); col_fun <- colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))
	Heatmap(cor_matrix, col=col_fun, column_title="diseases onset", row_title="blood biochemistry", heatmap_legend_param=list(title="r2", at=c(-1, 0, 1), labels=c("-1", "0", "1")))

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mediation桑基图 ⛵
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(networkD3, ggsankeyfier, ggplot2)
data("ecosystem_services")
	ggplot(ecosystem_services_pivot1, aes(x=stage, y=RCSES, group=node, connector=connector, edge_id=edge_id, fill=node)) + 
	geom_sankeynode(v_space="auto") + geom_sankeyedge(v_space="auto")
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 众多因素 vivid 交互作用分析 💃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(vivid)
	dat1$Y_yes <- as.factor(dat1$Y_yes); dat1$follow_years <- NULL
	fit.rforest <- randomForest(Y_yes ~ ., na.action=na.omit, data=dat1)
	#fit.rforest <- rfsrc(surv.obj ~ ., data=dat1)
	fit.vivi <- vivi(fit=fit.rforest, response="Y_yes", data=dat1); print(fit.vivi, digits=1)
	pdf(paste(Y,X,'heatmap.pdf',sep='.')); print(viviHeatmap(mat=fit.vivi)); dev.off()
	pdf(paste(Y,X,'network.pdf',sep='.')); print(viviNetwork(mat=fit.vivi)); dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 机器🤖学习📚
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 https://cran.r-project.org/web/packages/csmpv/vignettes/csmpv_vignette.html
# remotes::install_github("ajiangsfu/csmpv",force=TRUE)
library(csmpv)
setwd("D:/analysis/ai"); set.seed(12345) # temp_dir=tempdir(); knitr::opts_knit$set(root.dir=temp_dir)
data("datlist", package="csmpv")
	dat.t <- datlist$training
	dat.v <-datlist$validation
	Xs <- c("B.Symptoms","MYC.IHC","BCL2.IHC", "CD10.IHC","BCL6.IHC", "MUM1.IHC","Male","AgeOver60", "stage3_4","PS1","LDH.Ratio1", "Extranodal1","Bulk10cm","HANS_GCB", "DTI")
	AgeXvars <- setdiff(Xs, "AgeOver60")
	confirmVars(data=dat.t, biomks=Xs, Y="DZsig", outfile="confirmBinary"); bconfirm$fit; bconfirm$allplot[[2]]

modelout <- csmpvModelling(tdat=dat.t, vdat=dat.v, Ybinary="DZsig", varsBinary=Xs, Ycont="Age", varsCont=AgeXvars, time="FFP..Years.", event="Code.FFP", varsSurvival=Xs, outfileName= "All")

# 拟合 ☯
fit1 <- LASSO2(data=dat.t, biomks=Xs, Y="DZsig", outfile="out1"); fit$coef #LASSO2plus(), LASSO_plus() # topN=5
fit2 <- LASSO2_reg(data=dat.t, biomks=Xs, Y="DZsig", outfile="out2")
fit3 <- XGBtraining(data=dat.t, biomks=Xs, Y="DZsig", outfile="out3"); head(fit3$XGBoost_score) # LASSO2_XGBtraining(), LASSO_plus_XGBtraining(), LASSO2plus_XGBtraining()

# ⛅预测	
pred1 <- LASSO2_predict(fit1, newdata=dat.v, outfile="pred1"); head(pred1)
pred2 <- rms_model(fit2$fit, newdata=dat.v, outfile="pred2"); head(pred2)
pred3 <- XGBtraining_predict(fit3, newdata=dat.v, outfile="pred3"); head(pred3)

# 外部验证 ✔
vali1 <- LASSO2_predict(fit1, newdata=dat.v, newY=TRUE, outfile="vali1")
vali2 <- rms_model(fit2$fi, newdata=dat.v, newY=TRUE, outfile="vali2")
vali3 <- XGBtraining_predict(fit3, newdata=dat.v, newY=TRUE, outfile="vali3") 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 机器🤖学习📚
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 https://star-protocols.cell.com/protocols/3440#bib1, 
https://github.com/MobinKhoramjoo/Biomarker-identification-by-multi-omics-analysis
# pacman::p_load(impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, sva, limma, KEGGgraph, siggenes,BiocParallel, MSnbase, multtest, edgeR, fgsea, httr, RBGL, qs) 
# devtools::install_github("xia-lab/MetaboAnalystR", build=TRUE, build_vignettes=TRUE, build_manual=TRUE, force=TRUE)