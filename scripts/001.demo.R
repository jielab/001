pacman::p_load(data.table, tidyverse, lubridate, survival, rcssci)

dir0='D:'
source(paste0(dir0, '/scripts/f/main.f.R'))

dat0 <- readRDS(file=paste0(dir0, "/data/ukb/phe/Rdata/all.plus.rds")); dat0$bb_shbg.tp53.rs1042522_C <- NULL
covs_else <- "age sex bmi smoke_status alcohol_status PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()
covs_bald <- "age bmi smoke_status alcohol_status bb_TES bb_SHBG bb_VITD bb_IGF1 PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Table 1: 基本信息及Pheno关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=="White") %>% mutate(
	bald12=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==2,1,NA))),
	bald13=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==3,1,NA))),
	bald14=ifelse(sex==0, NA, ifelse(bald==1,0, ifelse(bald==4,1,NA))),
	sexf_cat=factor(ifelse(age_fsex <16, "early", ifelse(age_fsex<22, "average", "late")),levels=c("early", "average", "late")),
	sexp_cat=factor(ifelse(sex_partner <=1, "low", ifelse(sex_partner<=10, "average", "high")),levels=c("low", "average", "high"))
)
table(dat$sex, dat$bald12); table(dat$sexp_cat)
dat.male <- subset(dat, sex==1)
	library(compareGroups); descrTable(formula=as.formula(paste(c("bald ~ ", paste(covs_bald, collapse="+")), collapse="")), data=dat, digits=1)
	library(gtsummary); dat %>% filter(sex==1) %>% select(all_of(c("bald", covs_bald))) %>% tbl_summary(by=bald, missing="no") %>% add_p()
	Y="bald12"; fit <- glm(formula=as.formula(paste(c(Y, " ~ ", paste(covs_bald, collapse="+")), collapse="")), data=dat, family="binomial"); summary(fit)
	library(forestmodel); forest_model(fit)
	# fit <- lm(age_fsex ~ factor(bald) +..., data=dat); 
	# fit <- nnet::multinom(bald ~ covs_bald, data=dat)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteome🥚的表型关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ys <- c("sle", "pancre_cancer", "bald12", "bald13", "bald14") #  
res <- data.frame(Y=character(), X=character(), BETA=numeric(), SE=numeric(), P=numeric(), stringsAsFactors=FALSE)
for (Y in Ys) {
	print(Y); if (grepl("bald", Y)) covs=covs_bald else covs=covs_else
	for (X in grep("bb_|bc_|hla_|prot_|met_", names(dat), value=T)) {
	
	print(X)
	if (sum(!is.na(dat[[X]])) < 10000) next
	if (paste0("icdCt_", Y) %in% names(dat)) { # 存在发病日期，cox分析 🔫
		dat1 <- dat %>% mutate(
			Y_date=dat[[paste0('icdDate_',Y)]], Y_yes=ifelse(is.na(Y_date), 0,1),
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
# 核心结果展示🌳
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rcssci_cox(data=dat1, time="follow_years", y="Y_yes", x="X", covs=c("age","sex"), prob=0.1, filepath= 'D:/tmp') # 🐂

fit.cox %>% gtsummary::tbl_regression(exponentiate=TRUE) %>% plot()
	km.obj <- survfit(surv.obj ~ X_qt, data=dat1) # KM 是不能带协变量的，M会被作为分层变量
	# survdiff(surv.obj ~ X_qt, data=dat1) # log-rank test
	# plot(km.obj, fun=function(x) 1-x)
	survminer::ggsurvplot(km.obj, ylim=c(0,0.2), fun="event", pval=TRUE, risk.table=TRUE, ncensor.plot=TRUE, ggtheme=theme_bw(), palette=c("green","gray","orange"))
	survminer::ggforest(fit.cox, main="", cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat1)

library(forplo) #🏮
exdf <- data.frame(cbind(OR=c(1.21,0.90,1.02, 1.54,1.32,0.79,1.38,0.85,1.11, 1.58,1.80,2.27),LCI=c(0.82,0.61,0.66,1.08,0.91,0.48,1.15,0.39,0.91,0.99,1.48,0.92),UCI=c(1.79,1.34,1.57, 2.19,1.92,1.32,1.64,1.87,1.34, 2.54,2.19,5.59), groups=c(1,1,1, 2,2,2,2,2,2, 3,3,3)))
	rownames(exdf) <- c('Barry, 2005', 'Frances, 2000', 'Rowley, 1995', 'Biro, 2000', 'Crowe, 2010', 'Harvey, 1996', 'Johns, 2004', 'Parr, 2002', 'Zhang, 2011', 'Flint, 1989', 'Mac Vicar, 1993', 'Turnbull, 1996')
	forplo(exdf[,1:3], groups=exdf$groups, grouplabs=c('Low risk of bias', 'Some concerns', 'High risk of bias'), left.align=TRUE, col=2, ci.edge=FALSE, diamond=c(4,11,15), diamond.col='#b51b35')

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
# 表型关联和cisMr结果比较
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat.phe <- read.table("all.assoc.txt", header=TRUE) %>% 
	mutate(X = paste0("prot_", toupper(sub("prot_", "", X))), Analysis="phe") 
	Y="bald14"; valcano(Y, dat.phe[dat.phe$Y==Y,], "X", "BETA", "P", 0.05)
dat.mr <- read.table(paste0(dir0, '/analysis/assoc.sum/all.mr.log'), header=TRUE)[,c(1:4,6:8)] %>% rename(P=P.ivw)
dat.cis <- read.table(paste0(dir0, '/analysis/assoc.sum/all.cisMr.log'), header=TRUE)[,-5] 
dat.mr <- rbind(dat.mr, dat.cis) %>% mutate(X=paste0("prot_", X)) %>%
	mutate(Analysis = paste(Analysis, sprintf("%02d", p_t), sep='.')) %>% select(-p_t) %>% 
	mutate(Analysis = recode(Analysis, "gwIV.way1.08"="g1", "gwIV.way2.08"="g2", "cis.clump.06"="clp06", "cis.clump.08"="clp08", "Mr.08"="cml08", "Mr.06"="cml06"))
dat.mr <- dat.mr[, colnames(dat.phe)] # 🏮
dat <- rbind(dat.mr, dat.phe) %>% filter(Y %like% "bald") %>% 
	pivot_wider(names_from=c(Y, Analysis), values_from=c(BETA, SE, P), names_glue="{Y}.{.value}.{Analysis}") 
	names(dat)
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
	drug[duplicated(drug$X) | duplicated(drug$X, fromLast = TRUE), ]
	drug <- drug[!duplicated(drug$X), ]
dat <- merge(dat, drug, by="X", all.x=TRUE)
write.table(dat, "bald.analysis.txt", sep="\t", append=FALSE, quote=FALSE, row.names=FALSE)


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
# 参考 https://star-protocols.cell.com/protocols/3440#bib1, https://github.com/MobinKhoramjoo/Biomarker-identification-by-multi-omics-analysis
# pacman::p_load(impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, sva, limma, KEGGgraph, siggenes,BiocParallel, MSnbase, multtest, edgeR, fgsea, httr, RBGL, qs) 
# devtools::install_github("xia-lab/MetaboAnalystR", build=TRUE, build_vignettes=TRUE, build_manual=TRUE, force=TRUE)
pacman::p_load(readxl, MetaboAnalystR, ComplexHeatmap, circlize)
cytokines <- read_excel("Data.xlsx", sheet = 'Cytokines')
Proteins <- read_excel('Data.xlsx', sheet = 'Proteins')[,-(1:5)]
Metabolites <- read_excel('Data.xlsx', sheet = 'Metabolites')[,-(1:5)]
Data <- cbind(cytokines, Proteins, Metabolites)
	desc_cols =1:4; n_first_mol=5; Data <- cbind(Data[,desc_cols], as.data.frame(sapply(Data[,n_first_mol:length(Data)], function (x) as.numeric(as.character(x)))))
	Data <- replace(Data, Data == 0, NA)
	Transposed_Data <- as.data.frame(t(Data))
	p <-c(); for (i in 1:nrow(Transposed_Data)) { p[i] <- sum(is.na(Transposed_Data[i,]))/ncol(Transposed_Data)}
	Transposed_Data <- Transposed_Data %>% mutate(percent_of_missing_Values= p) %>% select(percent_of_missing_Values, everything())
	missing_value_cut_off = 0.5
	filtered_Transposed_Data <- Transposed_Data %>% filter(percent_of_missing_Values < missing_value_cut_off)
	cleaned_data <- as.data.frame(t(filtered_Transposed_Data))
	cleaned_data <- cleaned_data[-1,]
	row.names(cleaned_data) <- row.names(data)
	eliminated_data <- Transposed_Data %>% filter(percent_of_missing_Values > missing_value_cut_off)
	for (i in n_first_mol:ncol(cleaned_data)){
		for(j in 1:nrow(cleaned_data)){if (is.na(cleaned_data[j,i])=="TRUE"){cleaned_data[j,i] = min(cleaned_data[,i], na.rm = TRUE)}}
	}
	write.table(cleaned_data, "cleaned_data.txt", sep="\t", quote=FALSE, append=FALSE, row.names=FALSE, col.names=TRUE)
mSet <- InitDataObjects(data.type="conc", anal.type="stat", paired=FALSE)
	mSet <- Read.TextData(mSet, filePath="cleaned_data.txt", format="rowu", lbl.type="disc") # 🏮
	mSet <- SanityCheckData(mSet)
	mSet <- PreparePrenormData(mSet)
	mSet <- Normalization(mSet, "NULL", "LogNorm", "MeanCenter", "NULL", ratio=FALSE, ratioNum=20)
mSet <- PCA.Anal(mSet) #Perform PCA
	mSet <- PlotPCAPairSummary(mSet, "pca_pair_0_", format = "png", dpi = 72, width=NA, 5) # Create PCA overview
	mSet <- PlotPCAScree(mSet, "pca_scree_0_", "png", dpi = 72, width=NA, 5) # Create PCA scree plot
	mSet <- PlotPCA2DScore(mSet, "pca_score2d_0_", format = "png", dpi=300, width=NA, 1, 2, 0.95, 0, 0) # Create a 2D PCA score plot
	Fold_change_cut_off= 1.5; P_Value_cut_off = 0.05 
	mSet <- Volcano.Anal(mSet, FALSE, Fold_change_cut_off, 0, F, P_Value_cut_off, TRUE, "fdr")
	mSet <- PlotVolcano(mSet, "volcano_0_", 1, 0, format ="png", dpi=300, width=NA)
	top_mol = 100; mSet <- PlotSubHeatMap(mSet, "heatmap_1_", "png", 300, width=NA,"norm", "row", "euclidean", "ward.D","bwm", 8, "tanova", top_mol, "overview", F, T, T, F, T, T, T)

clinical_variables <- read.csv("Data.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 机器🤖学习📚
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 https://cran.r-project.org/web/packages/csmpv/vignettes/csmpv_vignette.html
# remotes::install_github("ajiangsfu/csmpv",force = TRUE)
library(csmpv)
set.seed(12345); temp_dir = tempdir(); knitr::opts_knit$set(root.dir = temp_dir)
data("datlist", package = "csmpv")
	tdat = datlist$training
	vdat = datlist$validation
	Xvars = c("B.Symptoms","MYC.IHC","BCL2.IHC", "CD10.IHC","BCL6.IHC", "MUM1.IHC","Male","AgeOver60", "stage3_4","PS1","LDH.Ratio1", "Extranodal1","Bulk10cm","HANS_GCB", "DTI")
bconfirm = confirmVars(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "confirmBinary")
	bconfirm$fit; bconfirm$allplot[[2]]
bl = LASSO2(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binaryLASSO2"); bl$coefs
	b2fit = LASSO2plus(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binaryLASSO2plus"); b2fit$fit$coefficients
	bfit = LASSO_plus(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binaryLASSO_plus", topN = 5); bfit$fit$coefficients
	blr = LASSO2_reg(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binaryLASSO2_reg"); blr$fit$coefficients
bxfit = XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binary_XGBoost"); head(bxfit$XGBoost_score)
	blxfit = LASSO2_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binary_LASSO2_XGBoost"); head(blxfit$XGBoost_score)
	blpxfit = LASSO_plus_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig", topN = 5,outfile = "binary_LASSO_plus_XGBoost"); head(blpxfit$XGBoost_score)
	bl2xfit = LASSO2plus_XGBtraining(data = tdat, biomks = Xvars, Y = "DZsig", outfile = "binary_LASSO2plus_XGBoost"); head(bl2xfit$XGBoost_score)

# ⛅预测	
pbl = LASSO2_predict(bl, newdata = vdat, outfile = "pred_LASSO2_binary"); head(pbl)
pblr = rms_model(blr$fit, newdata = vdat, outfile = "pred_LASSO2reg_binary"); head(pblr) # fit$fit b2fit$fit
pbxfit = XGBtraining_predict(bxfit, newdata = vdat, outfile = "pred_XGBoost_binary") blxfit blpxfit bl2xfit

# ✔ (External) Model Validation
vbl = LASSO2_predict(bl, newdata = vdat, newY = TRUE, outfile = "valid_LASSO2_binary")
	vblr = rms_model(blr$fit, newdata = vdat, newY = TRUE, outfile = "valid_LASSO2reg_binary") # bfit$fit 或 b2fit$fit
	vbxfit = XGBtraining_predict(bxfit, newdata = vdat, newY = TRUE, outfile = "valid_XGBoost_binary") # blxfit 或 blpxfit 或 bl2xfit

# 👨‍👩‍👧‍👦📚
modelout = csmpvModelling(tdat = tdat, vdat = vdat, Ybinary = "DZsig", varsBinary = Xvars, Ycont = "Age", varsCont = AgeXvars, time = "FFP..Years.", event = "Code.FFP", varsSurvival = Xvars, outfileName= "all_in_one")
xgobj = XGpred(data = tdat, varsIn = Xvars, nclass=3, vsMethod = "LASSO_plus", topN = 5, time = "FFP..Years.", event = "Code.FFP", outfile = "XGpred", selection = TRUE)
	tdat$XGpred_class2 = xgobj$XGpred_prob_class
	training_risk_confirm2 = confirmVars(data = tdat, biomks = "XGpred_class2", time = "FFP..Years.", event = "Code.FFP", outfile = "training2grps_riskSurvival", outcomeType = "time-to-event")
	training_risk_confirm2[[3]]
	xgNew = XGpred_predict(newdat = vdat, XGpredObj = xgobj)
	vdat$XGpred_class2 = xgNew$XGpred_prob_class
	risk_confirm2 = confirmVars(data = vdat, biomks = "XGpred_class2", time = "FFP..Years.", event = "Code.FFP", outfile = "riskSurvival2grps", outcomeType = "time-to-event")
	risk_confirm2[[3]]