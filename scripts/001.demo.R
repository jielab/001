pacman::p_load(data.table, tidyverse, lubridate, survival, rcssci, forestmodel)

dir0='D:'
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- readRDS(file=paste0(dir0, '/data/ukb/phe/Rdata/all.plus.rds'))
	covs <- 'age sex bmi smoke_status alcohol_status PC1 PC2 PC3 PC4' %>% strsplit(' ') %>% unlist()
	covs_bald <- "age bmi smoke_status alcohol_status deprivation income sleep_duration sexf_cat sexp_cat bb_TES bb_SHBG bb_VITD bb_IGF1 stiff dash_pts pa_pts smoke_pts sleep_pts bmi_pts nonhdl_pts hba1c_pts bp_pts prot_pts PC1 PC2 PC3 PC4" %>% strsplit(" ") %>% unlist()
	hist(dat0$date_attend, breaks='months', freq=TRUE)
	hist(dat0$icdDate_depress, breaks='months', freq=TRUE)
	subset(dat0, icdDate_depress < date_attend) %>% nrow()

library(naniar); gg_miss_var(subset(dat0, select=grep("sex|bb_", names(dat0), value=TRUE)), facet=sex)
dat <- dat0 %>% select(grep("bb_ALB|bb_APOB|bb_ALP|bb_CYS|bb_HDL|bb_LDL", names(dat0), value=TRUE)) %>% na.omit() %>% dplyr::sample_n(10000)
	car::scatterplotMatrix(dat, spread=FALSE, smoother.args=list(lty=0.1))
	psych::pairs.panels(dat) # 🏮


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 🏊‍ Survival 变量 Y 基本信息
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
date_adj <- function(date.in, date.attend.real, date.attend.fake) {
	return(date.in - (date.attend.real - date.attend.fake))
}
range2 <- function(dat.in) {
	X_mean <- mean(dat.in, na.rm=TRUE)
	X_sd <- sd(dat.in, na.rm=TRUE)
	return(c(X_mean-3*X_sd, X_mean+3*X_sd))
}

dat <- dat0 %>% filter(ethnic_cat=='White', prot.yes==1) 
summary(dat$date_attend); print(date_attend.med <- median(dat$date_attend, na.rm=TRUE)) # 2009-01-09
dat <- dat %>% mutate(
	X = prot_btn3a2, # abo, gast, lrrn1
	X_qt=cut(X, breaks=quantile(X, probs=seq(0,1,0.2), na.rm=T), include.lowest=T, labels=paste0('q',1:5)),
	X_qt=factor(ifelse(X_qt=='q1', 'low', ifelse(X_qt=='q5', 'high', 'middle'))),
	Y_date = icdDate_depress, 
	Y_date2 = date_adj(Y_date, date_attend, date_attend.med), # 🚣
	Y_yes = ifelse(is.na(Y_date), 0,1), 
	follow_date = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date('2021-12-31')))),
	follow_years = (as.numeric(follow_date) - as.numeric(date_attend)) / 365.25,
	before_date = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(birth_date), birth_date, NA)),
	before_years = (as.numeric(date_attend) - as.numeric(before_date)) / 365.25
) 
	dat %>% drop_na(date_attend, Y_date) %>% nrow()
	subset(dat, prot.yes==1 & Y_date <= date_attend) %>% nrow()
	hist(dat$Y_date, breaks="years", freq=TRUE)

par(mar = c(5, 4, 4, 5) + 0.1) # 画图举例🌰
	myhist <- hist(dat$Y_date, breaks="years", freq=TRUE, main="X", ylab="Y (freq)")
	X.avgs <- by(dat$X, cut(dat$Y_date, breaks=as.Date(myhist$breaks)), function(x) mean(x, na.rm = TRUE))
	X.sds <- by(dat$X, cut(dat$Y_date, breaks=as.Date(myhist$breaks)), function(x) sd(x, na.rm=TRUE))
	par(new=T)
	plot(as.Date(myhist$mids), X.avgs, xlim=range(as.Date(myhist$breaks)), ylim=range2(dat$X), pch=16, axes=FALSE, xlab=NA, ylab=NA, cex=1.2, col="red")
	arrows(as.Date(myhist$mids), X.avgs-X.sds, as.Date(myhist$mids), X.avgs+X.sds, angle=90, code=3, length=0.05, col="red")
	axis(side=4); mtext(side=4, line=3, "X", col="red")
	abline(v =date_attend.med, col="blue", lwd=3, lty=2)
	abline(h =mean(dat$X, na.rm=TRUE), col="green", lwd=3, lty=2)

surv.obj <- Surv(time=dat$follow_years, event=dat$Y_yes)
	rcssci_cox(data=dat, time='follow_years', y='Y_yes', x='X', covs=c('age','sex'), prob=0.1, filepath= 'D:/tmp') # 🐂
	km.obj <- survfit(surv.obj ~ smoke_status, data=dat) # KM 是不能带协变量的，M会被作为分层变量
	survdiff(surv.obj ~ smoke_status, data=dat) # log-rank test
	plot(km.obj, fun=function(x) 1-x)
	survminer::ggsurvplot(km.obj, ylim=c(0,0.2), fun='event', pval=TRUE, risk.table=TRUE, ncensor.plot=TRUE, ggtheme=theme_bw(), palette=c('green','gray','orange'))
	fit <- coxph(as.formula(paste('surv.obj ~ X +', paste(covs, collapse='+'))), data=dat)
	survminer::ggforest(fit, main='', cpositions=c(0, 0.1, 0.3), fontsize=1.2, data=dat)

dat1 <- dat %>% mutate(X=X_qt, M=age_cat) %>% drop_na(X, M) # ✝ 10年风险计算
	surv.obj <- Surv(time=dat1$follow_years, event=dat1$Y_yes)
	fit <- coxph(surv.obj ~ X + M, data=dat1) 
	dat2 <- expand.grid(X=levels(dat1$X), M=levels(dat1$M))
	surv.pred <- survfit(fit, newdata=dat2)
	surv.10y <- summary(surv.pred, times=10)
	dat2 <- dat2 %>% mutate(risk=1-t(surv.10y$surv), ci_lower=t(1-surv.10y$upper), ci_upper=t(1-surv.10y$lower))
	# library(nricens); dat1$risk <- nricens::get.risk.coxph(fit, 10) # 可用这一行替代上面4行 🏮
	ggplot(dat2, aes(x=M, y=risk, fill=X)) + geom_bar(stat='identity', position=position_dodge(width=0.7), width=0.7) +
		geom_errorbar(aes(ymin =ci_lower, ymax=ci_upper), width=0.2, position=position_dodge(width=0.7)) +
		geom_text(aes(label=sprintf('%.2f', risk), y=ci_upper), vjust=-1, position=position_dodge(width=0.7), size=2) +
		labs(x='X labels', y='10-year risk (%)', title='') + scale_fill_manual(values=c('green', 'gray', 'orange'), name='Legend name') + theme_minimal() 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 批量运行🏃‍Proteome🥚的表型关联
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Ys <- c('depress', 'bald12', 'bald13', 'bald14', 'bald134', 'bald1234') #  
for (Y in Ys) {
	print(Y)
	res <- data.frame(Analysis=character(), Y=character(), X=character(), BETA=numeric(), SE=numeric(), P=numeric())
	if (grepl('bald', Y)) {covs <- covs_bald}
	if (paste0('icdCt_', Y) %in% names(dat)) { # 存在发病日期，cox分析 🔫
		ana <- 'cox'
		dat1 <- dat %>% mutate(
			Y_date = dat[[paste0('icdDate_',Y)]], 
			Y_yes = ifelse(is.na(Y_date), 0,1),
			follow_date = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(date_lost), date_lost, fifelse(!is.na(date_death), date_death, as.Date('2021-12-31')))),
			follow_years = (as.numeric(follow_date) - as.numeric(date_attend)) / 365.25,
			before_date = fifelse(!is.na(Y_date), Y_date, fifelse(!is.na(birth_date), birth_date, NA)),
			before_years = (as.numeric(date_attend) - as.numeric(before_date)) / 365.25
		)
		dat1.a <- dat1 %>% filter(follow_years>0)
		dat1.b <- dat1 %>% filter(before_years>=0)
		surv.obj.a <- Surv(time=dat1.a$follow_years, event=dat1.a$Y_yes)
		surv.obj.b <- Surv(time=dat1.b$before_years, event=dat1.b$Y_yes)	
	} else if (all(unique(dat[[Y]] %in% c(0, 1, NA)))) {
		dat1 <- dat
		ana <- 'logit' # 0/1变量，logistic分析 🔫
	}

	for (X in grep('bb_|bc_|hla_|prot_|met_', names(dat), value=T)) {
		print(X)
		if (sum(!is.na(dat[[X]])) < 10000) next
		if (ana=='logit') {
			dat1$X = dat1[[X]]
			fit <- glm(as.formula(paste(Y, '~', X, '+', paste(covs, collapse='+'))), data=dat, family='binomial')
			if (!is.na(fit$coef[[X]])) {res <- rbind(res, data.frame(Analysis='logit', Y=Y, X=X, BETA=rb(summary(fit)$coef[2,1]), SE=rb(summary(fit)$coef[2,2]), P=rp(summary(fit)$coef[2,4])))}	
		} else if (ana=='cox') {
			dat1.a$X = dat1.a[[X]]
			dat1.b$X = dat1.b[[X]]
			fit.a <- coxph(as.formula(paste('surv.obj.a ~ ', X, '+', paste(covs, collapse='+'))), data=dat1.a)
			fit.b <- coxph(as.formula(paste('surv.obj.b ~ ', X, '+', paste(covs, collapse='+'))), data=dat1.b)
			if (!is.na(fit.a$coef[[X]])) {res <- rbind(res, data.frame(Analysis='cox.a', Y=Y, X=X, BETA=rb(summary(fit.a)$coef[1,1]), SE=rb(summary(fit.a)$coef[1,3]), P=rp(summary(fit.a)$coef[1,5])))}	
			if (!is.na(fit.b$coef[[X]])) {res <- rbind(res, data.frame(Analysis='cox.b', Y=Y, X=X, BETA=rb(summary(fit.b)$coef[1,1]), SE=rb(summary(fit.b)$coef[1,3]), P=rp(summary(fit.b)$coef[1,5])))}	
		}
	}
	write.table(res, paste0(Y,'.assoc.txt'), append=FALSE, quote=FALSE, row.names=FALSE)
	res1 <- read.table(paste0(Y,'.assoc.txt'), header=TRUE) %>% filter(X %like% 'prot_') %>%
		pivot_wider(names_from=Analysis, values_from=c(BETA, SE, P), names_glue='{.value}.{Analysis}'); names(res1)
		pdf(paste0(Y,'.pdf'))
		valcano(Y, res1[res1$X %like% 'prot_', ], 'X', 'BETA.cox.a', 'P.cox.a', 0.05, '', '') # 👀
		valcano(Y, res1[res1$X %like% 'prot_', ], 'X', 'BETA.cox.b', 'P.cox.b', 0.05, '', '') # 👀
		dev.off()
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 表型关联和cisMr结果比较
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
phes <- c('bald12', 'bald13', 'bald14', 'bald134', 'bald1234')
anas <- c('phe', 'clp08', 'cml08', 'g1', 'g2')
dat.phe <- read.table(paste0(dir0, '/analysis/assoc.sum/all.assoc.txt'), header=TRUE) %>%
	mutate(X=ifelse(X %like% 'prot_', gsub('PROT_', 'prot_', toupper(X)), X), Analysis='phe') 
	Y='bald12'; valcano(Y, dat.phe[dat.phe$Y==Y & dat.phe$X %like% 'prot_', ], 'X', 'BETA', 'P', 0.05, '', '') # 👀
dat.mr <- read.table(paste0(dir0, '/analysis/assoc.sum/all.mr.log'), header=TRUE)[,c(1:4,6:8)] %>% rename(P=P.ivw) %>% 
	mutate(X_new=ifelse(Analysis %like% 'way2', Y, X), Y_new=ifelse(Analysis %like% 'way2', X, Y), X=X_new, Y=Y_new) %>% dplyr::select(-X_new, -Y_new)
dat.cis <- read.table(paste0(dir0, '/analysis/assoc.sum/all.cisMr.log'), header=TRUE)[,-5] 
dat <- rbind(dat.mr, dat.cis) %>% mutate(X=paste0('prot_', X), Y=gsub('pancre', 'pancre_cancer', Y)) %>%
	mutate(Analysis=paste(Analysis, sprintf('%02d', p_t), sep='.')) %>% dplyr::select(-p_t) %>% 
	mutate(Analysis=recode(Analysis, 'gwIV.way1.08'='g1', 'gwIV.way2.08'='g2', 'cis.clump.06'='clp06', 'cis.clump.08'='clp08', 'Mr.08'='cml08', 'Mr.06'='cml06'))
dat <- dat[, colnames(dat.phe)] # 🏮
dat <- rbind(dat, dat.phe) 
dat <- dat %>% 
	pivot_wider(names_from=c(Y, Analysis), values_from=c(BETA, SE, P), names_glue='{Y}.{.value}.{Analysis}'); names(dat)
	dat <- dat %>% select(-matches('06$|SE')) 
	cols <- apply(expand.grid(c('BETA','P'), anas, phes)[,c(3,1,2)], 1, function(x) paste(x, collapse='.'))
	cols <- c('X', cols); dat <- dat[, cols] # 🏮
dat.coloc <- read.table(paste0(dir0, '/analysis/assoc.sum/all.coloc.log'), header=TRUE)[,c(1,2,8)] %>%
	mutate(X=paste0('prot_',X)) %>% 
	pivot_wider(names_from=Y, values_from=H4, names_glue='{Y}.{.value}'); names(dat.coloc)
dat <- merge(dat, dat.coloc, by='X', all=TRUE) 
	valcano('bald12', dat, 'X', 'bald12.BETA.g1', 'bald12.P.g1', 0.05, '', '') # 👀
	for (i in phes) { for (j in anas) {
		assign(paste0('plot.', i,'.',j), valcano(i, dat, 'X', paste0(i,'.BETA.',j), paste0(i,'.P.',j), 0.05, '', '')) # 最后两个空白表示 label_x, label_y
	}}
	plots <- mget( outer(phes, anas, function(i, j) paste0('plot.', i, '.', j)) |> as.vector() )
	gridExtra::grid.arrange(grobs=plots, ncol=length(phes), nrow=length(anas))
	bbplot('bald12 Clump^cML BETA', dat, 'X', 'bald12.BETA.clp08', 'bald12.BETA.cml08', 'yes', 'bald12.P.clp08', 'bald12.P.cml08', 'no', 0.05)	
	bbplot('bald12 Clump^cML P', dat, 'X', 'bald12.BETA.clp08', 'bald12.BETA.cml08', 'no', 'bald12.P.clp08', 'bald12.P.cml08', 'yes', 0.05)	


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prioritize 挑选🏮🏄‍
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ppp <- read.table(paste0(dir0, '/data/gwas/ppp/map_3k_v1.tsv'), sep='\t', quote='', header=TRUE, fill=TRUE)[,c(1,3,5)]
drug <- read.table(paste0(dir0, '/files/prot-drug.txt'), sep='\t', header=TRUE, fill=TRUE)
drug <- merge(ppp, drug, by.x='UniProt', by.y='UNIPROT_ACCESSION') %>% 
	mutate(drug=ifelse(APPROVED_DRUG_TARG_CONF !='', 'approved', ifelse(CLINICAL_DRUG_TARG_CONF !='', 'clinical', ifelse(DRUG_POTENTIAL_TARGET !='', 'potential', NA)))) %>%
	select(ID, olink_target_fullname, drug) %>% mutate(ID=paste0('prot_', ID)) %>% arrange('drug')
	names(drug) <- c('X', 'olink_target', 'drug')
	drug[duplicated(drug$X) | duplicated(drug$X, fromLast=TRUE), ]
	drug <- drug[!duplicated(drug$X), ]
dat <- merge(dat, drug, by='X', all.x=TRUE)
saveRDS(dat, "bald.rds")
write.table(dat, 'bald.tsv', sep='\t', append=FALSE, quote=FALSE, row.names=FALSE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# XY之间的热力图🌅
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat0 %>% filter(ethnic_cat=='White')
library(GGally)
	ggpairs(tips, mapping=aes(color=sex), columns=c('total_bill', 'time', 'tip'))
	data(tips, package='reshape'); ggpairs(tips[, c(1, 3, 4, 2)], upper=list(continuous='density', combo='box_no_facet'), lower=list(continuous='points', combo='dot_no_facet'))
library(ComplexHeatmap) # BiocManager::install(version='3.20')
	Ys <- c('asthma', 'breast_cancer', 'cad', 'cancer', 'circulatory', 'copd', 'covid', 'lung_cancer', 'metabolic', 'mi', 'pancreatic_cancer', 'pe', 'prostate_cancer', 'ra', 'stroke', 't1dm', 't2dm', 'varicose', 'vte')
		icdDate_Ys <- paste0('icdDate_', Ys)
		time2_Ys <- paste0('time2_', Ys)
	Xs <- grep('bb_', names(dat), value=TRUE) 
		dat1 <- subset(dat, select=Xs[1:10])[1:10000,]; ggpairs(dat1)+theme_bw()
	dat1 <- dat %>% mutate(across(.cols=intersect(names(dat)[matches('icdDate_', ignore.case=TRUE)], icdDate_Ys), .fns=~ as.numeric(. - date_attend), .names='time2_{col}')) %>%
		rename_with(~ gsub('icdDate_', '', .), .cols=starts_with('time2_icdDate_'))
		dat1 <- dat1 %>% select(all_of(grep('^(bb_|time2_)', names(dat1), value=TRUE)))
	cor_matrix <- matrix(NA, nrow=length(Xs), ncol=length(Ys))
		rownames(cor_matrix) <- Xs; colnames(cor_matrix) <- Ys
		for (i in 1:length(Xs)) {
		for (j in 1:length(Ys)) {
			cor_matrix[i,j] <- cor(dat1[[Xs[i]]], dat1[[time2_Ys[j]]], method='spearman', use='pairwise.complete.obs')
		}
		}
	library(circlize); col_fun <- colorRamp2(c(-0.5, 0, 0.5), c('blue', 'white', 'red'))
	Heatmap(cor_matrix, col=col_fun, column_title='diseases onset', row_title='blood biochemistry', heatmap_legend_param=list(title='r2', at=c(-1, 0, 1), labels=c('-1', '0', '1')))

  
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mediation桑基图 ⛵
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pacman::p_load(networkD3, ggsankeyfier, ggplot2)
data('ecosystem_services')
	ggplot(ecosystem_services_pivot1, aes(x=stage, y=RCSES, group=node, connector=connector, edge_id=edge_id, fill=node)) + 
	geom_sankeynode(v_space='auto') + geom_sankeyedge(v_space='auto')
nodes <- data.frame(name=c(
	'FA10', 'APOB', 'APOA1', 'APOA2', 'ITB3',
	'Retinoid metabolism', 'PPARA activation', 'Caspase activation', 'Chylomicron assembly', 'ROS Detoxification',
	'Lipid metabolism', 'Cellular stress response', 'Inflammation', 'Metabolism of vitamins', 'Post-translational modification',
	'Cardiovascular Disease', 'Metabolic Syndrome', 'Immune Dysfunction', 'Neurodegenerative Disease', 'General Disease'
))
links <- data.frame(
	source=c( 0, 0, 1, 2, 2, 3, 4, 4, 4,5, 5, 6, 7, 7, 8, 9, 9,10, 11, 12, 12, 13, 14, 14, 10, 11),
	target=c(5, 6, 6, 7, 8, 7, 8, 9, 9,10, 11, 11, 12, 12, 13, 13, 14, 15, 16, 16, 17, 18, 18, 19, 15, 19),
	value=c(10, 8, 6, 5, 7, 4, 6, 3, 8, 5, 6, 7, 8, 4, 5, 6, 4, 10, 9, 8, 7, 6, 5, 4, 3, 2) # Flow values: importance of connections
)
sankeyNetwork(Links=links, Nodes=nodes, Source='source', Target='target', Value='value', NodeID='name', units='Importance', fontSize=12, nodeWidth=30, sinksRight=TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 众多因素 vivid 交互作用分析 💃
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(vivid)
	dat1$Y_yes <- as.factor(dat1$Y_yes); dat1$follow_years <- NULL
	fit.rforest <- randomForest(Y_yes ~ ., na.action=na.omit, data=dat1)
	#fit.rforest <- rfsrc(surv.obj ~ ., data=dat1)
	fit.vivi <- vivi(fit=fit.rforest, response='Y_yes', data=dat1); print(fit.vivi, digits=1)
	pdf(paste(Y,X,'heatmap.pdf',sep='.')); print(viviHeatmap(mat=fit.vivi)); dev.off()
	pdf(paste(Y,X,'network.pdf',sep='.')); print(viviNetwork(mat=fit.vivi)); dev.off()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 机器🤖学习📚
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 https://cran.r-project.org/web/packages/csmpv/vignettes/csmpv_vignette.html
# remotes::install_github('ajiangsfu/csmpv',force=TRUE)
library(csmpv)
setwd('D:/analysis/ai'); set.seed(12345) # temp_dir=tempdir(); knitr::opts_knit$set(root.dir=temp_dir)
data('datlist', package='csmpv')
	dat.t <- datlist$training
	dat.v <-datlist$validation
	Xs <- c('B.Symptoms','MYC.IHC','BCL2.IHC', 'CD10.IHC','BCL6.IHC', 'MUM1.IHC','Male','AgeOver60', 'stage3_4','PS1','LDH.Ratio1', 'Extranodal1','Bulk10cm','HANS_GCB', 'DTI')
	AgeXvars <- setdiff(Xs, 'AgeOver60')
	confirmVars(data=dat.t, biomks=Xs, Y='DZsig', outfile='confirmBinary'); bconfirm$fit; bconfirm$allplot[[2]]

modelout <- csmpvModelling(tdat=dat.t, vdat=dat.v, Ybinary='DZsig', varsBinary=Xs, Ycont='Age', varsCont=AgeXvars, time='FFP..Years.', event='Code.FFP', varsSurvival=Xs, outfileName= 'All')

# 拟合 ☯
fit1 <- LASSO2(data=dat.t, biomks=Xs, Y='DZsig', outfile='out1'); fit$coef #LASSO2plus(), LASSO_plus() # topN=5
fit2 <- LASSO2_reg(data=dat.t, biomks=Xs, Y='DZsig', outfile='out2')
fit3 <- XGBtraining(data=dat.t, biomks=Xs, Y='DZsig', outfile='out3'); head(fit3$XGBoost_score) # LASSO2_XGBtraining(), LASSO_plus_XGBtraining(), LASSO2plus_XGBtraining()

# ⛅预测	
pred1 <- LASSO2_predict(fit1, newdata=dat.v, outfile='pred1'); head(pred1)
pred2 <- rms_model(fit2$fit, newdata=dat.v, outfile='pred2'); head(pred2)
pred3 <- XGBtraining_predict(fit3, newdata=dat.v, outfile='pred3'); head(pred3)

# 外部验证 ✔
vali1 <- LASSO2_predict(fit1, newdata=dat.v, newY=TRUE, outfile='vali1')
vali2 <- rms_model(fit2$fi, newdata=dat.v, newY=TRUE, outfile='vali2')
vali3 <- XGBtraining_predict(fit3, newdata=dat.v, newY=TRUE, outfile='vali3') 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 机器🤖学习📚
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 参考 https://star-protocols.cell.com/protocols/3440#bib1, 
https://github.com/MobinKhoramjoo/Biomarker-identification-by-multi-omics-analysis
# pacman::p_load(impute, pcaMethods, globaltest, GlobalAncova, Rgraphviz, preprocessCore, genefilter, sva, limma, KEGGgraph, siggenes,BiocParallel, MSnbase, multtest, edgeR, fgsea, httr, RBGL, qs) 
# devtools::install_github('xia-lab/MetaboAnalystR', build=TRUE, build_vignettes=TRUE, build_manual=TRUE, force=TRUE)