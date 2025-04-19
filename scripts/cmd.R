pacman::p_load(tidyverse, MASS, rrcov, rlang, openxlsx, ggpubr, reshape2, survey, grid, egg, patchwork)

dir0='D:'
source(paste0(dir0, '/scripts/f/cmd.cluster.f.R'))
set.seed(1234)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. UKB cluster 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
var_plot <- c("sbp","dbp","bb_HDL", "bb_LDL", "bb_TC", "bb_TG", "bmi", "whr", "bb_GLU","bb_HBA1C")
	var_level <-  c("Low risk","Low HDL, Low BMI","Obesity, Low cholesterol","Severe obesity","High TG","High cholesterol", "High blood pressure","High blood pressure, High cholesterol","High hyperglycemia","High hyperglycemia, Obesity")
	var_reverse <- c('height', 'bb_HDL', 'bb_eGFR')

arr_labels <- list(
  c("Mid-risk", "Low risk", "Severe hyperglycemia", "High blood pressure", "Severe obesity"),
  c("Low risk", "High blood pressure", "Low DBP, low eGFR", "Mid-risk", "Severe obesity", "Severe hyperglycemia"),
  c("Low DBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk", "Severe hyperglycemia", "Low risk", "Low BMI, high HDL"),
  c("Mid-risk short", "Mid-risk tall", "High blood pressure", "Severe obesity", "Severe hyperglycemia", "Low DBP, low eGFR", "Low risk", "Low BMI, high HDL"),
  c("Mid-risk tall", "Low risk", "High blood pressure", "Mid-risk short", "High cholesterol", "Severe hyperglycemia", "Low BMI, high HDL", "Low DBP, low eGFR", "Severe obesity"),
  c("Low risk", "High blood pressure", "Low DBP, low eGFR", "High heart rate", "Severe hyperglycemia", "Low BMI, high HDL", "Mid-risk tall", "Severe obesity", "High cholesterol", "Mid-risk short"),
  c("Low risk", "Low DBP, low eGFR", "High heart rate", "Mid-risk tall", "Low BMI, high HDL", "Severe hyperglycemia", "High blood pressure", "Mid-risk short", "Severe obesity", "Overweight", "High cholesterol"),
  c("Mid-risk tall", "Overweight", "High cholesterol", "Low risk", "High SBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk short", "Low DBP, low eGFR", "High heart rate", "Low BMI, high HDL", "Severe hyperglycemia")
)

dat0 <- readRDS(file="../data/all.rds")
dat <- dat0 %>% filter(ethnic_cat=="White", !is.na(icdDate_is.2), complete.cases(across(all_of(var_plot)))) 
dat.s <- scale(dat[, var_plot])

for (k in 5:12) {
  km <- kmeans(dat.s, centers = k, iter.max = 10000, nstart = 50, algorithm = "Lloyd")
  dat[[paste0("cluster_", k)]] <- km$cluster
}

plot_clusters(dat, clusters = "cluster_10")
dat$cluster_10 <- factor(dat$cluster_10, levels = var_level, labels =var_label)
dat <- dat[, c("eid", paste0("cluster_", 5:10))]
saveRDS(sankey, file = "cmd.ukb.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. NHANES cluster 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
var_plot <- c("height","bmi","whtr","hba1c","hdl","non_hdl","sbp","dbp","eGFR","pulse")
	var_label <- c("(low) \n Height","BMI","WHtR","HbA1c","(low) \n HDL","Non-HDL","SBP","DBP","Pulse \n rate","(low) \n eGFR")
	var_reverse <- c('height', 'hdl', 'eGFR')

arr_labels <- list(
  c("Mid-risk", "Low risk", "Severe hyperglycemia", "High blood pressure", "Severe obesity"),
  c("Low risk", "High blood pressure", "Low DBP, low eGFR", "Mid-risk", "Severe obesity", "Severe hyperglycemia"),
  c("Low DBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk", "Severe hyperglycemia", "Low risk", "Low BMI, high HDL"),
  c("Mid-risk short", "Mid-risk tall", "High blood pressure", "Severe obesity", "Severe hyperglycemia", "Low DBP, low eGFR", "Low risk", "Low BMI, high HDL"),
  c("Mid-risk tall", "Low risk", "High blood pressure", "Mid-risk short", "High cholesterol", "Severe hyperglycemia", "Low BMI, high HDL", "Low DBP, low eGFR", "Severe obesity"),
  c("Low risk", "High blood pressure", "Low DBP, low eGFR", "High heart rate", "Severe hyperglycemia", "Low BMI, high HDL", "Mid-risk tall", "Severe obesity", "High cholesterol", "Mid-risk short"),
  c("Low risk", "Low DBP, low eGFR", "High heart rate", "Mid-risk tall", "Low BMI, high HDL", "Severe hyperglycemia", "High blood pressure", "Mid-risk short", "Severe obesity", "Overweight", "High cholesterol"),
  c("Mid-risk tall", "Overweight", "High cholesterol", "Low risk", "High SBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk short", "Low DBP, low eGFR", "High heart rate", "Low BMI, high HDL", "Severe hyperglycemia")
)
	
dat0 <- readRDS(file=paste0(dir0, '/data/ukb/phe/Rdata/NHANES_Cleaned_single.RDS')) %>% 
	rename(sbp=sbp_final, dbp=dbp_final) %>% filter(age>=20)
	vars <- grep("clean", names(dat0), value=T); vars1 <- gsub("_clean$", "", vars)
	dat0[,vars1] <- dat0[,vars]; dat0[,vars] <- NULL
dat <- dat0 %>% filter(sex==1) %>% mutate(non_hdl = tc - hdl, whtr = waist/height) %>% drop_na(all_of(vars)) 
	dat.s <- scale(dat[,vars])

km.10 <- kmeans(dat.s, 10, iter.max=10000, nstart=50, algorithm="Lloyd")
	names_w_10 <- c("Low risk","High blood pressure","Low DBP,low eGFR", "High heart rate","Severe hyperglycemia","Low BMI,high HDL","Mid-risk tall","Severe obesity","High cholesterol", "Mid-risk short")
	dat <- cbind(dat, cluster_10=km.10$cluster)
	dat$cluster_10 <- factor(dat$cluster_10, 1:10, names_w_10)
	dat$cluster_10 <- factor(dat$cluster_10, levels=c("Low risk","Mid-risk short","Mid-risk tall","Low BMI,high HDL","High heart rate","High cholesterol","High blood pressure", "Severe obesity","Severe hyperglycemia","Low DBP,low eGFR"))
	plot_clusters(dat, var_plot, var_label, var_reverse, clusters="cluster_10")

km_res_list <- list()
for (k in 5:12) {
	km_res_list[[k]] <- kmeans(dat.s, k, iter.max = 10000, nstart = 50, algorithm = "Lloyd")
	cluster_column_name <- paste("cluster", k, sep = "_")
	dat[[cluster_column_name]] <- km_res_list[[k]]$cluster
	dat[[cluster_column_name]] <- factor(dat[[cluster_column_name]], levels = 1:k, labels = arr_labels[[k - 4]])
}
saveRDS(dat, "men.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 3. Adjust sample weight 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat %>% mutate(midage = floor(age/5) * 5 + 2.5, midage=ifelse(age >= 80, 82.5, midage), age_grp = factor(midage), mid_year = as.factor(mid_year), samplewt_hba1c_adj = samplewt_hba1c)
	ss_year <- dat %>% group_by(mid_year) %>% summarise(n1=sum(samplewt_hba1c))
	tot <- sum(ss_year$n1); ss_year$adj <- (tot/ss_year$n1)/nrow(ss_year) 
	for (year in c(1991, 2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014, 2016, 2018)) {
		dat[dat$mid_year == year, "samplewt_year_adj"] <- dat[dat$mid_year == year, "samplewt_hba1c"] * as.numeric(ss_year[ss_year$mid_year == year, "adj"])
	}
	pop <- read.xlsx("D:/data/ukb/phe/Rdata/US_census_2020_formatted.xlsx") # from https://data.census.gov/
	pop <- pop[(pop$AgeGrp >= 20) & (pop$AgeGrp <= 90),]; pop$midage <- floor(pop$AgeGrp/5) * 5 + 2.5; pop$midage[pop$AgeGrp >= 80] <- 82.5
adjust <- function(data, group_var) {
	adjust <- list()
	data[,group_var] <- as.factor(data[,group_var])
	k <- 0
	for (group in levels(data[,group_var])) {
    d <- data %>% filter((!!sym(group_var)) == group)
    d_age <- d %>% count(midage,wt=samplewt_hba1c) %>% mutate(per=n/sum(n))
    true_pop_age <- pop %>% dplyr::select(midage,PopTotal) %>% count(midage,wt=PopTotal) %>% mutate(per=n/sum(n))
    d_age$adjust <- true_pop_age[true_pop_age$midage %in% d_age$midage,]$per/d_age$per
    adjust[[group]] <- d_age
	}
	return(adjust)
}
adjust_smplwt <- function(midage, midyear, samplewt_hba1c) {
	adj <- ratios[[midyear]][ratios[[midyear]]$midage==midage,"adjust"]; return(samplewt_hba1c * adj)
}
ratios <- adjust(data=dat, group_var="mid_year")
dat$samplewt_year_age_adj <- as.vector(mapply(adjust_smplwt, dat$midage, dat$mid_year, dat$samplewt_year_adj))
dat$mid_year <- as.numeric(as.character(dat$mid_year))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 4. Trends significance 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cluster="cluster_10"
clusters <- levels(dat[,cluster])
res <- data.frame(beta_year=rep(-1,length(clusters)), p_val_year=rep(-1,length(clusters)))
	row.names(res)<-clusters
	for (clust in clusters){
		dat$y <- ifelse(dat$cluster_10==clust,1,0)
		design <- svydesign(ids=dat$psu, strata=dat$stratum, data=dat, weights=~samplewt_year_age_adj, nest=TRUE)   
		mylogit <- svyglm(y ~ age_grp + mid_year, data=dat, design=design, family=quasibinomial(link="logit"))
		a<-summary(mylogit)
		res[clust,"beta_year"] <- a$coefficients["mid_year","Estimate"]
		res[clust,"p_val_year"] <- a$coefficients["mid_year","Pr(>|t|)"]
  }
write.csv(res, file="age_std_logit.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 5. Predict cluster membership
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- dat %>% mutate(smoker2 = ifelse(smoker==1,"Current","Never"), smoker2 =ifelse((smoke_ever==1)&(smoker==0),"Former",smoker2), smoker2 = factor(smoker2,levels= c("Never","Former","Current")))
clusters <- levels(dat[,"cluster_10"])
res <- data.frame(cluster=character(), var=character(), OR=integer(), p_val=integer())
for (c in clusters){
    dat$y <- ifelse(dat$cluster_10==c, 1, 0)
    design <- svydesign(ids=dat$psu, strata=dat$stratum, data=dat, weights=~samplewt_year_age_adj, nest=TRUE)   
	covs=c("age_grp ethnicity education drug_hyper drug_diab_insu drug_diab_pill drug_chol smoker2 self_mi self_stroke self_heart_failure") %>% strsplit(" ") %>% unlist()
	setdiff(covs, names(design$variables))
	mylogit <- svyglm(as.formula(paste('y ~ ', paste(covs, collapse='+'))), data=dat, design=design, family=quasibinomial(link='logit'))   
	a <- data.frame(summary(mylogit)$coefficients)
    a$cluster <- c
    a$p_val <- a$Pr...t..
    a$OR <- exp(a$Estimate)
    a$var <- row.names(a)
    res <- rbind(res,a[,names(res)])
}
write.csv(res, file="predictors.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 图2 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- dat
source(paste0(dir0, "/scripts/f/cmd.border_themes.R"))
p_val <- read.csv("age_std_logit.csv",row.names=c(1))
colors <- c("Low risk"="#B1CD7F", "Mid-risk short"="#A1C3AB", "Mid-risk tall"="gray83", "Low BMI,high HDL"="#DACCE6", "High cholesterol"= "#9C6FAE", "High heart rate"="#E4AFAA", "High blood pressure"="#004B96", "Severe obesity"="#EE7651", "Severe hyperglycemia"="#B20832", "Low DBP,low eGFR"="#ACD9EA")
s_letter=16 
s_numbers= 6 
s_axis_1=18 
s_side_title=4 
s_leg=16
p_val$P_val_formated <- formatC(p_val$p_val_year,format="e",digits=2) 
dat1$mid_year <- ifelse(dat1$mid_year==1991,1992,floor((dat1$mid_year-2000)/4)*4+2000) 
dat1$mid_year <- as.factor(dat1$mid_year)
d0 <- dat1 %>% dplyr::select(mid_year,cluster_10,samplewt_year_age_adj)
	d1 <- d0 %>% group_by(mid_year) %>% summarise(n=sum(samplewt_year_age_adj))
	d2 <- d0 %>% group_by(mid_year,cluster_10) %>% summarise(n1=sum(samplewt_year_age_adj)) %>% left_join(d1) %>% mutate(freq=n1/n) %>%  dplyr::select(-n1,-n)
	d2$cluster_10 <- factor(d2$cluster_10)
	d2$freq <- round(d2$freq,3)
	d2$percent <- d2$freq*100
	d2$percent <- sprintf("%.1f",d2$percent) 
	d2$mid_year <- c(rep(1991,10),rep(2000,10),rep(2004,10),rep(2008,10),rep(2012,10),rep(2016,10))

p1 <- ggplot(d2,aes(x=mid_year,y=freq,fill=cluster_10,label=percent)) + 
  geom_col(position="fill")+ scale_fill_manual(values=colors,guide=guide_legend(reverse=FALSE)) +
  labs(y="Percent of population (%)",x="Year",fill="Cluster")+ scale_y_continuous(labels=scales::percent,expand=c(0,0))+theme(plot.title=element_text(hjust=0.5,size=15),plot.subtitle=element_text(hjust=0.5,size=10))+ theme_bw()+
  theme(axis.text.x=element_text(size=s_axis_1), axis.title.x=element_text(size=s_axis_1), axis.text.y=element_text(size=s_axis_1), axis.title.y=element_text(size=s_axis_1),legend.title=element_text(size=s_leg), legend.text=element_text(size=s_leg), panel.grid.minor=element_blank(), panel.grid.major=element_blank())+
  theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.position="right")+ geom_text(size=s_numbers,position=position_stack(vjust=0.5),colour="white") +theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(),legend.position="None") 

cluster_10.labs <- c()
for (cluster in (row.names(p_val))){
  cluster_10.labs <- c(cluster_10.labs,paste0(cluster," P for trend=(",p_val[cluster,"P_val_formated"],")"))
}
names(cluster_10.labs) <- row.names(p_val)

p <- list()
i <- 1
clusters <- levels(d2$cluster_10)
d2$percent <- d2$freq*100
bwn_plots <- -0.4 
for (cluster in clusters[1:(length(clusters)-1)]){ 
	if(cluster=="Low risk"){margin=c(0,0,bwn_plots,0); border <- c("top","bottom","left","right")}else{ margin=c(bwn_plots,0,bwn_plots,0); border <- c("bottom","left","right")}
	col <- colors[[cluster]]
	data <- d2[d2$cluster_10==cluster,]
	title <- unname(cluster_10.labs[cluster]) 
	M <- max(data$percent)*1.25
	p4 <- ggplot(data,aes(x=mid_year,y=percent,group=1)) +  geom_line(size=2,colour=col)+ labs(y="",x="")
	if(cluster!="Severe hyperglycemia"){ p4 <- p4+scale_y_continuous(limits=c(0,M),breaks=c(0,5,10,15,20),expand=c(0,1))} else {p4 <- p4+scale_y_continuous(limits=c(0,M),breaks=c(0,1,2,3,4),expand=c(0,0.25))}
	p4 <- p4+scale_color_manual(values=colors) + guides(colour=guide_legend(title="",title.theme=element_text(size=12, face="plain",colour="black", angle=0)))+
    theme_bw()+ theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), panel.grid.minor=element_blank(), strip.background=element_blank(), strip.text=element_blank(),panel.grid.major=element_blank(), panel.border=theme_border(type=border,size=0.5), plot.margin= unit(margin,"lines"))+annotate("text",label=title,x=2003,y=M*0.95,size=s_side_title)
	p[[i]] <- p4
	i <- i+1
}
cluster="Low DBP,low eGFR"
col <- colors[[cluster]]
data <- d2[d2$cluster_10==cluster,]
border <- c("bottom","left","right")
title <- unname(cluster_10.labs[cluster])
M <- max(data$percent)*1.25 
p2 <- ggplot(data,aes(x=mid_year,y=percent,group=1)) + geom_line(size=2,colour=col)+ labs(y="",x="Year")+
	scale_y_continuous(limits=c(0,M),breaks=c(0,5,10,15,20),expand=c(0,1))+ scale_color_manual(values=colors) + guides(colour=guide_legend(title="",title.theme=element_text(size=12, face="plain", colour="black", angle=0)))+
	theme_bw()+ theme(axis.ticks.x=element_blank(), panel.grid.minor=element_blank(), strip.background=element_blank(), strip.text=element_blank(), panel.grid.major=element_blank(), panel.border=theme_border(type=border,size=0.5), plot.margin= unit(c(bwn_plots,0,-1.7,0),"lines") )+annotate("text",label=title,x=2003,y=M*0.95,size=s_side_title)
p[[i]] <- p2
side <- egg::ggarrange(plots=p,ncol=1)
a <- annotate_figure(side, left=textGrob("Percent of population (%)",rot=90,vjust=0.5,gp=gpar(cex=1.3))) 
p <- list()
p[[1]] <- p1
p[[2]] <- a
plot1 <- ggarrange(plots=p, nrow=1, legend="none", widths=c(3,1))
plot1 <- annotate_figure(plot1,fig.lab="",fig.lab.size=s_letter)
plot1 <- plot1+plot_annotation(title="",theme=theme(plot.title=element_text(size=20)))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 图3 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- dat
dat1$midage <- floor(dat1$age/5)*5 + 2.5
dat1$midage[dat1$age>=80]<-82.5 
dat1$mid_year <- ifelse(dat1$mid_year==1991,1991,floor((dat1$mid_year-2000)/4)*4+2000)
dat1$cohort <- dat1$mid_year - dat1$midage
d1 <- dat1 %>% dplyr::select(mid_year,midage,cluster_10,samplewt_hba1c) %>% group_by(mid_year,midage) %>% summarise(n=sum(samplewt_hba1c))
d01 <- dat1 %>% dplyr::select(mid_year,midage,cluster_10,samplewt_hba1c) %>% group_by(mid_year,midage,cluster_10) %>% summarise(n1=sum(samplewt_hba1c)) %>% left_join(d1) %>% mutate(freq=n1/n) %>% ungroup() %>%dplyr::select(-n1,-n)
	d01$mid_year<-as_factor(d01$mid_year)
colors_years <- c("1991"="saddlebrown", "2000"="darkorange3", "2004"="goldenrod", "2008"="lightgoldenrod3", "2012"="gray75", "2016"="gray53")
p1 <- ggplot(d01,aes(midage,freq,group=paste(cluster_10,mid_year),colour=mid_year)) + theme_bw() + geom_point(shape=16,size=2) + labs(y="Prevalence (%)",x="Age")+ scale_y_continuous(labels=scales::percent) + scale_colour_manual(values=colors_years)+
	guides(colour=guide_legend(title="Year", title.theme=element_text( size=14, face="plain", colour="black", angle=0)))+
	theme(axis.text.x=element_text(color="black",size=12), axis.title.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=12),text=element_text(family="sans"), legend.title=element_text(size=12), legend.position="bottom", legend.text=element_text(size=12), plot.title=element_text(hjust=0.5,size=14,face="bold"), strip.text=element_text(size=12))+ facet_wrap(~cluster_10,nrow=2)
p1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 图4-5 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- read_csv(file="predictors.csv"); dat1$...1 <- NULL
reverselog_trans <- function(base=exp(1)) { # for the y-scale of the plot
  trans <- function(x) -log(x,base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-",format(base)),trans,inv, log_breaks(base=base), domain=c(1e-120,Inf))
}
library(scales); volcano_s(data=dat1, threshold=0.05, title="Men")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 图6 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(networkD3); library(ggplot2)

cls_map <- c(
  "Low risk" = "#B1CD7F", "Severe hyperglycemia" = "#B20832", "Overweight" = "#E4AFAA", "Severe obesity" = "#EE7651", "Mid-risk short" = "#A1C3AB", 
  "Mid-risk" = "#EEC591", "Mid-risk tall" = "#828282", "Low DBP, low eGFR" = "#ACD9EA", "Low BMI, high HDL" = "#DACCE6",
  "High blood pressure" = "#004B96", "High cholesterol" = "#9C6FAE", "High heart rate" = "#8B3E2F", "High SBP, low eGFR" = "#609CBF"
)

dat1 <- readRDS("men.rds")
dat1 <- dat1[, c("cluster_5", "cluster_6", "cluster_7", "cluster_8", "cluster_9", "cluster_10", "cluster_11", "cluster_12")]
links <- data.frame(
  source = rep(1:nrow(dat1), each = ncol(dat1)),
  target = rep(1:ncol(dat1), times = nrow(dat1)),
  value = as.vector(as.matrix(dat1))
)

sankey_plot <- sankeyNetwork(
  Links = links,
  Source = "source",
  Target = "target",
  Value = "value",
  NodeID = "target",  # Or any column that contains node names
  fontSize = 15,
  nodeWidth = 30,
  colourScale = cls_map
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 附图5 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
variables <- c("height","bmi","whtr","hba1c","hdl","non_hdl","sbp","dbp","eGFR","pulse")
dat1 = dat
	inter_intra <- cls.scatt.data(data=scale(dat1[,variables]),clust=as.numeric(dat1[,"cluster_10"]))
	a <- inter_intra$cluster.center
	dist <- inter_intra$intercls.average
	a <- inter_intra$intracls.average
	intra <- as.numeric(a)
	diag(dist) <- intra
	melted_res <- melt(dist)
	remove <- c(11,21,22,seq(31,33),seq(41,44),seq(51,55),seq(61,66),seq(71,77),seq(81,88),seq(91,99))
	melted_res <- melted_res[!row.names(melted_res)%in%remove,]
	max <- 7
	min <- 2.6
	names_cluster <- levels(dat1$cluster_10)
	melted_res$Var1 <- factor(melted_res$Var1,c("c1","c2","c3","c4","c5","c6","c7","c8","c9","c10"),names_cluster)
	melted_res$Var2 <- factor(melted_res$Var2,c("c1","c2","c3","c4","c5","c6","c7","c8","c9","c10"),names_cluster)
p1 <- ggplot(data=melted_res,aes(x=Var1,y=Var2,fill=value)) + geom_tile()+
	scale_fill_gradient2(low="white",high="blue",mid="lightblue", midpoint=1.0*(max+min)/2,limit=c(min,max),space="Lab", name="Average\ndistance")+
	geom_text(aes(Var1,Var2,label=round(value,2)),color="black",size=2.9)+
	xlab("")+ ylab("")+ theme_minimal()+ 
	theme(axis.text.x=element_text(angle=45,vjust=1, size=10,hjust=1), axis.text.y=element_text(size=10,hjust=1), legend.position="right", text=element_text(family="sans"),legend.key.size=unit(1,'cm'), legend.text=element_text(size=10),plot.title=element_text(hjust=0.5,size=20),legend.title=element_text(size=10), panel.grid=element_blank()) + plot_annotation(title="Men", theme=theme(plot.title=element_text(size=16)))
p1