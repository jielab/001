pacman::p_load(tidyverse, MASS, ggpubr, reshape2, survey, grid, egg, patchwork)

dir0='D:'
set.seed(1234)

plot_clusters <- function(data, var_plot, var_reverse, var_10levles, clusters, title="", legend="bottom", add_values=FALSE){
	cluster_names <- levels(as.factor(data[,clusters]))
	plot_df <- data %>% group_by_at(clusters) %>% summarize_at(.vars=var_plot,.funs=median) %>% as.data.frame(plot_df)
	quantile_cluster <- plot_df
	row.names(quantile_cluster) <- cluster_names
	for (i in (cluster_names)){
	for (var in var_plot){
		percentile <- ecdf(data[,var])
		quantile_cluster[which(quantile_cluster[,clusters]==i),var] <- round(percentile(plot_df[which(plot_df[,clusters]==i),var]),2)
	}
	}
	df <- melt(quantile_cluster, id.vars=c(clusters))
	a  <- as.data.frame(plot_df[,var_plot])
	row.names(a) <- cluster_names
	for (var in var_reverse) {df[which(df$variable == var),]$value <- 1 - df[which(df$variable == var),]$value}
	filling="quantile"
	df$quantile <- df$value*100 
	df$variable <- factor(df$variable, levels=var_10levles)
	p <- list()
	for (i in cluster_names){
		b <- as.matrix(a[i,var_plot])
		val <- c(round(b[1],0),round(b[2],1),round(b[3],2),round(b[4],1),round(b[5],1),round(b[6],1),round(b[7],1),round(b[8],1),round(b[9],0),round(b[10],0))
		start_angle <- pi/10 
		p1 <- ggplot(df[which(df[,clusters]==i),],aes_string(x="variable",y="value",fill=filling)) +
		geom_col(width=0.5)+
		geom_hline(yintercept=c(0,0.25,0.75,1),colour="grey44",size=0.1,alpha=0.4,linetype="solid")+        
		geom_hline(yintercept=c(0.5),colour="black",size=0.3,alpha=0.4,linetype="solid") +
		scale_y_continuous(limits=c(-0.2,1.1), expand=c(0,0), breaks=c(0,0.25,0.5,0.75,1),position="right")	
		p1 <- p1+scale_fill_gradientn(name=stringr::str_wrap("Risk factor percentile"), colours=c("forestgreen","chartreuse","gold","orange","red"), breaks=c(10,30,50,70,90),limits=c(0,100), labels=c("more optimal  ","","","","  more hazardous"), guide=guide_colourbar(title.position="left",label.hjust=0.5,title.hjust=30))
		if(add_values){ p1 <- p1+ geom_text(label=val,vjust=-0.25,size=4.8)} 
		p1 <- p1 + theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(color="black",size=12), legend.text=element_text(size=10),
        text=element_text(family="sans"), panel.background=element_rect(fill="white",color="white"), plot.title=element_text(hjust=0.5,size=14,face="bold"), plot.subtitle=element_text( hjust=0.5))+ 
		labs(title=paste0(i," (",round((table(data[,clusters])/length(data[,clusters]))[i]*100,0),"%)"))+
		coord_polar(theta='x',start=start_angle,clip='off') +
		geom_segment(aes(x=seq(1,10),xend=seq(1,10),y=0,yend=1),linetype="dotted",size=0.2)
		p[[i]] <- p1
	}
	if((length(cluster_names))<13){ n_row <- 3 }else{ n_row <- (length(cluster_names)%/%4)+1}
	q <- ggpubr::ggarrange(plotlist=p, common.legend=TRUE,ncol=4,nrow=n_row, font.label=list(color="black"),legend=legend)
	q <- q+plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5)))
	return(q)
}


var_include <- c("height", "bmi", "whtr", "hba1c", "hdl", "non_hdl", "sbp", "dbp", "eGFR", "pulse")
var_reverse <- c('height', 'hdl', 'eGFR')
var_labels <- list( # cluster 5-12
	c("Mid-risk", "Low risk", "Severe hyperglycemia", "High blood pressure", "Severe obesity"),
	c("Low risk", "High blood pressure", "Low DBP, low eGFR", "Mid-risk", "Severe obesity", "Severe hyperglycemia"),
	c("Low DBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk", "Severe hyperglycemia", "Low risk", "Low BMI, high HDL"),
	c("Mid-risk short", "Mid-risk tall", "High blood pressure", "Severe obesity", "Severe hyperglycemia", "Low DBP, low eGFR", "Low risk", "Low BMI, high HDL"),
	c("Mid-risk tall", "Low risk", "High blood pressure", "Mid-risk short", "High cholesterol", "Severe hyperglycemia", "Low BMI, high HDL", "Low DBP, low eGFR", "Severe obesity"),
	c("Low risk", "High blood pressure", "Low DBP,low eGFR", "High heart rate", "Severe hyperglycemia", "Low BMI,high HDL", "Mid-risk tall", "Severe obesity", "High cholesterol", "Mid-risk short"),
	c("Low risk", "Low DBP, low eGFR", "High heart rate", "Mid-risk tall", "Low BMI, high HDL", "Severe hyperglycemia", "High blood pressure", "Mid-risk short", "Severe obesity", "Overweight", "High cholesterol"),
	c("Mid-risk tall", "Overweight", "High cholesterol", "Low risk", "High SBP, low eGFR", "High blood pressure", "Severe obesity", "Mid-risk short", "Low DBP, low eGFR", "High heart rate", "Low BMI, high HDL", "Severe hyperglycemia")
)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1. UKB cluster 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS(file="../data/all.rds")
dat <- dat0 %>% filter(ethnic_cat=="White", !is.na(icdDate_is.2), complete.cases(across(all_of(var_include)))) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# NHANES cluster 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat0 <- readRDS(file=paste0(dir0, '/data/ukb/phe/Rdata/NHANES_Cleaned_single.RDS')) %>% 
	rename(sbp=sbp_final, dbp=dbp_final) %>% filter(age>=20)
	vars0 <- grep("clean", names(dat0), value=T); vars1 <- gsub("_clean$", "", vars0)
	dat0[,vars1] <- dat0[,vars0]; dat0[,vars0] <- NULL
	
dat <- dat0 %>% filter(sex==1) %>% mutate(non_hdl = tc - hdl, whtr = waist/height) %>% drop_na(all_of(var_include)) 
	dat.s <- scale(dat[,var_include])

for (k in 5:12) {
	km_res <- kmeans(dat.s, k, iter.max=10000, nstart=50, algorithm = "Lloyd")
	cluster_col <- paste("cluster", k, sep = "_")
	dat[[cluster_col]] <- km_res$cluster
	dat[[cluster_col]] <- factor(dat[[cluster_col]], levels=1:k, labels=var_labels[[k - 4]])
}
plot_clusters(dat, var_include, var_reverse, clusters="cluster_8")
saveRDS(dat, "dat.rds")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Predict cluster membership
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat <- readRDS("dat.rds") %>% mutate(
	smoker2 = ifelse(smoker==1,"Current","Never"), smoker2 =ifelse((smoke_ever==1)&(smoker==0),"Former",smoker2), smoker2 = factor(smoker2,levels= c("Never","Former","Current"))
)
clusters <- levels(dat[,"cluster_10"])
res <- data.frame(cluster=character(), var=character(), OR=integer(), p_val=integer())
for (c in clusters){
    dat$y <- ifelse(dat$cluster_10==c, 1, 0)
    design <- svydesign(ids=dat$psu, strata=dat$stratum, data=dat, weights=~age, nest=TRUE) # 🏮
	covs=c("ethnicity education drug_hyper drug_diab_insu drug_diab_pill drug_chol smoker2 self_mi self_stroke self_heart_failure") %>% strsplit(" ") %>% unlist()
	setdiff(covs, names(design$variables)) # 必须显示没有差别
	mylogit <- svyglm(as.formula(paste('y ~ ', paste(covs, collapse='+'))), data=dat, design=design, family=quasibinomial(link='logit'))   
	a <- data.frame(summary(mylogit)$coef)
    a$cluster <- c
    a$p_val <- a$Pr...t..
    a$OR <- exp(a$Estimate)
    a$var <- row.names(a)
    res <- rbind(res,a[,names(res)])
}
write.csv(res, file="predictors.csv")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 附图5 🎇
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(clv)
	inter_intra <- clv::cls.scatt.data(data=scale(dat[,var_include]),clust=as.numeric(dat[,"cluster_10"]))
	a <- inter_intra$cluster.center
	dist <- inter_intra$intercls.average
	a <- inter_intra$intracls.average
	intra <- as.numeric(a)
	diag(dist) <- intra
	melted_res <- melt(dist)
	remove <- c(11,21,22,seq(31,33),seq(41,44),seq(51,55),seq(61,66),seq(71,77),seq(81,88),seq(91,99))
	melted_res <- melted_res[!row.names(melted_res)%in%remove,]
	max <- 7; min <- 2.6
	names_cluster <- levels(dat$cluster_10)
	melted_res$Var1 <- factor(melted_res$Var1, paste0("c", 1:10), names_cluster)
	melted_res$Var2 <- factor(melted_res$Var2, paste0("c", 1:10), names_cluster)
ggplot(data=melted_res, aes(x=Var1,y=Var2,fill=value)) + geom_tile() +
	scale_fill_gradient2(low="white",high="blue",mid="lightblue", midpoint=1.0*(max+min)/2,limit=c(min,max),space="Lab", name="Average\ndistance")+
	geom_text(aes(Var1,Var2,label=round(value,2)),color="black",size=2.9)+
	xlab("")+ ylab("")+ theme_minimal()+ 
	theme(axis.text.x=element_text(angle=45,vjust=1, size=10,hjust=1), axis.text.y=element_text(size=10,hjust=1), legend.position="right", text=element_text(family="sans"),legend.key.size=unit(1,'cm'), legend.text=element_text(size=10),plot.title=element_text(hjust=0.5,size=20),legend.title=element_text(size=10), panel.grid=element_blank()) + 
	plot_annotation(title="Men", theme=theme(plot.title=element_text(size=16)))