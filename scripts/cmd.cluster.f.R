
plot_clusters <- function(data, var_plot, var_label, var_reverse, clusters, title="", legend="bottom", add_values=FALSE){
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
	a <-as.data.frame(plot_df[,var_plot])
	row.names(a) <- cluster_names
  
	# reversing axis for protective Risk factors
	for (var in var_reverse) {df[df$variable == var, "value"] <- 1 - df[df$variable == var, "value"]}
	filling="quantile"
	df$quantile <- df$value*100 
  
	levels(df$variable) <- var_plot
	df$variable <- factor(df$variable,levels=var_label)
  
	p <- list()
	for (i in cluster_names){
    b <- as.matrix(a[i,var_plot])
    val <- c(round(b[1],0),round(b[2],1),round(b[3],2),round(b[4],1),round(b[5],1),round(b[6],1),round(b[7],1),round(b[8],1),round(b[9],0),round(b[10],0))
    start_angle <- pi/10 
    p1 <- ggplot(df[which(df[,clusters]==i),],aes_string(x="variable",y="value",fill=filling)) +
		geom_col(width=0.5)+
		geom_hline(yintercept=c(0,0.25,0.75,1),colour="grey44",size=0.1,alpha=0.4,linetype="solid") +        
		geom_hline(yintercept=c(0.5),colour="black",size=0.3,alpha=0.4,linetype="solid") +
		scale_y_continuous(limits=c(-0.2,1.1), expand=c(0,0), breaks=c(0,0.25,0.5,0.75,1),position="right") +
		scale_fill_gradientn(name=stringr::str_wrap("Risk factor percentile"), colours=c("forestgreen","chartreuse","gold","orange","red"), breaks=c(10,30,50,70,90),limits=c(0,100), labels=c("more optimal  ","","","","  more hazardous"), guide=guide_colourbar(title.position="left",label.hjust=0.5,title.hjust=30)) +
		theme(axis.title=element_blank(), axis.ticks=element_blank(), axis.text.y=element_blank(), axis.text.x=element_text(color="black",size=12), legend.text=element_text(size=10), text=element_text(family="sans"), panel.background=element_rect(fill="white",color="white"), plot.title=element_text(hjust=0.5,size=14,face="bold"), plot.subtitle=element_text(hjust=0.5))+ 
		labs(title=paste0(i," (",round((table(data[,clusters])/length(data[,clusters]))[i]*100,0),"%)") ) +
		coord_polar(theta='x',start=start_angle,clip='off') +
		geom_segment(aes(x=seq(1,10),xend=seq(1,10),y=0,yend=1),linetype="dotted",size=0.2)
		p[[i]] <- p1
	}
  
	if((length(cluster_names))<13){n_row <- 3} else {n_row <- (length(cluster_names)%/%4)+1}
	q <- ggpubr::ggarrange(plotlist=p, common.legend=TRUE, ncol=4,nrow=n_row, font.label=list(color="black"), legend=legend)
	q <- q + plot_annotation(title=title, theme=theme(plot.title=element_text(hjust=0.5)))
	return(q)
}

volcano_s <- function(data,threshold=0.05,title="A)"){
  names(data) <- c("cluster","var","mean","p")
  # not plotting the intercept
  data <- data[!grepl("(Intercept)",data$var),]
  # Add labels for the plots
  data$label <- ""
  data$label <- ifelse(data$p<=threshold,data$var,"")
  data[grepl("age_grp27.5",data$label),"label"] <- "age[25,30)"
  data[grepl("age_grp32.5",data$label),"label"] <- "age[30,35)"
  data[grepl("age_grp37.5",data$label),"label"] <- "age[35,40)"
  data[grepl("age_grp42.5",data$label),"label"] <- "age[40,45)"
  data[grepl("age_grp47.5",data$label),"label"] <- "age[45,50)"
  data[grepl("age_grp52.5",data$label),"label"] <- "age[50,55)"
  data[grepl("age_grp57.5",data$label),"label"] <- "age[55,60)"
  data[grepl("age_grp62.5",data$label),"label"] <- "age[60,65)"
  data[grepl("age_grp67.5",data$label),"label"] <- "age[65,70)"
  data[grepl("age_grp72.5",data$label),"label"] <- "age[70,75)"
  data[grepl("age_grp77.5",data$label),"label"] <- "age[75,80)"
  data[grepl("age_grp82.5",data$label),"label"] <- "age80+" 
  data[grepl("smoker2Former",data$label),"label"] <- "Smoker (former)"
  data[grepl("smoker2Current",data$label),"label"] <- "Smoker (current)"
  data[grepl("ethnicityOther",data$label),"label"] <- "Other ethnicity"
  data[grepl("ethnicityNon-Hispanic Black",data$label),"label"] <- "non-Hispanic Black"
  data[grepl("ethnicityHispanic",data$label),"label"] <- "Hispanic"
  data[grepl("educationHigh School",data$label),"label"] <- "High school"
  data[grepl("educationUniversity or College",data$label),"label"] <- "University or college"
  data[grepl("drug_hyper",data$label),"label"] <- "Antihypertensive"
  data[grepl("drug_chol",data$label),"label"] <- "Statins"
  data[grepl("drug_diab_pill",data$label),"label"] <- "Oral hypoglycemic"
  data[grepl("drug_diab_insu",data$label),"label"] <- "Insulin"
  data[grepl("mid_year",data$label),"label"] <- "Year"
  data[grepl("self_mi",data$label),"label"] <- "MI history"
  data[grepl("self_stroke",data$label),"label"] <- "Stroke history"
  data[grepl("self_heart_failure",data$label),"label"] <- "CHF history"
  data$color <- "P>=0.05"
  data$color <- ifelse((data$p<threshold),"P<0.05",data$color)
  cols <- c("P<0.05"="black","P>=0.05"="grey")

  res <- list()
  for (clust in unique(data$cluster)){
    d <- data[data$cluster==clust,]
    p_min <- min(d$p)
    m <- min(d$mean)
    s <- max(d$mean)-min(d$mean)
    center <- m+s*0.9
    volcano <- ggplot(d,aes(x=mean,y=p,color=color,label=label)) +
      geom_point() + 
      labs(title =clust)+
      scale_color_manual(values=cols,labels =c(expression("P"<"0.05"),expression("P">="0.05"))) +
      geom_hline(aes(yintercept=threshold,linetype="Custom Threshold"),col="cadetblue3",size=0.5) +
      scale_linetype_manual(
        name="Threshold",
        values=c("Custom Threshold"=2), 
        labels=c(paste0("P=",round(threshold,5))) 
      ) +
      guides(
        color=guide_legend(title="Odds ratios",title.position="top",title.hjust=0.5,override.aes=aes(label="")),
        linetype=guide_legend(title="",title.position="top",title.hjust=0.5)
      ) +
      scale_y_continuous(trans=reverselog_trans(10))+
      scale_x_log10()+
      xlab("Odds ratio relative to the whole population")+
      ylab("P value")+
      geom_text_repel(size=4,max.overlaps=100) +
      theme(
        axis.line=element_line(size=0.5,colour="black",linetype=1),
        axis.text.x=element_text(color="black",size=12),
        axis.title=element_text(color="black",size=12),
        legend.position="right",
        legend.key.size=unit(1.8,'cm'),
        legend.background= element_rect(fill="white",color="white"),
        legend.key=element_rect(fill="white"),
        legend.title.align=0.5,
        panel.background=element_rect(fill="white",color="white"),
        plot.title=element_text(hjust=0.5,size=14,face="bold"))
    res[[clust]] <- volcano
  }
  leg <- get_legend(res[["Low risk"]])
  res[["legend"]] <- as_ggplot(leg)
  res_all <- ggpubr::ggarrange(plotlist=res, common.legend=TRUE, ncol=4,nrow=3,font.label=list(size=20,color="black"), legend="none")
  res_all <- res_all+plot_annotation(
    title=title,
    theme=theme(plot.title=element_text(size=25,face="plain",vjust=0.5)))
  return(res_all)
}

alpha_q2 <- function(q){ # coloring for extended data figure 1.
  if((q>45)&(q<55)){
    res <- 1
  }else if((q>25)&(q<75)){
    res <- 1-(abs(50-q))/29
  }else{
    if(q<=25){
      res <- 0.1724138-min(((25-q)/100),0.1724138)
    }else(
      res <- 0.1724138-min(((q-75)/100),0.1724138)
    )
  }
  return(res)
}

plot_shaded_percentile <- function(data,color_bar='blue',title="",size_bar=1.2){
  cluster_names <- levels(data[,"cluster_10"])
  res <- list()
  for (clust in cluster_names){
    d <- data %>% 
      subset(cluster_10==clust)
    var_plot <- c("height","bmi","whtr","hba1c","hdl","non_hdl","sbp","dbp","eGFR","pulse")
    p <- d%>%
      group_by(cluster_10)
    df <- data.frame(variables=factor(),ymin=double(),ymax=double(),percentile=double())
    for(percent in seq(0,99)){
      for (var in var_plot){
        q <- percent/100
        percentile <- ecdf(data[,var])
        min <- quantile(d[,var],q)
        max <- quantile(d[,var],q+0.01)
        dtemp <- data.frame(variables=var,ymin=percentile(min)*100,ymax=percentile(max)*100,percentile=q)
        df <- rbind(df,dtemp)
      }
    }

    df$variables <- as.factor(df$variables)
    levels(df$variables) <- c("BMI","DBP","(low) \n eGFR","HbA1c","(low) \n HDL","(low) \n Height","Non-HDL","Pulse \n rate","SBP","WHtR")   
    df$variables <- factor(df$variables,levels=c("(low) \n Height","BMI","WHtR","HbA1c","(low) \n HDL","Non-HDL","SBP","DBP","Pulse \n rate","(low) \n eGFR"))
    temp <- df[df$variables=="(low) \n Height",'ymin']
    df[df$variables=="(low) \n Height",'ymin'] <- 100-df[df$variables=="(low) \n Height",'ymax']
    df[df$variables=="(low) \n Height",'ymax'] <- 100-temp   
    temp <- df[df$variables=="(low) \n HDL",'ymin']
    df[df$variables=="(low) \n HDL",'ymin'] <- 100-df[df$variables=="(low) \n HDL",'ymax']
    df[df$variables=="(low) \n HDL",'ymax'] <- 100-temp   
    temp <- df[df$variables=="(low) \n eGFR",'ymin']
    df[df$variables=="(low) \n eGFR",'ymin'] <- 100-df[df$variables=="(low) \n eGFR",'ymax']
    df[df$variables=="(low) \n eGFR",'ymax'] <- 100-temp
    
    p <- ggplot(df,aes(x=variables,y=percentile,group=1)) +
      geom_linerange(d=df[df$percentile==0,] %>% group_by(variables),aes(ymin=ymin,ymax=ymax),color="navy", alpha=0.001**1.5, size=size_bar) 
    
    for(percent in seq(1,99)){
      p <- p+geom_linerange(d=df[df$percentile==percent/100,]%>% group_by(variables),aes(ymin=ymin,ymax=ymax),color="navy", alpha=alpha_q2(percent), size=size_bar) 
    }
    p <- p+geom_hline(yintercept=c(0,25,75,100),colour="grey44",size=0.1,alpha=0.4,linetype="solid")+       
      geom_hline(yintercept=c(50),colour="black",size=0.3,alpha=0.4,linetype="solid") +
      scale_y_continuous(
        limits=c(-20,110),
        expand=c(0,0),
        breaks=c(0,25,50,75,100),position="right"
      )+
      labs(title =clust)+
      theme(
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(color="black",size=10),
        legend.position="none",
        text=element_text(family="sans"),
        panel.background=element_rect(fill="white",color="white"),
        plot.title=element_text(hjust=0.5,size=14,face="bold"),
        plot.subtitle=element_text( hjust=0.5)
      )+ 
      coord_polar(theta='x',start=pi/10,clip='off')
    res[[clust]] <- p
  }
  n_row <- 3 
  legend <- data.frame(percentile=(seq(0,99))/100,y="T")
  l <- ggplot(legend,aes(x=y,y=percentile))+
    geom_linerange(d=legend[legend$percentile==0,]%>%
                     group_by(y),aes(ymin=0,ymax=1),color="navy",alpha=0.001**1.5,
                   size=300)+theme_bw()
  for(percent in seq(1,99)){
    l <- l+geom_linerange(d=legend[legend$percentile==percent/100,]%>%
                          group_by(y),
                        aes(ymin=percentile,
                            ymax=percentile+0.01),
                        color="navy",alpha=alpha_q2(percent),
                        size=300) 
  }
  l <- l+xlab("")+
    coord_flip()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),
                       aspect.ratio=0.15,
                       text=element_text(family="sans"),
                       panel.background=element_blank(),
                       axis.line=element_line(colour="black"),
                       axis.text.y =element_blank(),
                       axis.text.x=element_text(color="black",size=10),
                       axis.title.x=element_text(color="black",size=10),
                       axis.ticks.y=element_blank(),
                       plot.margin=margin(0,0,0,0,"cm"))
  
  res[["legend"]] <- l
  
  q <- ggpubr::ggarrange(plotlist=res,
                       ncol=4,nrow=3,font.label=list(color="black"))
  
  q <- q+plot_annotation(
    title=title,
    theme=theme(
      plot.title=element_text(hjust=0.5)
    ))+
  return(q)
}

plot_circular_age_group <- function(data,clusters,title=""){
  cluster_names <- levels(data[,clusters])
  p <- list()
  for (cluster_name in cluster_names){
    
    data$age_group <- ifelse(data$age<40,"20-39",0)
    data$age_group <- ifelse((data$age>=40)&(data$age<60),"40-59",data$age_group)
    data$age_group <- ifelse((data$age>=40)&(data$age>=60),"60+",data$age_group)
    data$age_group <- as.factor(data$age_group)
    var_plot <- c("height","bmi","whtr","hba1c","hdl","non_hdl","sbp","dbp","eGFR","pulse")
    group_cols  <-  c(clusters,"age_group")
    
    plot_df <- data[which(data[,clusters]==cluster_name),] %>%
      group_by(across(all_of(group_cols))) %>% 
      summarize_at(.vars=var_plot,.funs=c(median))
    plot_df <- as.data.frame(plot_df)
    quantile_cluster <- plot_df
    
    for (group in c("20-39","40-59","60+")){
      for (var in var_plot){
        percentile <- ecdf(data[,var])
        quantile_cluster[which((quantile_cluster[,clusters]==cluster_name)&(quantile_cluster[,"age_group"]==group)),var] <- 
          round(percentile(plot_df[which((quantile_cluster[,clusters]==cluster_name)&(quantile_cluster[,"age_group"]==group)),var]),2)
        
      }
    }
    quantile_med <- quantile_cluster[,c(clusters,"age_group",var_plot)]
    
    df_med <- melt(quantile_med,id.vars=c(clusters,"age_group"))
    
    df_med[which(df_med$variable=='height'),]$value <- 1-df_med[which(df_med$variable=='height'),]$value
    df_med[which(df_med$variable=='hdl'),]$value <- 1-df_med[which(df_med$variable=='hdl'),]$value
    df_med[which(df_med$variable=='eGFR'),]$value <- 1-df_med[which(df_med$variable=='eGFR'),]$value
    df <- df_med
    colors <- c("60+"="navy","40-59"="orange","20-39"="darkgreen")
    sizes <- c("60+"=1,"40-59"=1,"40-59"=1)
    levels(df$variable) <- c("(low) \n Height","BMI","WHtR","HbA1c","(low) \n HDL","Non-HDL","SBP","DBP","(low) \n eGFR","Pulse \n rate")
    df$variable <- factor(df$variable,levels=c("(low) \n Height","BMI","WHtR","HbA1c","(low) \n HDL","Non-HDL","SBP","DBP","Pulse \n rate","(low) \n eGFR"))
    
    bridge <- df[which(df$variable=="(low) \n Height"),] # The "bridge" is to connect start and end of the circular plot when going polar
    bridge$variable <- NA 
    
    df$age_group <- factor(df$age_group,levels=c("20-39","40-59","60+"))
    p1 <- ggplot(rbind(df,bridge),aes(x=variable,y=value,group=age_group,colour=age_group,fill=age_group)) +
      geom_line(size=0.5)+
      scale_fill_manual(values=colors)+
      geom_hline(yintercept=c(0,0.25,0.75,1),colour="grey44",size=0.1,alpha=0.4,linetype="solid")+         
      geom_hline(yintercept=c(0.5),colour="black",size=0.3,alpha=0.4,linetype="solid") +
      geom_segment(aes(x=variable,xend=variable,y=0,yend=1),linetype="dotted",size=0.2)+
      scale_y_continuous(
        limits=c(-0.2,1.1),
        expand=c(0,0),
        breaks=c(0,0.25,0.5,0.75,1),position="right"
      )+scale_color_manual(name="Age groups \n (years)",values=colors)+
      scale_x_discrete(expand=c(0,0),breaks=levels(df$variable))+
      guides(size=FALSE)+
      theme(
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        axis.text.y=element_blank(),
        axis.text.x=element_text(color="black",size=10),
        legend.position="right",
        text=element_text(family="sans"),
        panel.background=element_rect(fill="white",color="white"),
        plot.title=element_text(hjust=0.5,size=14,face="bold"),
        plot.subtitle=element_text( hjust=0.5)
      )+ 
      labs(title=cluster_name )+
      coord_polar(theta='x',start=2*pi/10,clip='off') 
    
    p[[cluster_name]]=p1
  }
  
  leg <- get_legend(p[[1]])
  leg <- as_ggplot(leg)
  
  p[[11]] <- leg
  q <- ggpubr::ggarrange(plotlist=p,
                       ncol=4,nrow=3,font.label=list(color="black"),legend="none")
  q <- q+plot_annotation(
    title=title,
    theme=theme(
      plot.title=element_text(hjust=0.5)
    ))
  return(q)
}
