mrMed <- function(dat_mrMed, method_list=c("Diff_IVW","Prod_IVW","Prod_Median")){

	#dat_mrMed <- form_dat(dat_mrMed)
	
	res <- lapply(method_list, function(meth){get(meth)(dat_mrMed)})

	res_tab <- list(
	TE = cbind(method_list,plyr::rbind.fill(lapply(res, function(x) x$TE))),
	DE = cbind(method_list,plyr::rbind.fill(lapply(res, function(x) x$DE))),
	IE = cbind(method_list,plyr::rbind.fill(lapply(res, function(x) x$IE))),
	rho = cbind(method_list,plyr::rbind.fill(lapply(res, function(x) x$rho)))
	)

	res_tab <- lapply(res_tab, function(x){colnames(x)=c("methods","b","se","pval","CI_lower","CI_upper");return(x)})

	return(res_tab)
}


Mtd1 <- function(b_TE,se_TE,b_DE,se_DE,cov_TEDE=0,dat_mrMed,gamma=0.05){ 

	#critical value
	z=qnorm(1-gamma/2)

	#IE
	b_IE <- b_TE-b_DE
	se_IE <- sqrt(se_TE^2 + se_DE^2 - 2*cov_TEDE)

	#rho
	b_rho <- b_IE/b_TE
	se_rho <- sqrt(max(se_DE^2/b_TE^2 + b_DE^2*se_TE^2/b_TE^4 - 2*b_DE*cov_TEDE/b_TE^3, 0))

	return(list(
	TE=data.frame(b=b_TE,se=se_TE,pval=2*pnorm(-abs(b_TE/se_TE)),CI_lower=b_TE-z*se_TE,CI_upper=b_TE+z*se_TE),
	DE=data.frame(b=b_DE,se=se_DE,pval=2*pnorm(-abs(b_DE/se_DE)),CI_lower=b_DE-z*se_DE,CI_upper=b_DE+z*se_DE),
	IE=data.frame(b=b_IE,se=se_IE,pval=2*pnorm(-abs(b_IE/se_IE)),CI_lower=b_IE-z*se_IE,CI_upper=b_IE+z*se_IE),
	rho=data.frame(b=b_rho,se=se_rho,pval=2*pnorm(-abs(b_rho/se_rho)),CI_lower=b_rho-z*se_rho,CI_upper=b_rho+z*se_rho)
	))
}


Diff_IVW_0 <- function(dat_mrMed,gamma=0.05){


	dat_mrMed <- form_dat(dat_mrMed)

	#TE using Gx as IVs based on mr_ivw in TwoSampleMR::mr 
	indx_Gxy <- !(dat_mrMed$Gx==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	#DE 
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_DE <- res$coefficients[1,1]
	se_DE <- res$coefficients[1,2]/min(1,res$sigma)

	cov_TEDE <- 0

	return(Mtd1(res_3$b,res_3$se,b_DE,se_DE,cov_TEDE,dat_mrMed,gamma))
}

Diff_IVW <- function(dat_mrMed,gamma=0.05){


	dat_mrMed <- form_dat(dat_mrMed)

	#TE using Gx_plum as IVs based on mr_ivw in TwoSampleMR::mr 
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	#DE 
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[1,1]
	se_X <- res$coefficients[1,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	#Design matrix for cov_TEDE
	Gamma_mr <- replace(dat_mrMed$beta.X,which(indx_Gxy==0),0)
	W <- diag(replace(1/dat_mrMed$se.Y^2,which(indx_mvmr==0),0))
	Gamma_mvmr <- cbind(replace(dat_mrMed$beta.X,which(indx_mvmr==0),0),replace(dat_mrMed$beta.M,which(indx_mvmr==0),0))

	#UVMR and MVMR mutiplicative rnd effect error std
	sigma_uvmr <- summary(lm(dat_mrMed$beta.Y[indx_Gxy]~0+dat_mrMed$beta.X[indx_Gxy],weights=1/dat_mrMed$se.Y[indx_Gxy]^2))$sigma
	sigma_mvmr <- res$sigma

	cov_TEDE <- solve(t(Gamma_mvmr)%*%W%*%Gamma_mvmr)%*%t(Gamma_mvmr)%*%W%*%Gamma_mr/sum(t(Gamma_mr)%*%W%*%Gamma_mr)
	cov_TEDE <- cov_TEDE*max(1,min(sigma_uvmr^2,sigma_mvmr^2))

	return(Mtd1(res_3$b,res_3$se,b_X,se_X,cov_TEDE[1],dat_mrMed,gamma))
}

Diff_Egger <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	#orient direction based on beta of exposure
	dat_mrMed$beta.Y <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.Y
	dat_mrMed$beta.M <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.M
	dat_mrMed$beta.X <- sign(dat_mrMed$beta.X)*dat_mrMed$beta.X

	#TE using Gx_plum as IVs based on mr_egger_regression in TwoSampleMR::mr 
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <-TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_egger_regression"))

	#DE 
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	b_X <- res$coefficients[2,1]
	se_X <- res$coefficients[2,2]/min(1,res$sigma)
	pval_X <- 2*(1-pnorm(abs(b_X/se_X)))

	#Design matrix for cov_TEDE
	Gamma_mr <- cbind(replace(rep(1,dim(dat_mrMed)[1]),which(indx_Gxy==0),0), replace(dat_mrMed$beta.X,which(indx_Gxy==0),0))
	W <- diag(replace(1/dat_mrMed$se.Y^2,which(indx_mvmr==0),0))
	Gamma_mvmr <- cbind(replace(rep(1,dim(dat_mrMed)[1]),which(indx_mvmr==0),0),replace(dat_mrMed$beta.X,which(indx_mvmr==0),0),replace(dat_mrMed$beta.M,which(indx_mvmr==0),0))

	#UVMR and MVMR mutiplicative rnd effect error std
	sigma_uvmr <- summary(lm(dat_mrMed$beta.Y[indx_Gxy]~dat_mrMed$beta.X[indx_Gxy],weights=1/dat_mrMed$se.Y[indx_Gxy]^2))$sigma
	sigma_mvmr <- res$sigma	

	cov_TEDE <- solve(t(Gamma_mvmr)%*%W%*%Gamma_mvmr)%*%t(Gamma_mvmr)%*%W%*%Gamma_mr%*%solve(t(Gamma_mr)%*%W%*%Gamma_mr)
	cov_TEDE <- cov_TEDE*max(1,min(sigma_uvmr^2,sigma_mvmr^2))

	return(Mtd1(res_3$b,res_3$se,b_X,se_X,cov_TEDE[2,2],dat_mrMed,gamma))
}


Diff_Median <- function(dat_mrMed,gamma=0.05,Nboot=1000){

	dat_mrMed <- form_dat(dat_mrMed)

	#critical value
	z=qnorm(1-gamma/2)

	#TE using Gx_plum as IVs based on mr_weighted_median in TwoSampleMR::mr 
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_weighted_median"))
	b_TE <- res_3$b

	#point estimates of DE, IE and rho
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	b_DE <- quantreg::rq(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2)$coefficients[[1]]
	b_IE <- b_TE - b_DE
	b_rho <- b_IE/b_TE

	#===CI and se via Bootstrap
	DE_se=c(NA);IE_se=c(NA);rho_se=c(NA)
	IE_CI_boot=c(NA,NA);DE_CI_boot=c(NA,NA);rho_CI_boot=c(NA,NA)

	res_boot <- data.frame(TE_i=rep(0,Nboot),DE_i=rep(0,Nboot),IE_i=rep(0,Nboot),rho_i=rep(0,Nboot))
	dat_mrMed_i <- dat_mrMed
		for(i in 1:Nboot){
			indx_X <- !is.na(dat_mrMed_i$beta.X)
			indx_M <- !is.na(dat_mrMed_i$beta.M)
			indx_Y <- !is.na(dat_mrMed_i$beta.Y)
			dat_mrMed_i$beta.X[indx_X] <- as.vector(mvtnorm::rmvnorm(1,dat_mrMed$beta.X[indx_X],diag(dat_mrMed$se.X[indx_X]^2)))
			dat_mrMed_i$beta.M[indx_M] <- as.vector(mvtnorm::rmvnorm(1,dat_mrMed$beta.M[indx_M],diag(dat_mrMed$se.M[indx_M]^2)))
			dat_mrMed_i$beta.Y[indx_Y] <- as.vector(mvtnorm::rmvnorm(1,dat_mrMed$beta.Y[indx_Y],diag(dat_mrMed$se.Y[indx_Y]^2)))


			b_iv <-dat_mrMed_i$beta.Y[indx_Gxy]/dat_mrMed_i$beta.X[indx_Gxy]
			VBj <- (dat_mrMed_i$se.Y[indx_Gxy]^2)/(dat_mrMed_i$beta.X[indx_Gxy]^2) + (dat_mrMed_i$beta.Y[indx_Gxy]^2)*(dat_mrMed_i$se.X[indx_Gxy]^2)/(dat_mrMed_i$beta.X[indx_Gxy]^4)
			b_TE_i  <- TwoSampleMR::weighted_median(b_iv, 1 / VBj)
			#tau_mvmr <- mr(cbind(dat_XY_i,mr_keep),method_list=c("mr_weighted_median"))$b
			#tau_mvmr <- rq(gwasY_i[indx_Gx]~0+gwasX_i[indx_Gx],weights=1/gwasY_rmna$se.outcome[indx_Gx]^2)$coefficients[[1]]

			b_DE_i <- quantreg::rq(dat_mrMed_i$beta.Y[indx_mvmr]~0+dat_mrMed_i$beta.X[indx_mvmr]+dat_mrMed_i$beta.M[indx_mvmr],weights=1/dat_mrMed_i$se.Y[indx_mvmr]^2)$coefficients[[1]]

			res_boot[i,1] <- b_TE_i
			res_boot[i,2] <- b_DE_i
			res_boot[i,3] <- b_TE_i-b_DE_i
			res_boot[i,4] <- 1-b_DE_i/b_TE_i
		}
		
	se_TE <- sd(res_boot[,c("TE_i")])
	se_DE <- sd(res_boot[,c("DE_i")])
	se_IE <- sd(res_boot[,c("IE_i")])
	se_rho <- sd(res_boot[,c("rho_i")])
	TE_CI_boot=quantile(res_boot[,c("TE_i")],c(gamma/2,1-gamma/2))
	DE_CI_boot=quantile(res_boot[,c("DE_i")],c(gamma/2,1-gamma/2))
	IE_CI_boot=quantile(res_boot[,c("IE_i")],c(gamma/2,1-gamma/2))
	rho_CI_boot=quantile(res_boot[,c("rho_i")],c(gamma/2,1-gamma/2))
	names(IE_CI_boot) <- c("lower","upper")
	names(DE_CI_boot) <- c("lower","upper")
	names(rho_CI_boot) <- c("lower","upper")


	return(list(
	TE=data.frame(b=b_TE,se=se_TE,pval=2*pnorm(-abs(b_TE/se_TE)),CI_lower=TE_CI_boot[1],CI_upper=TE_CI_boot[2]),
	DE=data.frame(b=b_DE,se=se_DE,pval=2*pnorm(-abs(b_DE/se_DE)),CI_lower=DE_CI_boot[1],CI_upper=DE_CI_boot[2]),
	IE=data.frame(b=b_IE,se=se_IE,pval=2*pnorm(-abs(b_IE/se_IE)),CI_lower=IE_CI_boot[1],CI_upper=IE_CI_boot[2]),
	rho=data.frame(b=b_rho,se=se_rho,pval=2*pnorm(-abs(b_rho/se_rho)),CI_lower=rho_CI_boot[1],CI_upper=rho_CI_boot[2])
	))
}


Mtd2 <- function(b_alpha,se_alpha,b_beta,se_beta,b_tau,se_tau,gamma=0.05){ 

	#critical value
	z=qnorm(1-gamma/2)

	#TE
	b_TE <- b_tau
	se_TE <- se_tau

	#IE
	b_IE <- b_alpha*b_beta
	se_IE <- sqrt(b_alpha^2*se_beta^2 + b_beta^2*se_alpha^2)
	
	#DE
	b_DE <- b_tau-b_alpha*b_beta
	se_DE <- sqrt(se_tau^2 + se_IE^2)
	
	#rho
	b_rho <- b_IE/b_tau
	se_rho <- sqrt(se_IE^2/b_tau^2 + b_IE^2*se_tau^2/b_tau^4)	

	return(list(
	TE=data.frame(b=b_TE,se=se_TE,pval=2*pnorm(-abs(b_TE/se_TE)),CI_lower=b_TE-z*se_TE,CI_upper=b_TE+z*se_TE),
	DE=data.frame(b=b_DE,se=se_DE,pval=2*pnorm(-abs(b_DE/se_DE)),CI_lower=b_DE-z*se_DE,CI_upper=b_DE+z*se_DE),
	IE=data.frame(b=b_IE,se=se_IE,pval=2*pnorm(-abs(b_IE/se_IE)),CI_lower=b_IE-z*se_IE,CI_upper=b_IE+z*se_IE),
	rho=data.frame(b=b_rho,se=se_rho,pval=2*pnorm(-abs(b_rho/se_rho)),CI_lower=b_rho-z*se_rho,CI_upper=b_rho+z*se_rho)
	))
	
}



Prod_IVW_0 <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	#TE using Gx as IVs based on mr_ivw in TwoSampleMR::mr 
	indx_Gxy <- !(dat_mrMed$Gx==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))
	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))
	
	
	#exposure-mediator(alpha) using Gx as IVs based on mr_ivw in TwoSampleMR::mr 
	indx_Gxm <- !(dat_mrMed$Gx==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_ivw"))

	#mediator-outcome effect (beta) using MVMR
	indx_mvmr <- !(is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y)|dat_mrMed$G_mvmr==0)
	res <- summary(lm(dat_mrMed$beta.Y[indx_mvmr]~0+dat_mrMed$beta.X[indx_mvmr]+dat_mrMed$beta.M[indx_mvmr],weights=1/dat_mrMed$se.Y[indx_mvmr]^2))
	#b_DE <- res$coefficients[1,1]
	#se_DE <- res$coefficients[1,2]/min(1,res$sigma)

	return(Mtd2(res_1$b,res_1$se,res$coefficients[2,1],res$coefficients[2,2],res_3$b,res_3$se))
}



Prod_IVW <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)
 
	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha based on mr_ivw in TwoSampleMR::mr 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_ivw"))

	#estimate beta based on mr_ivw in TwoSampleMR::mr 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_ivw"))

	#estimate (TE) based on mr_ivw in TwoSampleMR::mr 
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_ivw"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se,gamma=0.05))
} 

Prod_Egger <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha based on mr_egger_regression in TwoSampleMR::mr 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_egger_regression"))

	#estimate beta based on mr_egger_regression in TwoSampleMR::mr 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_egger_regression"))

	#estimate (TE) based on mr_egger_regression in TwoSampleMR::mr 
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_egger_regression"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se,gamma=0.05))
} 


Prod_Median <- function(dat_mrMed,gamma=0.05){

	dat_mrMed <- form_dat(dat_mrMed)

	indx_Gxm <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.M))
	indx_Gmy <- !(dat_mrMed$Gm_plum==0|is.na(dat_mrMed$beta.M)|is.na(dat_mrMed$beta.Y))
	indx_Gxy <- !(dat_mrMed$Gx_plum==0|is.na(dat_mrMed$beta.X)|is.na(dat_mrMed$beta.Y))


	dat_XM <- dat_mrMed[indx_Gxm,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M")]

	dat_MY <- dat_mrMed[indx_Gmy,c("SNP","M","id.M","effect_allele.M","other_allele.M","eaf.M","beta.M","se.M","pval.M",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	dat_XY <- dat_mrMed[indx_Gxy,c("SNP","X","id.X","effect_allele.X","other_allele.X","eaf.X","beta.X","se.X","pval.X",
	"Y","id.Y","effect_allele.Y","other_allele.Y","eaf.Y","beta.Y","se.Y","pval.Y")]

	#estimate alpha based on mr_weighted_median in TwoSampleMR::mr 
	names(dat_XM) <- gsub("X","exposure",names(dat_XM))
	names(dat_XM) <- gsub("M","outcome",names(dat_XM))
	mr_keep <- rep(TRUE,dim(dat_XM)[1])
	res_1 <- TwoSampleMR::mr(cbind(dat_XM,mr_keep),method_list=c("mr_weighted_median"))

	#estimate beta based on mr_weighted_median in TwoSampleMR::mr 
	names(dat_MY) <- gsub("M","exposure",names(dat_MY))
	names(dat_MY) <- gsub("Y","outcome",names(dat_MY))
	mr_keep <- rep(TRUE,dim(dat_MY)[1])
	res_2 <- TwoSampleMR::mr(cbind(dat_MY,mr_keep),method_list=c("mr_weighted_median"))

	#estimate (TE) based on mr_weighted_median in TwoSampleMR::mr 
	names(dat_XY) <- gsub("X","exposure",names(dat_XY))
	names(dat_XY) <- gsub("Y","outcome",names(dat_XY))
	mr_keep <- rep(TRUE,dim(dat_XY)[1])
	res_3 <- TwoSampleMR::mr(cbind(dat_XY,mr_keep),method_list=c("mr_weighted_median"))

	return(Mtd2(res_1$b,res_1$se,res_2$b,res_2$se,res_3$b,res_3$se))
}