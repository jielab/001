library(mvtnorm)

ACTION <- function(dir,pheno,ld,pval) {
	PHENO <- read.table(paste(dir,"files",pheno,sep="/"),row.names=1,header=T)
	PHENO.matrix <- as.matrix(PHENO,nrow=nrow(PHENO))
	
	LD.matrix <- as.matrix(read.table(paste(dir,"files",ld,sep="/"),header=F))
	diag(LD.matrix) <- 0
	mc <- apply(LD.matrix,2,max)
	
    #Removes perfect proxies
	while(max(mc) >= 1) {
		w <- (1:length(mc))[mc==1][1]
		LD.matrix <- LD.matrix[1:nrow(LD.matrix) != w,1:nrow(LD.matrix) !=w]
		mc <- apply(LD.matrix,2,max)
	}
	diag(LD.matrix) <- 1
	
    #Fix matrix to make it positive definite
	e0 <- eigen(LD.matrix)$values
	e <- e0
	V <- eigen(LD.matrix)$vectors
	Vinv <- solve(V)
	smallest.pos.e <- 0.0000000001
	while(min(e) < 0) {
		neg <- (e < smallest.pos.e)
		smallest.pos.e <- min(e[!neg])
		e[neg] <- smallest.pos.e
		D <- diag(e)
		LD.matrix <- V %*% D %*% Vinv
		e <- eigen(LD.matrix)$values
	}


	sets <- 1e+6
	sims <- max(1,0.001/pval) ## 10 sims for p=1e-4, 100 sims for p=1e-5, etc
	final <- data.frame(count=rep(0,ncol(PHENO.matrix)+2)) # first row for 0, last row for total # of simulations
	rownames(final) <- 0:(ncol(PHENO.matrix)+1)
	for (i in 1:sims) {
		Z <- rmvnorm(sets, rep(0,ncol(PHENO.matrix)*ncol(LD.matrix)), kronecker(PHENO.matrix,LD.matrix))
		P <- 2*pnorm(abs(Z),lower=F)
		rm(Z)
		P <- apply(P,2,function(x) x<=pval)
		P_byPheno <- matrix(0,nrow=sets,ncol=ncol(PHENO.matrix)) ## must put outside of the following loop
		for(j in 1:ncol(PHENO.matrix)) {
			P_byPheno[,j] <- apply( P[,((j-1)*ncol(LD.matrix)+1):(j*ncol(LD.matrix))], 1, function(x) max(x) ) 
		}
		P_cnt <- apply(P_byPheno,1,function(x) sum(x))
		rm(P_byPheno)
		final[names(table(P_cnt)),1] <- final[names(table(P_cnt)),1] + table(P_cnt)
		print(c(i,"finished"))
	}
	final[ncol(PHENO.matrix)+2,1] <- sims*sets
	write.table(final,file=paste(dir,"/",pheno,".",ld,".p",pval,".cnt",sep=""),quote=F)
}
# binom <- dbinom(0:max_traits,nrow(PHENO.matrix),1-(1-pval)^indep_SNPs)*sims
# binom <- -log10(binom/sims) 




dir <- "/data/home/jiehuang/prime/sim"
pheno_num <-9 
rounds <- 10
pvals <- c(1e-04,1e-05,1e-06,1e-07,5e-08)
phenos <- c("blood.cor.SH2B3")
for (pheno in 1:length(phenos)) {
	for (pval in 1:length(pvals)) {
		final <- data.frame(count=rep(0,pheno_num+2))
		rownames(final) <- 0:(pheno_num+1)
		for (r in 1:rounds) {	
			each_cnt <- read.table(paste(dir,"/",phenos[pheno],".ld.",pvals[pval],".",r,".cnt",sep=""),row.names=1,header=T)
			final[rownames(each_cnt),1] <- final[rownames(each_cnt),1] + each_cnt
		}
		write.table(final,file=paste(phenos[pheno],".ld.",pvals[pval],".cnt",sep=""),quote=F)
	}
}

