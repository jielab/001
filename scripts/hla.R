# 来源 https://github.com/Ben-JWL/DISH/tree/master/DISH
library(tidyverse)
dir1="D:/scripts/main"
output ="hla.imputed.txt"
lambda = 0.15

gwas <- read.table(paste0(dir1,"/W.EUR_sample.txt"), header=T) # 需要有Z值 🏮
ref <- readRDS(paste0(dir1,"/W.European_hg19.Rds"))
	sigma = as.matrix(ref[[1]]) # 一个 N*N 的 matrix
	ref_map = ref[[2]] # matrix里面SNP的基本信息
	merged <- merge(x=ref_map, y=gwas, by='SNP') 

typed_idx <- which(ref_map$SNP %in% merged.SNP) 
	sigma_typ = sigma[typed_idx, typed_idx]
	sigma_typ = sigma_typ + lambda * diag(ncol(sigma_typ)); det(sigma_typ)
	sigma_imp = sigma[-typed_idx, typed_idx]
	imputed = ref_map[-typed_idx,]

sigma_typ_inv = solve(sigma_typ) # make an inverse matrix of sigma_typ
	weight = sigma_imp %*% sigma_typ_inv
	Zi = weight %*% Z
	var = matrix(NA, nrow(weight), 1)
	for(i in 1:nrow(weight)) { var[i] = sum(weight[i,] %o% weight[i,] * sigma_typ ) }

results = cbind(imputed$SNP, imputed$POS, imputed$Minor, imputed$Majo, Zi, var, 2*pnorm(-abs(Zi)))
	colnames(results) = c("SNP", "POS", "EA", "NEA", "Z.imputed",  "r2pred", "P.imputed")
	write.table(results, output, row.names=F, col.names=T, quote=F, sep="\t")