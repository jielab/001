# Source https://github.com/Ben-JWL/DISH/tree/master/DISH
library(tidyverse)

dir1="D:/scripts/main"
output ="hla.imputed.txt"; lambda=0.15

dat = read.table(paste0(dir1,"/EUR_sample.txt"), header=T) %>% mutate(
	Z =ifelse(Tval >=0, abs(qnorm((2*pt(-abs(Tval), df=2, log.p=T))-log(2),lower.tail=FALSE, log.p=TRUE)), -abs(qnorm((2*pt(-abs(Tval), df=2, log.p=TRUE))-log(2), lower.tail=FALSE, log.p=TRUE)))
)
ref <- readRDS(paste0(dir1,"/European_hg19.Rds"))
	ref_map = ref[[2]] %>% rename(SNP=SNP_id, POS=POS19); ref_map$POS18 <- NULL
	filtered_idx = which(ref_map$maf >= 0.005)
	ref_map = ref_map[filtered_idx,]
	sigma = ref[[1]][filtered_idx, filtered_idx]; sigma=as.matrix(sigma)
	# nrow(ref_map); nrow(sigma); rownames(sigma)
	
dat1 <- merge(x=ref_map, y=dat, by='POS')
	dat1.pos <- dat1[(dat1$Major==dat1$EA & dat1$Minor==dat1$NEA) | (dat1$Major==dat1$NEA & dat1$Minor==dat1$EA),"POS"]
	dat1.pos <- sort(dat1.pos) 

typed_idx <- which(ref_map$POS %in% dat1.pos)
	sigma_tt = sigma[typed_idx, typed_idx]
	sigma_tt = sigma_tt + lambda*diag(ncol(sigma_tt))
	sigma_it = sigma[-typed_idx, typed_idx]
	det(sigma_tt)
	sigma_tt_inv = solve(sigma_tt)
	weight = sigma_it %*% sigma_tt_inv
	Zi = weight %*% Zt
	var = matrix(NA,nrow(weight),1)
	for(i in 1:nrow(weight)) {
		var[i] = sum(weight[i,] %o% weight[i,] * sigma_tt )
	}

typed_in_both_idx <- which(dat$POS %in% dat1.pos)
	dat = dat[typed_in_both_idx,]
	Zt = rep(NA, length(dat$Z))
	ref_map_t = ref_map[typed_idx,]
	ref_map_i = ref_map[-typed_idx,]
	for (i in 1:length(Zt)) {
		if(as.character(dat$EA[i]) != as.character(ref_map_t$Minor[which(ref_map_t$POS == dat$POS[i])]) ) {
			Zt[which(ref_map_t$POS == dat$POS[i])] = -dat$Z[i]
		} else {
			Zt[which(ref_map_t$POS == dat$POS[i])] = dat$Z[i]
		}
	}
	Zt = as.matrix(Zt)

results = cbind(as.character(ref_map_i$SNP), as.character(ref_map_i$POS), as.character(ref_map_i$Minor), as.character(ref_map_i$Majo), Zi, var, 2*pnorm(-abs(Zi)))
	colnames(results) = c("Marker_id", "Marker_pos", "EA", "NEA", "Imputed_Z",  "r2pred", "imputed_P")
	write.table(results, output, row.names=F, col.names=T, quote=F, sep="\t")