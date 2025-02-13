library(cisMRcML); library(tidyverse); library(data.table)

p_t="08" # p_t="08" 

chr=1; pos0=55039447; pos1=55064852; flank=100000
	X.cojo <- read.table(paste0("5e-", p_t, ".jma.cojo"), header=T)
	Y.cojo <- read.table(paste0("Y.5e-", p_t, ".jma.cojo"), header=T)
	IV = union(X.cojo$SNP, Y.cojo$SNP)
		
LD_mat.fn = paste0("5e-", p_t, '.ldr.cojo')
	LD_mat = fread(LD_mat.fn)
	LD_mat = LD_mat[,2:(ncol(LD_mat)-1)]; LD_mat = as.matrix(LD_mat); rownames(LD_mat) = colnames(LD_mat)
	
dat.X <- read.table("X.4gcta", header=T) %>% filter(SNP %in% colnames(LD_mat))
	dat.X$cor = dat.X$BETA / sqrt(dat.X$BETA^2 + (dat.X$N-2) * dat.X$SE^2)
	dat.X$se_cor = dat.X$SE * dat.X$cor/dat.X$BETA
	dat.X$bJ = solve(LD_mat) %*% dat.X$cor
	
dat.Y = read.table("Y.4gcta", header=T) %>% filter(SNP %in% colnames(LD_mat))
    dat.Y$cor = dat.Y$BETA / sqrt(dat.Y$BETA^2 + (dat.Y$N-2) * dat.Y$SE^2)
    dat.Y$se_cor = dat.Y$SE * dat.Y$cor/dat.Y$BETA
    dat.Y$bJ = solve(LD_mat) %*% dat.Y$cor
	
dat.mr = list(b_exp=dat.X$bJ, b_out=dat.Y$bJ, N1=median(dat.X$N), N2=median(dat.Y$N), LD_mat=LD_mat, exp_df=dat.X, out_df=dat.Y, exp_IV=X.cojo$SNP, out_IV=setdiff(IV, X.cojo$SNP))
	Sig_exp1 = solve(dat.mr$LD_mat) %*% (dat.mr$exp_df$se_cor %o% dat.mr$exp_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat); Sig_exp_inv=solve(Sig_exp1)
	Sig_out1 = solve(dat.mr$LD_mat) %*% (dat.mr$out_df$se_cor %o% dat.mr$out_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat); Sig_out_inv=solve(Sig_out1)
	res = cismr_cML_DP(
		b_exp=dat.mr$exp_df$bJ, b_out=dat.mr$out_df$bJ, Sig_exp_inv=Sig_exp_inv, Sig_out_inv=Sig_out_inv, maxit=200, n=dat.mr$N1, 
		random_start=5, min_theta_range=-0.1, max_theta_range=0.1, num_pert=100, random_start_pert=5, random_seed=12345
	)

print(paste(p_t, length(IV), res$BIC_DP_theta, res$BIC_DP_se, res$BIC_DP_p))