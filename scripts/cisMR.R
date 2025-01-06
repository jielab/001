pacman::p_load(data.table, tidyverse, stringi, TwoSampleMR, MendelianRandomization, cisMRcML, coloc, mrMed) # foreach, doParallel

allele.qc = function( # 🏮
  a1,a2,ref1,ref2) {
	strand_flip = function(ref) {flip = ref; flip[ref == "A"] = "T"; flip[ref == "T"] = "A"; flip[ref == "G"] = "C"; flip[ref == "C"] = "G"; flip}
	flip1 = strand_flip(ref1); flip2 = strand_flip(ref2)
	snp = list()
	snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C")) # remove strand ambiguous
	snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F # remove non-ATCG coding
	snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
	snp[["sign_flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
	snp[["strand_flip"]] = (a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)
	exact_match = (a1 == ref1 & a2 == ref2) # remove tri-allelic
	snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] = F
	return(snp)
}

label="ABO-CAD"; log_file=paste0(label, ".cisMr.log")
chr=9; flank=100000; pos0=136108665-flank; pos1=136151440+flank
dir0="/work/sph-huangj/data/gwas"
X="ABO"; dir.X=paste0(dir0,"/ppp/clean"); dir.X.cojo=paste0(dir0,"/ppp/cis") 
Y="y.cad"; dir.Y=paste0(dir0,"/main/clean")
dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)

header <- c('SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N')
header.y <- c(header[1], paste0(header[-1], ".y"))
bfile <- paste("--bfile", paste0(dir.X.cis, "/", X), "--chr", chr) 

dat.X <- dat.X.raw %>% filter(CHR==chr, POS>=pos0, POS<=pos1) %>% select(SNP,EA,NEA,EAF,BETA,SE,P,N); names(dat.X) <- header
dat.Y <- dat.Y.raw %>% filter(SNP %in% dat.X$SNP) %>% select(SNP,EA,NEA,EAF,BETA,SE,P); dat.Y$N <- 100000; names(dat.Y) <- header.y
dat.XY <- merge(dat.X, dat.Y, by='SNP')
	remove_flip = allele.qc(dat.XY$A1, dat.XY$A2, dat.XY$A1.y, dat.XY$A2.y)
	dat.XY$b[which(remove_flip$sign_flip)] = -dat.XY$b[which(remove_flip$sign_flip)]
	dat.XY$freq[which(remove_flip$sign_flip)] = 1-dat.XY$freq[which(remove_flip$sign_flip)]
	dat.XY <- dat.XY[remove_flip$keep,]
	dat.X <- dat.XY %>% select(all_of(header))
	dat.Y <- dat.XY %>% select(all_of(header.y)); names(dat.Y) <- header
	fwrite(dat.X, paste0(X,'.4gcta'), sep='\t')
	fwrite(dat.Y, paste0(Y,'.4gcta'), sep='\t')
				
X.cojo.fn <- paste0(dir.X.cojo, '/', X, '.jma.cojo')
    if(!file.exists(X.cojo.fn)) { print(paste('Skip', X, 'jma.cojo files does not exist')); next }
	X.cojo <- fread(X.cojo.fn); if(length(X.cojo$SNP)<3) { print(paste('Skip', X, 'pqtal_cojo_res less than 3 records')); next }

Y.cojo.fn <- paste0(Y, '.jma.cojo')
	Y.cojo.cmd <- paste("gcta", bfile, "--cojo-file", paste0(Y,'.4gcta'), "--cojo-slct --cojo-p 5e-8 --out", Y); system(Y.cojo.cmd)
    if(file.exists(Y.cojo.fn)) {Y.cojo = fread(Y.cojo.fn); IV = union(X.cojo$SNP, Y.cojo$SNP)} else {IV = X.cojo$SNP}

write.table(IV, paste0(X, '.iv'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=FALSE)

XY.cojo.cmd = paste("gcta", bfile, "--cojo-file", paste0(X,'.4gcta'), "--extract", paste0(X, ".iv --cojo-joint --out"), X) # 🏮
    p = system(XY.cojo.cmd, intern=F); if(p==1) { print(paste('Collinear-', X, Y)); next }

LD_mat = fread(paste0(X, '.ldr.cojo'))
LD_mat = LD_mat[,2:(ncol(LD_mat)-1)]; LD_mat = as.matrix(LD_mat); rownames(LD_mat) = colnames(LD_mat)

dat.X <- dat.X %>% filter(SNP %in% IV)
	dat.X <- dat.X[match(colnames(LD_mat), dat.X$SNP),]
	dat.X <- dat.X %>% mutate(p = as.numeric(p), z = b/se, cor.X = b / sqrt(b^2 + (N - 2) * se^2), se_cor = se * cor.X/b, bJ = solve(LD_mat) %*% cor.X)

dat.Y = dat.Y %>% filter(SNP %in% IV)
	dat.Y = dat.Y[match(colnames(LD_mat), dat.Y$SNP),]
    dat.Y <- dat.Y %>% mutate(p = as.numeric(p), z = b/se, cor.Y = b / sqrt(b^2 + (N - 2) * se^2), se_cor = se * cor.Y/b, bJ = solve(LD_mat) %*% cor.Y)

dat.mr = list(b_exp=dat.X$bJ, b_out=dat.Y$bJ, N1=median(dat.X$N), N2=median(dat.Y$N), LD_mat = LD_mat, exp_df = dat.X, out_df = dat.Y, exp_IV=X.cojo$SNP, out_IV=setdiff(IV, X.cojo$SNP))
	data("mr_dat"); str(mr_dat) # 与软件自带的数据比较📏

save(dat.mr, file=paste0(X, '.RData'))

b_exp_cond=dat.mr$exp_df$bJ
b_out_cond=dat.mr$out_df$bJ 
Sig_exp1 = solve(dat.mr$LD_mat) %*% (dat.mr$exp_df$se_cor %o% dat.mr$exp_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat)
Sig_out1 = solve(dat.mr$LD_mat) %*% (dat.mr$out_df$se_cor %o% dat.mr$out_df$se_cor * dat.mr$LD_mat) %*% solve(dat.mr$LD_mat)
Sig_exp_inv=solve(Sig_exp1)
Sig_out_inv=solve(Sig_out1)
res = cismr_cML_DP(
	b_exp=b_exp_cond, b_out=b_out_cond, Sig_exp_inv=Sig_exp_inv, Sig_out_inv=Sig_out_inv, maxit=200, n = dat.mr$N1, 
	random_start = 5, min_theta_range=-0.1,max_theta_range=0.1, num_pert=100,random_start_pert=5,random_seed = 12345
)

res$BIC_DP_theta; res$BIC_DP_se; res$BIC_DP_p # 最终结果 🎇