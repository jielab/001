# devtools::install_github("scllin/mrMed") 
pacman::p_load(data.table, tidyverse, stringi, TwoSampleMR, MendelianRandomization, mrMed, cisMRcML) # foreach, doParallel
dir0='/work/sph-huangj/'
source(paste0(dir0, 'scripts/main/f.R'))
dir=paste0(dir0, 'data/gwas/')

n_cores=40
folder.X = '?'
folder.Y = '?'
folder.M = '?' # 有可能有多个folders，用|分割
Xs = '?' # 一个X，或 ALL 表示全部
Ys = '?' # 一个Y，或 ALL 表示全部
label='?'; log_mr=paste0(label,'.mr.log')

cis=?
cisMRcML=?; log_cisMrcML=paste0(label,'.cisMrcML.log')
coloc=?; log_coloc=paste0(label,'.coloc.log')
mrMed=?; log_mrMed=paste0(label,'.mrMed.log')


dir.X = paste0(dir,folder.X,'/clean')
dir.Y = paste0(dir,folder.Y,'/clean')
if (grepl("ALL", Xs)) {Xs = gsub('.gz', '', list.files(path=dir.X, pattern='.gz$'))}
if (grepl("ALL", Ys)) {Ys = gsub('.gz', '', list.files(path=dir.Y, pattern='.gz$'))}
Ms.full = NULL # 如果 folder.M 不是 NA，后面会自动给 Ms 赋值
if (!grepl("NA", folder.M)) {
	for (folder.m in str_split_1(folder.M, ' ')) {
		Ms.full <- append(Ms.full, list.files(path=paste0(dir,folder.m,'/clean'), pattern='.gz$', full.names=TRUE))	
	}
}

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 进行分析
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

	for (X in Xs) { # 🍷
		dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
		names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
		if (cis==1) { # 下面所有的分析，只保留X 的 cis-SNP 数据。 所以，只有protein作为X的时候，才设置 cis=1
			dat.X.raw <- subset(dat.X.raw, CHR==chr & POS>=(pos_begin-flank) & POS<=(pos_end+flank)) 
		}
		
		##🗡 基于显著的SNP，进行常规MR分析 🏮🏮
		if (file.exists(paste0(dir.X, '/', X, '.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.top.snp'), header=T); names(dat.X.iv) <- "SNP" 
		} else if (file.exists(paste0(dir.X, '/', X, '.NEW.top.snp'))) {
			dat.X.iv <- read.table(paste0(dir.X, '/', X, '.NEW.top.snp'), header=T)
		} else {
			dat.X.sig <- dat.X.raw %>% filter(P<=5e-08) %>% mutate(mb=ceiling(POS/1e+05))
			dat.X.iv <- dat.X.sig %>% group_by(mb) %>% slice(which.min(P)) %>% ungroup() %>% select("SNP")
			write.table(dat.X.iv, paste0(dir.X, '/', X, '.NEW.top.snp'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
		}
		if(nrow(dat.X.iv) ==0) {write(paste("SKIP:", X,"NA",Y, "X no IV !!"), file=log_mr, append=TRUE); next}	
		dat.X <- dat.X.raw %>% merge(dat.X.iv) %>% format_data(type='exposure', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.exposure=X)
		dat.Y <- dat.Y.raw %>% merge(dat.X.iv, by="SNP") %>% format_data(type='outcome', snp_col='SNP', effect_allele_col='EA', other_allele_col='NEA', beta_col='BETA', se_col='SE', pval_col='P') %>% mutate(id.outcome =Y)		
		dat.XY <- harmonise_data(dat.X, dat.Y, action=1) 
		write.table(dat.XY, paste0(label, '.dat'), append=FALSE, quote=FALSE, row.names=FALSE, col.names=TRUE)
		fit.X2Y <- mr(dat.XY, method_list=c("mr_wald_ratio", "mr_ivw")); fit.X2Y 
			beta.X2Y <- fit.X2Y$b; se.X2Y <- fit.X2Y$se; p.X2Y <- fit.X2Y$pval
			p.hetero <- mr_heterogeneity(dat.XY)$Q_pval
            p.pleio <- mr_pleiotropy_test(dat.XY)$pval
			write("X Y M BETA SE P.ivw p.hetero.egger p.hetero.weighted p.pleio", file=log_mr, append=FALSE)
            write(paste(X,Y, nrow(dat.XY), rb(beta.X2Y), rb(se.X2Y), rp(p.X2Y), rp(p.hetero[1]), rp(p.hetero[2]), rp(p.pleio)), file=log_mr, append=TRUE)

		##🗡 基于某个locus的所有SNP，进行更加精准的cisMR-cML分析 🏮🏮
		if (cisMRcML==1) {
			dat$exp_df$cor = dat$exp_df$b / sqrt(dat$exp_df$b^2 + (dat$exp_df$N - 2) * dat$exp_df$se^2)
			dat$exp_df$bJ = solve(dat$LD_mat) %*% dat$exp_df$cor
			dat$out_df$cor = dat$out_df$b / sqrt(dat$out_df$b^2 + (dat$out_df$N - 2) * dat$out_df$se^2)
			dat$out_df$bJ = solve(dat$LD_mat) %*% dat$out_df$cor
			b_exp_cond=dat$exp_df$bJ
			b_out_cond=dat$out_df$bJ 
			Sig_exp1 = solve(dat$LD_mat) %*% (dat$exp_df$se_cor %o% dat$exp_df$se_cor * dat$LD_mat) %*% solve(dat$LD_mat)
			Sig_out1 = solve(dat$LD_mat) %*% (dat$out_df$se_cor %o% dat$out_df$se_cor * dat$LD_mat) %*% solve(dat$LD_mat)
			Sig_exp_inv=solve(Sig_exp1)
			Sig_out_inv=solve(Sig_out1)
			res = cismr_cML_DP(
				b_exp=b_exp_cond, b_out=b_out_cond, Sig_exp_inv=Sig_exp_inv, Sig_out_inv=Sig_out_inv, maxit=200,
				n=dat$N1,random_start=5, min_theta_range=-0.1, max_theta_range=0.1, num_pert=100, random_start_pert=5, random_seed=12345
			)			  
			write("X Y m BETA SE P", file=log_cisMrcML, append=FALSE)
            write(paste(X,Y, nrow(dat), rb(res$BIC_DP_theta), rb(es$BIC_DP_se), rp(res$BIC_DP_p)), file=log_cisMrcML, append=TRUE)
		}
		
		##🗡 对MR显著的locus，进行coloc分析 🏮🏮
		if (coloc==1 & p.X2Y <=0.05) {
			dat <- harmonise_data(dat.X, dat.Y, action=1) 
			dat.coloc.X <- dat %>% {list(type="quant", snp=.$SNP, position=.$pos.exposure, MAF=.$eaf.exposure, N=.$samplesize.exposure, beta=.$beta.exposure, varbeta =.$se.exposure^2, pvalues=.$pval.exposure)}
			dat.coloc.Y <- dat %>% {list(type="cc", snp=.$SNP, position=.$pos.outcome, MAF=.$eaf.outcome, N=100000, beta=.$beta.outcome, varbeta =.$se.outcome^2, pvalues=.$pval.outcome)}
			fit.coloc <- coloc::coloc.abf(dat.coloc.X, dat.coloc.Y)
			t(as.data.frame(fit.coloc[["summary"]])) 
			par(mfrow = c(2,1)); plot_dataset(dat.coloc.X); plot_dataset(dat.coloc.Y)
			sensitivity(fit.coloc,"H4>0.9")
		}
		
		##🗡 对MR显著的locus，进行中介分析 🏮🏮
		if(mrMed==1 & p.X2Y <=0.05) {
			# numCores <- detectCores()-1; cl <- makeCluster(numCores); registerDoParallel(cl)
			# foreach::foreach(M = Ms) %dopar% {
			for (file.M in Ms.full) {
				mediation(X, dat.X.raw, dat.X.iv, file.M, Y, dat.Y.raw, log_mrMed) 
			}
			# stopCluster(cl)
		}	
	}
}

