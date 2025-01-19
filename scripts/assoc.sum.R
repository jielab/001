pacman::p_load(data.table, tidyverse, stringi, TwoSampleMR, MendelianRandomization, cisMRcML, coloc, mrMed) # foreach, doParallel

dir0='/work/sph-huangj'
source(paste0(dir0, '/scripts/f/f.R'))
dir=paste0(dir0, '/data/gwas/')

n_cores=40; flank=500000
folder.X = '?'
folder.Y = '?'
folder.M = '?' # 有可能有多个folders，用|分割
Xs = '?' # 一个X，或 ALL 表示全部
Ys = '?' # 一个Y，或 ALL 表示全部

IV_filter = ?
mr = ?
cisMr = ?
coloc = ?
mrMed = ?


dir.X = paste0(dir,folder.X, '/clean')
dir.X.cojo = paste0(dir, 'ppp/cojo')
dir.Y = paste0(dir,folder.Y, '/clean')
if (grepl("ALL", Xs)) {Xs = gsub('.gz', '', list.files(path=dir.X, pattern='.gz$')); Xs=head(Xs)}
if (grepl("ALL", Ys)) {Ys = gsub('.gz', '', list.files(path=dir.Y, pattern='.gz$'))}
Ms.full = NULL # 如果 folder.M 不是 NA，后面会自动给 Ms 赋值
if (mrMed==1 & !grepl("NA", folder.M)) {
	for (folder.m in str_split_1(folder.M, ' ')) {
		Ms.full <- append(Ms.full, list.files(path=paste0(dir,folder.m,'/clean'), pattern='.gz$', full.names=TRUE))	
	}
}
glist <- read.table(paste0(dir0,'/data/gwas/ppp/ppp.b38.bed'), header=FALSE)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 开始 🏃‍
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	dat.Y.raw <- dat.Y.raw %>% mutate(CHR= ifelse(CHR==23, "X", CHR))
	if (!"N" %in% colnames(dat.Y.raw)) { dat.Y.raw$N <- 100000 }
	
#	cl <- makeCluster(detectCores()-1); registerDoParallel(cl); 
#	foreach(X = Xs) %dopar% {
	for (X in Xs) { # 🍷

	label=paste0(X, '-', Y)
	dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)
	dat.X.raw <- dat.X.raw %>% mutate(CHR= ifelse(CHR==23, "X", CHR))
	if (!"N" %in% colnames(dat.X.raw)) { dat.X.raw$N <- 100000 }

	if (mr==1) { # TwoSampleMR 🏮
		log_file=paste0(label,".mr.log")
		write("Analysis X Y n-IV BETA SE P.ivw p.hetero.egger p.hetero.weighted p.pleio P.outcome.min", file=paste0(label,'.mr.log'), append=FALSE)
		run_mr2s('genome-wide-IV.1way', label, X, dir.X, dat.X.raw, "NO", IV_filter, Y, dir.Y, dat.Y.raw)
		run_mr2s('genome-wide-IV.2way', label, Y, dir.Y, dat.Y.raw, "NO", IV_filter, X, dir.X, dat.X.raw) 
		clump08.fn=paste0(dir.X.cojo, "/", X, ".5e-08.clump")
		clump06.fn=paste0(dir.X.cojo, "/", X, ".5e-06.clump")
		if (file.exists(clump08.fn) && file.info(clump08.fn)$size > 0) {run_mr2s("cis.5e-08.clump", label, X, dir.X, dat.X.raw, clump08.fn, IV_filter, Y, dir.Y, dat.Y.raw)} else {write(paste(label, 'no clump.5e-08 SNPs'), file=log_file, append=TRUE)}
		if (file.exists(clump06.fn) && file.info(clump06.fn)$size > 0) {run_mr2s("cis.5e-06.clump", label, X, dir.X, dat.X.raw, clump06.fn, IV_filter, Y, dir.Y, dat.Y.raw)} else {write(paste(label, 'no clump.5e-06 SNPs'), file=log_file, append=TRUE)}
	}

	if (cisMr==1) { # cisMr 🏮
		log_file=paste0(label,".cisMr.log")
		chrPos <- subset(glist, V4==X)
		chr=chrPos[1,1]; flank=100000; pos0=chrPos[1,2]; pos1=chrPos[1,3]
		write("X Y n-IV BETA SE P", file=log_file, append=FALSE)
		if (nrow(chrPos)!=1) { write(paste(X,Y,nrow(dat.X.raw), "SKIP:", X, "does not have positions"), file=log_file, append=FALSE); next }
		run_cisMr(label, chr, flank, pos0, pos1, X, dat.X.raw, dir.X.cojo, Y, dat.Y.raw)
	}
	
	if(mrMed==1) { # mrMed 🏮
		log_file=paste0(label,'.mrMed.log')
		write("X M Y n-IV BETA.xy P.xy TE P.TE ACME P.ACME ADE P.ADE Prop P.Prop", file=log_file, append=FALSE)
		for (file.M in Ms.full) { run_mrMed(label, X, dat.X.raw, file.M, Y, dat.Y.raw) }
	}
	
	}
	# stopCluster(cl)
}

