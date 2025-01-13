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

mr = ? # 1或0
cisMr = ? # 1或0
coloc = ? # 1或0
mrMed = ? # 1或0

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
glist <- read.table(paste0(dir0,'/data/gwas/ppp/ppp.b19.bed'), header=FALSE)

	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 开始 🏃‍
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (Y in Ys) { # 🙍
#	writeLines(paste('\n\n-->Run:', Y))
	dat.Y.raw <- read.table(paste0(dir.Y, '/', Y, '.gz'), header=T)
	names(dat.Y.raw) <- stri_replace_all_regex(toupper(names(dat.Y.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

#	cl <- makeCluster(detectCores()-1); registerDoParallel(cl); 
#	foreach(X = Xs) %dopar% {
	for (X in Xs) { # 🍷

	label=paste0(X, '-', Y)
	dat.X.raw <- read.table(paste0(dir.X, '/', X, '.gz'), header=T)
	names(dat.X.raw) <- stri_replace_all_regex(toupper(names(dat.X.raw)), pattern=toupper(pattern), replacement=replacement, vectorize_all=FALSE)

	if (mr==1) { # TwoSampleMR 🏮
		write("X Y n-IV BETA SE P.ivw p.hetero.egger p.hetero.weighted p.pleio P.outcome.min", file=paste0(label,'.mr.log'), append=FALSE)
		run_mr2s(label, X, dir.X, dat.X.raw, Y, dir.Y, dat.Y.raw)
		run_mr2s(label, Y, dir.Y, dat.Y.raw, X, dir.X, dat.X.raw) 
	}
	p.X2Y <- read.table(paste0(label,".mr.log"), header=TRUE)[1,"P.ivw"]

	if (cisMr==1) { # cisMr 🏮
		chrPos <- subset(glist, V4==X)
		if (nrow(chrPos)!=1) { write(paste(X,Y,nrow(dat.X.raw), "SKIP:", X, "does not have positions"), file=paste0(label,'.cisMr.log'), append=FALSE); next }
		chr=chrPos[1,1]; flank=100000; pos0=chrPos[1,2]; pos1=chrPos[1,3]
		write("X Y n-IV BETA SE P", file=paste0(label,'.cisMr.log'), append=FALSE)
		run_cisMr(label, chr, flank, pos0, pos1, X, dat.X.raw, dir.X.cojo, Y, dat.Y.raw)
	}
	
	if (coloc==1 & p.X2Y <=5e-05) { # coloc 🏮
		chrPos <- subset(glist, V4==X)
		if (nrow(chrPos)!=1) { write(paste(X,Y,nrow(dat.X.raw), "SKIP:", X, "does not have positions"), file=paste0(label,'.cisMr.log'), append=FALSE); next }
		chr=chrPos[1,1]; flank=100000; pos0=chrPos[1,2]; pos1=chrPos[1,3]
		write("X Y n-SNP H0 H1 H2 H3 H4", file=paste0(label,'.coloc.log'), append=FALSE)
		run_coloc(label, chr, flank, pos0, pos1, X, dat.X.raw, Y, dat.Y.raw)
	}
		
	if(mrMed==1 & p.X2Y <=5e-02) { # mrMed 🏮
		write("X M Y n-IV BETA.xy P.xy TE P.TE ACME P.ACME ADE P.ADE Prop P.Prop", file=paste0(label,'.mrMed.log'), append=FALSE)
		for (file.M in Ms.full) { run_mrMed(label, X, dat.X.raw, file.M, Y, dat.Y.raw) }
	}
	
	}
	# stopCluster(cl)
}

