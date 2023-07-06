ACTION <- function(chr, locus_name, build, ld, highlight) {
	locus <- read.table(paste(locus_name,".region.sorted",sep=""), header=F, as.is=T)
	names(locus) <- c("CHR","SNP","POS","P")
	locus$P[locus$P <1E-300] <- 1E-300
	locus$P <- -log10(locus$P)
	maxP <- max(locus$P)
	locus$P[locus$P>40] <- 20+ (locus$P[locus$P>40]-20) * 20/(maxP-20) ## when logP > 40, re-scaling after 20
	if (ld =="N") {
		locus$COL <- "FFFFFF"
	} else {
                if (file.info(ld)$size ==0) {
                        locus$COL <- "#FFFFFF"
                } else {
                        LD <- read.table(ld, header=F, as.is=T)
                        names(LD) =c("SNP","LD")
                        LD$COL <- hsv(0,round(LD$LD,digits=2),1)
                        locus <- merge(locus,LD,all.x=T)
                        locus$COL[is.na(locus$COL)] <- "#FFFFFF" ## for SNPs with no LD, use no color
                }
        }
	locus$PCH <- 1
	locus$TEXT <- ""
	hit_index <- which(locus$P ==minP)
	locus$COL[hit_index] <- "blue" ## for the top SNP
	locus$PCH[hit_index]] <- 18
	locus$TEXT[hit_index]] <- paste(locus$SNP[hit_index],"(",signif(locus$P[hit_index],3),")")
	high_index <- which(locus$SNP ==highlight)
	locus$PCH[high_index] <- 18
	locus$TEXT[high_index] <- paste(locus$SNP[high_index], "(",signif(locus$P[high_index],3),")")
	annotation <- read.table(paste("/nfs/team151_data03/references_panel/1000g/chr",chr,".annotation_v2.20101123.nsSNPs",sep=""), header=F, as.is=T)
        names(annotation) <- c("CHR","POS","GENE","PROTEIN")
	locus <- merge(locus,annotation,all.x=T)
	ns_index <- which(!is.na(locus$GENE))
	locus$PCH[ns_index] <- 18 ## nsSNP
	ns_index2 <- which(!is.na(locus$GENE) & locus$P <5e-8)
	locus$TEXT[ns_index2] <- paste(locus$SNP[ns_index2], "(",signif(locus$P[ns_index2],3), ",", locus$GENE[ns_index2], ",", locus$PROTEIN[ns_index2], ")")
	min.pos <- min(locus$POS)
	max.pos <- max(locus$POS)
	size.pos <- max.pos - min.pos
	center.pos <- min.pos + ( size.pos / 2 )
	center.pos.kb <- round(center.pos / 1000)
	offset.pos.kb <- round(size.pos/3000)
	range = ceiling(maxP)
	locus <- locus[2:nrow(locus),]
	offset <- range * 0.2 # dedicates to the gene annotations
	big.range <- range + offset 
	ystart.gene <- - offset
	ystart.recomb <- - offset + (big.range / 8)

	# write.table(locus,file=paste("chr",chr,".",locus_name,".merged",sep=""),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE)
	pdf(paste(locus_name,".pdf",sep=""), w=12,h=8)
#	bitmap(paste("chr",chr,".",locus_name,".bmp",sep=""), w=10,h=5,res=256)
	par(mar=c(4,4,3,4))
	### Plot recombination rate and draw the graph frame
	recomb <- read.table(paste("/nfs/team151/jh21/data/hapmap/genetic_map/genetic_map_chr", chr, "_",build,".txt", sep=""), header=T)
	keep.recomb <- subset(recomb, recomb[,1] > min.pos & recomb[,1] < max.pos)
#	plot(keep.recomb[,1], ( keep.recomb[,2] / 60 ) * big.range,  type="l", col="lightblue", xlim=c(min.pos, max.pos), ylim=c(-offset,range), main=locus_name, xlab="", ylab="", axes=F)
	plot(keep.recomb[,1], ( keep.recomb[,2] / 100 ) * range,  type="l", col="lightblue", xlim=c(min.pos, max.pos), ylim=c(-offset,range), main=locus_name, xlab="", ylab="", axes=F)
	mtext(paste("Chromosome", chr, "position (kb)", sep=" "), side=1, line=2.5)
	axis(1, at=c(center.pos.kb*1000 - offset.pos.kb*1000, center.pos.kb*1000, center.pos.kb*1000 + offset.pos.kb*1000), labels=c((center.pos.kb - offset.pos.kb), center.pos.kb, center.pos.kb + offset.pos.kb), las=1) 
	if (maxP <=40 ) {
		at_pos <- c(seq(from=0,to=min(40,2*floor(range/2)),by=2))
                at_label <- at_pos
        } else {
                at_pos <- c(seq(from=0,to=20,by=2), range)
                at_label <- c(seq(from=0,to=20,by=2), maxP)
	}
	axis(2, at=at_pos, labels=at_label, cex.axis=1.5)
	axis(4, at=c( 0, range/4, range*2/4, range*3/4 ), labels=c("0","25","50","75"), las=1)
	box()
	lines(c(min.pos, max.pos), c(0,0), lty="dotted", lwd=1, col="black")
	mtext("-log10(P)", side=2, line=2, cex=1)
	mtext("Recombination rate (cM/Mb)", side=4, at=(-offset+(range+offset)/2), line=2, cex=1)
	## Plot the data
	points(locus$POS, locus$P, pch=locus$pch, col=locus$col, cex=1)
	if (as.numeric(hit$P) < 5E-08) {
		segments(min.pos,-log10(5E-08),max.pos,-log10(5E-08),lty=2,col="orange")
	}
	locus_text <- locus[which(locus$Text != ""),]
	text(locus_text$POS, -(log10(locus_text$P)), labels=locus_text$Text, pos=4, offset=1, col=locus_text$col)
	### Plot the genes
	genelist <- read.table(paste("/nfs/team151/jh21/files/refGene.",build,".txt",sep=""),as.is=T, header=F)[,1:13]
        names(genelist)=c("bin","name","CHR","STRAND","START","STOP","cdsSTART","cdsEnd","NEXONS","EXSTART","EXEND","score","GENE")
	chr <- paste("chr",chr,sep="") 
        genelist <- genelist[genelist$CHR==chr,]
	genelist <- genelist[order(genelist$START),]
	gene.before.locus <- genelist[genelist$STOP < min.pos,]
	gene.before.locus <- gene.before.locus[length(gene.before.locus[,1]),]
	dist.to.gene.before <- min.pos - gene.before.locus$STOP
	name.of.gene.before <- gene.before.locus$GENE
	gene.after.locus <- genelist[genelist$START > max.pos,]
	gene.after.locus <- gene.after.locus[1,]
	dist.to.gene.after <- gene.after.locus$START - max.pos
	name.of.gene.after <- gene.after.locus$GENE
	genes.in.locus <- subset(genelist, ( genelist$START > min.pos & genelist$START < max.pos ) | ( genelist$STOP > min.pos & genelist$STOP < max.pos) | ( genelist$START < min.pos & genelist$STOP > max.pos ) )
	genes.in.locus = genes.in.locus[order(genes.in.locus$START),]
	write.table(genes.in.locus,file=paste(locus_name,".gene",sep=""),append=FALSE,quote=FALSE,col.names=TRUE,row.names=FALSE)
	if ( nrow(genes.in.locus) >0) {
		for ( i in 1:nrow(genes.in.locus) ) { 
			starty <- -(offset/3) * ((i+2) %%3) -offset/10
		##	segments( max(genes.in.locus[i,]$START,min.pos), starty, min(genes.in.locus[i,]$STOP, max.pos), starty, lwd=2, lty="solid", col="darkgreen")
			if ( genes.in.locus[i,]$STRAND == "+" ) {
				arrows(max(genes.in.locus[i,]$START,min.pos), starty, min(genes.in.locus[i,]$STOP, max.pos), starty, length=0.08, lwd=2, code=2, lty="solid", col="darkgreen")
			} else {
				arrows(max(genes.in.locus[i,]$START,min.pos), starty, min(genes.in.locus[i,]$STOP, max.pos), starty, length=0.08, lwd=2, code=1, lty="solid", col="darkgreen")
			}	
			exon.n=genes.in.locus[i,]$NEXONS
			exon.start = as.numeric(strsplit(genes.in.locus[i,]$EXSTART,split=",")[[1]])
			exon.end = as.numeric(strsplit(genes.in.locus[i,]$EXEND,split=",")[[1]])
			if (!is.na(exon.n)) {
				for (j in 1:exon.n) {
					if (exon.start[j] >= min.pos & exon.end[j] <= max.pos) {
					##	segments(exon.start[j],starty,exon.end[j],starty,lwd=4,lty="solid",col="black")
					}
				}
			}
			if(!is.na(genes.in.locus[i,]$GENE)) {
				text(max(genes.in.locus[i,]$START,min.pos)+size.pos/50, starty - offset/8, labels=genes.in.locus[i,]$GENE, font=2, cex=1)
			}
		}
	} else {
		savefont <- par(font=2)
		legend("bottomleft", paste("<--",dist.to.gene.before,"from",name.of.gene.before,sep=" "), cex=1, bty="n")
		legend("bottomright", paste(dist.to.gene.after,"to",name.of.gene.after,"-->",sep=" "), cex=1, bty="n")
		par(savefont)
	}
	dev.off()
}
