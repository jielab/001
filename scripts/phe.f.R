rb <- function(x) (if(is.null(x)) return(NA) else round(x,3))
rp <- function(x) (if(is.null(x)) return(NA) else signif(x,2))
inormal <- function(x) qnorm((rank(x, na.last='keep') - 0.5) / sum(!is.na(x)))
std <- function(x) (x - mean(x,na.rm=T)) / sd(x,na.rm=T)
remove_outlier <- function(x) (ifelse((x > (mean(x,na.rm=TRUE) + 3*sd(x,na.rm=TRUE)) | x < (mean(x,na.rm=TRUE) - 3*sd(x,na.rm=TRUE))), NA, x))
hardcall <- function(x) ifelse(x<0.5, 0, ifelse(x<1.5, 1, 2))
expo <- function(x) 1.1^x

scale10 <- function(x) {
  scale_f <- function(y) {scaled <- scale(y); scaled <- 2 + (scaled - min(scaled)) * (10-2) / (max(scaled) - min(scaled)); return(scaled[-1])}
  if (max(x) < 10) {return(x)} else {x_scaled <- ifelse(x<=2, x, scale_f(c(2,x))); return(x_scaled)}
}

bbplot <- function( # 🏄‍ 
  label, dat, X, BETA1_col, BETA2_col, BETA_compare, P1_col, P2_col, P_compare, sig.level) {
	dat$X <- dat[[X]]; dat$BETA1 <- scale(dat[[BETA1_col]]); dat$BETA2 <- scale(dat[[BETA2_col]])
	dat$P1 <- dat[[P1_col]]; dat$P2 <- dat[[P2_col]]
	dat1 <- dat %>% filter(P1 <=sig.level, !is.na(P1), !is.na(P2)) %>% 
		mutate(Plog10.1=scale10(-log10(P1)), Plog10.2=scale10(-log10(P2)))
	top5 <- dat1 %>% arrange(desc(Plog10.1)) %>% head(5)
	if (BETA_compare=="yes") { 
		gg <- ggplot(dat1, aes(x=BETA1, y=BETA2)) 
	} else if (P_compare=="yes") { 
		gg <- ggplot(dat1, aes(x=Plog10.1, y=Plog10.2))
	}
	gg + geom_point(size=2, shape=1) + ggtitle(label) + geom_text(data=top5, aes(label=X), size=3, color='darkblue', fontface='bold') + 
		theme_minimal() + theme(axis.text = element_text(size=14, face='bold'), axis.title = element_text(size=14, face='bold'), axis.line=element_line(linewidth=1.2), legend.position='none', plot.title=element_text(size=16, face='bold', hjust=0.5))
}

valcano <- function( # 🌋
  label, dat, X_col, BETA_col, P_col, sig.level, label_x, label_y) {
	dat$X <- dat[[X_col]]
	dat$BETA <- scale(dat[[BETA_col]])
	dat$P <- dat[[P_col]]
	dat <- dat %>% mutate(color=ifelse(P < sig.level & BETA > 0, 'positive', ifelse(P < sig.level & BETA < 0, 'negative', 'NS')))
	top5 <- dat %>% arrange(P) %>% head(10); top5$X=toupper(gsub('prot_', '', top5$X))
	ggplot(dat, aes(x=BETA, y=-log10(P), color=color)) + geom_point(size=2) +
		scale_color_manual(values=c('positive'='purple', 'negative'='green', 'NS'='gray')) +
		geom_text(data=top5, aes(label=X), hjust=0, vjust=-1.5, size=3, color='darkblue', fontface='bold') +
	#	geom_segment(data=top5, aes(x=BETA + 0.05, y=-log10(P) + 1, xend=BETA, yend=-log10(P)), color='black') +
		geom_hline(yintercept=-log10(sig.level), linetype='dashed', color='red', linewidth=1.2) +
		geom_vline(xintercept=0, linetype='dotted', linewidth=1.2) + 
		labs(x=label_x, y=label_y, title=label) +
		theme_minimal() + theme(axis.text=element_text(size=12, face='bold'), axis.title=element_text(size=14, face='bold'), axis.line=element_line(linewidth=1.2), legend.position='none', plot.title=element_text(size=16, face='bold', hjust=0.5))
}

splitMatch <- function(x, arr, sep = "|") {
  sapply(x, function(val) any(as.numeric(strsplit(as.character(val), split = sep, fixed = TRUE)[[1]]) %in% arr))
}

rowMeans2 <- function(x) {
	row_means <- rowMeans(x, na.rm=TRUE)
	row_means[is.nan(row_means)] <- NA #🏮如果所有的row都是NA，会给出NaN或者0⚡
	row_means[row_means==0 & rowMeans(is.na(x))==ncol(x)] <- NA
	return(row_means)
}

rowSums2 <- function(x) { 
  row_sums <- rowSums(x, na.rm=TRUE)
  row_sums[is.nan(row_sums)] <- NA 
  row_sums[row_sums==0 & rowSums(is.na(x))==ncol(x)] <- NA
  return(row_sums)
}

bulk_rowMeans <- function(dat, fields_names) {
	dat %>% mutate( purrr::imap_dfc(fields_names, function(x, name) {
	prefix <- strsplit(x, "\\|")[[1]][1]
	new_col <- rowMeans2(select(dat, starts_with(prefix)))
	setNames(as.data.frame(new_col), strsplit(x, "\\|")[[1]][2]) }))
}


