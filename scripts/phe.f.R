pick_first_number <- function(x) {
	as.numeric(strsplit(x, "\\|")[[1]][1]) # for strings such as "1|2|3"
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
replace_if_equal <- function(dat, var_list, num) {
  sapply(var_list, function(pair) { fields <- strsplit(pair, "\\|")[[1]]; ifelse(dat[[fields[1]]]==num, dat[[fields[2]]], 0) })
}
