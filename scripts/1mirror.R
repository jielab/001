pacman::p_load(data.table, tidyverse, patchwork, geosphere, splines, sf, leaflet, vegan)  # , RColorBrewer

dir0 = "D:"
source(paste0(dir0, '/scripts/f/phe.f.R'))

dat0 <- readRDS(file = paste0(dir0, '/data/ukb/phe/Rdata/all.plus.rds'))
dat <- dat0 %>% filter(prot.yes == 1, !is.na(home_east), !is.na(PC1), !is.na(deprive)) %>% 
	st_as_sf(., coords = c("home_east", "home_north"), crs = 27700) %>% st_transform(., crs = 4326) %>% 
	mutate(lon = st_coordinates(geometry)[,1], lat = st_coordinates(geometry)[,2]) %>% st_drop_geometry(.)

leaflet(data = dat) %>% addTiles() %>% addCircleMarkers( lng = ~lon, lat = ~lat,
    color = ~colorNumeric(palette = c("green", "red"), domain  = dat$center)(center), radius = 3,
	popup = ~paste("ID:", eid, "<br>", "Center:", center)
)

dat1 <- dat %>% mutate(deprive_scaled = rescale(deprive, to = c(0.3, 1)))
	centers <- sort(unique(dat1$center)); n_centers <- length(centers)
	hues <- seq(15, 375, length.out = n_centers + 1)[1:n_centers]; names(hues) <- centers
	dat1 <- dat1 %>% mutate(hue = hues[center], color = hcl(h = hue, c = 80, l = 100 - 60 * deprive_scaled))

pal <- colorFactor(palette = "Set1", domain = dat1$center)
leaflet(data = dat1) %>% addTiles() %>%
	addCircleMarkers(lng = ~lon, lat = ~lat, color = ~color, radius = 3, popup = ~paste("ID:", eid, "<br>", "Center:", center)) %>%
	addLegend("bottomright", pal = pal, values = ~center, title = "Center")

analyze_center <- function(center_id, dat, outcome_var) {
	df <- dat %>% filter(center == center_id); if (nrow(df) < 100) return(NULL)
	city_point <- c(mean(df$lon, na.rm = TRUE), mean(df$lat, na.rm = TRUE))
	df <- df %>% mutate(dist2 = distHaversine(matrix(c(lon, lat), ncol = 2), city_point))
	y <- df[[outcome_var]]; v <- var(y, na.rm = TRUE); if (is.na(v) || v == 0) return(NULL)
	lm_model <- lm(as.formula(paste(outcome_var, "~ dist2")), data = df)
	lm_summary <- tryCatch(tidy(lm_model, conf.int = TRUE) %>% filter(term == "dist2"), error = function(e) return(NULL))
	if (is.null(lm_summary)) return(NULL)
	spline_model <- lm(as.formula(paste(outcome_var, "~ ns(dist2, df = 3)")), data = df)
	spline_rsq <- summary(spline_model)$r.squared
	spline_p <- tryCatch(anova(spline_model)[["Pr(>F)"]][1], error = function(e) NA_real_)
	tibble(
		outcome = outcome_var, center = center_id, n = nrow(df),
		lm_slope = lm_summary$estimate, lm_p = lm_summary$p.value, lm_rsq = summary(lm_model)$r.squared,
		conf_low = lm_summary$conf.low, conf_high = lm_summary$conf.high, spline_rsq = spline_rsq, spline_p = spline_p
	)
}


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Genes (PC1-40) 🪞 geography
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1 <- scale(dat[, c("lon", "lat", paste0("PC", 1:40))])
res <- list(); cors <- c()
for (i in 1:39) {
	proc <- procrustes(dat1[, c("lon", "lat")], dat1[, paste0("PC", i:(i+1))])
	r <- cor(proc$X[,1], proc$Yrot[,1]) + cor(proc$X[,2], proc$Yrot[,2])
	res[[i]] <- list(pair = paste0("PC", i, "-PC", i+1), proc = proc, r = r, pc_index = i)
	cors[i] <- r
}
results <- data.frame(pair = sapply(res, function(x) x$pair), r = sapply(res, function(x) x$r))
head(results[order(-results$r), ])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Proteins 🪞 geography
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
prots <- grep("^prot_", names(dat), value = TRUE)
prots.good <- prots %>% purrr::keep(~ { x <- dat[[.x]]; sum(!is.na(x) & is.finite(x)) > 10000 && var(x, na.rm = TRUE) > 0 })
dat1 <- dat %>% dplyr::select(all_of(prots.good)) %>%
	mutate(across(all_of(prots.good), ~ ifelse(is.na(.) | !is.finite(.), mean(., na.rm = TRUE), .))) %>% 
	scale()# 🏮

library(irlba); svd_res <- irlba(dat1, nv = 40) # 需要 matrix
	pc40 <- as.data.frame(svd_res$u %*% diag(svd_res$d))
	colnames(pc40) <- paste0("PC", 1:40); 
	pc40 <- pc40 %>% bind_cols(eid = dat$eid, lon = dat$lon, lat = dat$lat, .)
	loadings <- svd_res$v[, 1]; names(loadings) <- colnames(dat1)
	sort(loadings, decreasing = TRUE)[1:20]; sort(loadings, decreasing = FALSE)[1:20]

results <- list(); interim_file = "prot.geo.log"
Ys <- grep("deprive|^prot_", names(dat), value = TRUE) # 🏮
	centers <- dat %>% count(center) %>% filter(n >= 100) %>% pull(center)
	pairs <- expand.grid(outcome = Ys, center = centers, stringsAsFactors = FALSE)
	if (file.exists(interim_file)) file.remove(interim_file)
	for (i in seq_len(nrow(pairs))) {
		outcome <- pairs$outcome[i]; center <- pairs$center[i]
		result <- tryCatch( analyze_center(center_id = center, dat = dat, outcome_var = outcome),
			error = function(e) {message("Error in ", outcome, " @ ", center, ": ", e$message); NULL}
		)
		if (!is.null(result)) {results[[length(results) + 1]] <- result}
		if (i %% 500 == 0 || i == nrow(pairs)) {
			message("Saving interim results at record ", i)
			interim_df <- bind_rows(results)
			fwrite(interim_df, interim_file, sep = "\t", append = TRUE, col.names = !file.exists(interim_file))
			results <- list()
		}
	}
	fread(file=interim_file, header=TRUE) %>% arrange(lm_p) %>% head()


# library(paran); paran(dat1, iterations = 100, graph = TRUE, centile = 95)

pca <- prcomp(dat1, center = TRUE, scale. = TRUE); max_pc <- min(40, ncol(pca$x))
	pc40 <- as.data.frame(pca$x[, 1:max_pc]) %>% rename_with(~ paste0("PPC", seq_along(.))) %>% bind_cols(eid = dat$eid, .)
	sdev <- pca$sdev; explained_var <- sdev^2 / sum(sdev^2); cumulative <- summary(pca)$importance[3, 1:40]
	plot(explained_var[1:50], type = "b", xlab = "PC", ylab = "Variance Explained")
	plot(cumulative, type = "b", ylab = "Cumulative variance explained", xlab = "Number of PCs")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#🏮下面是 jupyter notebook 代码	
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import pandas as pd
dat = pd.read_csv("/mnt/d/files/120.txt", sep="\t"); # dat.head()
import keplergl
map = keplergl.KeplerGl(height=500); map # 需要把这个作为最后一行命令
map.add_data(data=dat.copy(), name="house")
# %run 120.config.py
# with open('120.config.py', 'w') as f: f.write('config = {}'.format(map.config))
