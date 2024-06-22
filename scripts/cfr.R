# 🏮在 LINUX系统里面运行！
# https://github.com/thimotei/CFR_calculation

pacman::p_load(data.table, dplyr, tidyr, reshape2, lubridate, countrycode)
# 必须先 conda create -n r-reticulate python
library(reticulate); use_condaenv('r-reticulate', required=TRUE)
library(greta); library(greta.gp) # 不能跟之前的libraries一起load
# 🏮必须先 pip install tensorflow
greta::install_tensorflow(method="conda", extra_packages="tensorflow-probability")

dir='/work/sph-huangj'
sapply(list.files(path=paste0(dir,'/scripts/covid/library'), pattern="\\.R$", full.names=TRUE), source)

cfr_baseline <- 1.4; cfr_range <- c(1.2, 1.7)
mean <- 13; median <- 9.1; mu_cfr <- log(median)
sigma_cfr <- sqrt(2*(log(mean) - mu_cfr))
hospitalisation_to_death_truncated <- function(x) { plnorm(x + 1, mu_cfr, sigma_cfr) - plnorm(x, mu_cfr, sigma_cfr)}

dat.case <- import_jhu_data(paste0(dir,'/data/covid/time_series_covid19_confirmed_global.csv')) #"https://raw.githubusercontent.com/cssegisanddata/covid-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
dat.death <- import_jhu_data(paste0(dir,'/data/covid/time_series_covid19_deaths_global.csv')) %>% rename(new_deaths=new_cases) #"https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
dat0 <- merge.data.table(dat.case, dat.death, by=c("country", "iso3c", "date"), all=FALSE) 

iso.code <- "GBR"
dat <- dat0 %>% filter(iso3c==iso.code) %>% dplyr::mutate(date = date - mean)
death_threshold_date <- dat %>% mutate(death_cum_sum = cumsum(new_deaths)) %>% filter(death_cum_sum >= 10) %>% pull(date) %>% min() #date where cumulative deaths passed 10
reporting_spike <- dat %>% filter(new_deaths < 0) %>% pull(date) %>% min() # date where negative death reporting spikes occur
# return adjusted date and reporting estimate
cfr <- scale_cfr_temporal(dat) %>% as_tibble() %>% mutate(reporting_estimate = cfr_baseline/ccfr) %>% 
	mutate(reporting_estimate = pmin(reporting_estimate, 1), country = dat$country, date = dat$date, date_num = as.numeric(dat$date), deaths = dat$new_deaths, cases_known = cum_known_t) %>% 
    filter(date >= death_threshold_date, date < reporting_spike) %>% 
    dplyr::select(country, date, date_num, reporting_estimate, deaths, cases_known)
# run_bayesian_model
prediction <- run_bayesian_model(cfr, n_inducing=5, cfr_baseline=cfr_baseline, cfr_range=cfr_range, cfr_trend=NULL, verbose=TRUE)
ci_poly <- tibble::tibble(x = c(plot_data$date, rev(plot_data$date)), y = c(prediction$upper, rev(prediction$lower)))