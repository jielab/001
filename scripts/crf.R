# 🏮在 LINUX系统里面运行！
# https://github.com/thimotei/CFR_calculation
# install.packages(c("reticulate", "greta", "greta.gp"))
pacman::p_load(data.table, dplyr, tidyr, reshape2, lubridate, reticulate, greta, greta.gp)
#install_greta_deps(); library(greta)
use_condaenv('r-reticulate', required=TRUE) # 需要先 conda_create("r-reticulate"); conda_install("r-reticulate", "python")
#greta::install_tensorflow(method="conda", extra_packages = "tensorflow-probability")

setwd('/mnt/d/scripts/covid')
source('library/scale_cfr_temporal.R')
source('library/run_bayesian_model.R')
source('library/delay_distributions.R')

cfr_baseline <- 1.4
cfr_range <- c(1.2, 1.7)
mean <- 13; median <- 9.1
mu_cfr <- log(median)
sigma_cfr <- sqrt(2*(log(mean) - mu_cfr))
hospitalisation_to_death_truncated <- function(x) { plnorm(x + 1, mu_cfr, sigma_cfr) - plnorm(x, mu_cfr, sigma_cfr)}

# read JHU data
file.case <- '/mnt/d/data/covid/time_series_covid19_confirmed_global.csv' #"https://raw.githubusercontent.com/cssegisanddata/covid-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv"
file.death <- '/mnt/d/data/covid/time_series_covid19_deaths_global.csv' #"https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv"
dat.case <- data.table::fread(file.case) %>% 
    .[, c("Province/State", "Lat", "Long") := NULL] %>%
    melt(., id.vars="Country/Region", variable.name="date", value.name="new_cases") %>% 
    setnames(., old = "Country/Region", new = "country") %>% as.data.table() %>% # 🏮 必须加上这个
    .[, date := mdy(date)] %>%
	.[, new_cases := sum(as.numeric(new_cases)), by = c("date", "country")] %>% unique() %>% 
    .[, new_cases := new_cases - shift(new_cases), by = "country"] %>% 
    .[new_cases < 0, new_cases := 0] %>% .[is.na(new_cases) == FALSE] %>%
	.[order(country, date)] %>% unique()
dat.death <- data.table::fread(file.death) %>% 
    .[, c("Province/State", "Lat", "Long") := NULL] %>%
    melt(., id.vars="Country/Region", variable.name="date", value.name="new_deaths") %>% 
    setnames(., old = "Country/Region", new = "country") %>% as.data.table() %>%
    .[, date := mdy(date)] %>%
	.[, new_deaths := sum(as.numeric(new_deaths)), by = c("date", "country")] %>% unique() %>% 
    .[, new_deaths := new_deaths - shift(new_deaths), by = "country"] %>% 
    .[new_deaths < 0, new_deaths := 0] %>% .[is.na(new_deaths) == FALSE] %>%
	.[order(country, date)] %>% unique()
dat0 <- merge.data.table(dat.case, dat.death, by=c("country", "date"), all=FALSE) 

name.country <- "United Kingdom"
dat <- dat0 %>% filter(country==name.country) %>% dplyr::mutate(date = date - mean)
death_threshold_date <- dat %>% mutate(death_cum_sum = cumsum(new_deaths)) %>% filter(death_cum_sum >= 10) %>% pull(date) %>% min() #date where cumulative deaths passed 10
reporting_spike <- dat %>% filter(new_deaths < 0) %>% pull(date) %>% min() # date where negative death reporting spikes occur
# return adjusted date and reporting estimate
cfr <- scale_cfr_temporal(dat) %>% as_tibble() %>% mutate(reporting_estimate = cfr_baseline/ccfr) %>% 
	mutate(reporting_estimate = pmin(reporting_estimate, 1), country = dat$country, date = dat$date, date_num = as.numeric(dat$date), deaths = dat$new_deaths, cases_known = cum_known_t) %>% 
    filter(date >= death_threshold_date, date < reporting_spike) %>% 
    dplyr::select(country, date, date_num, reporting_estimate, deaths, cases_known)
# run_bayesian_model
prediction <- run_bayesian_model(
	cfr, n_inducing=5, cfr_baseline=cfr_baseline, cfr_range=cfr_range, cfr_trend=NULL, verbose=TRUE
)
ci_poly <- tibble::tibble(x = c(plot_data$date, rev(plot_data$date)), y = c(prediction$upper, rev(prediction$lower)))