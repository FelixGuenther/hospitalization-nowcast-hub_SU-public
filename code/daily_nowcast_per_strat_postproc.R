source("./daily_nowcast_per_strat_fun.R")
dat = read_csv("../hospitalization-nowcast-hub/data-truth/COVID-19/COVID-19_hospitalizations_preprocessed.csv")
now = max(dat$date)
max_delay = 35
per_len = 56

if (now==Sys.Date()) {
  files = list.files(paste0("../results/daily_results/", now), full.names = T)
  smry_list = lapply(files, function(x) {summarise_fit(x, 
                                                       now = now, 
                                                       start_date = now-per_len+1, 
                                                       max_delay = max_delay)})
  write_csv(do.call(rbind, smry_list), 
            file = paste0("../hospitalization-nowcast-hub/data-processed/SU-hier_bayes/",now,"-SU-hier_bayes.csv"))
}