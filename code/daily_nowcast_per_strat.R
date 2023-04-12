source("./daily_nowcast_per_strat_fun.R")

dat = read_csv("../hospitalization-nowcast-hub/data-truth/COVID-19/COVID-19_hospitalizations_preprocessed.csv")
now = max(dat$date)

print(paste0("Calculations for day: ", now))

max_delay = 35
per_len = 56

dat_prep = prepare_data(dat=dat, 
                        now = now, 
                        start_date = now-per_len+1, 
                        max_delay = max_delay)

if(!dir.exists(paste0("../results/daily_results/", now)))
  dir.create(paste0("../results/daily_results/", now))  
# Load two STAN models
mod <- cmdstan_model("./stan/mod_per_strat.stan")
mod_red <- cmdstan_model("./stan/mod_per_strat_red.stan")

# German-wide nowcast
fit_ger = mod_red$sample(
  data = list(T=per_len,
              D=max_delay,
              S=1,
              r=array(dat_prep$r_ov, dim=c(per_len, max_delay+1)),
              k_wd_haz = dim(dat_prep$x_haz_wd)[3],
              W_wd = dat_prep$x_haz_wd,
              n_week=per_len/7,
              week_ind = rep(1:(per_len/7), each=7),
              alpha=rep(1, max_delay+1)),
  seed = 1241222,
  chains = 4,
  parallel_chains = 4,
  # Reduced number of (warmup) samples for illustrational purpose, increase!
  iter_warmup = 100,
  iter_sampling = 400,
  refresh = 250,
  max_treedepth = 12,
  adapt_delta=0.95,
  show_messages = FALSE)

fit_ger$save_object(paste0("../results/daily_results/", now, "/", now, "_fit_00+_day.rds"))
rm(fit_ger)

# Nowcast per age-group

res_age = mclapply(seq_along(dat_prep$r_age_list), function(ind, r, names, now) {
  print(names[[ind]])
  # Select model per age-group (here: reduced complexity for all)
  if (names[[ind]] %in% c("00-04", "05-14", "15-34")) {
    mod_use = mod_red  
  } else {
    mod_use = mod_red
  }
  fit = mod_use$sample(
    data = list(T=per_len,
                D=max_delay,
                r=r[[ind]],
                k_wd_haz = dim(dat_prep$x_haz_wd)[3],
                W_wd = dat_prep$x_haz_wd,
                n_week=per_len/7,
                week_ind = rep(1:(per_len/7), each=7),
                alpha=rep(1, max_delay+1)),
    seed = 129,
    chains = 4,
    parallel_chains = 4,
    # Reduced number of (warmup) samples for illustrational purpose, increase!
    iter_warmup = 100,
    iter_sampling = 400,
    refresh = 250,
    max_treedepth = 12,
    adapt_delta=0.95,
    show_messages = FALSE)
  fit$save_object(paste0("../results/daily_results/", now, "/", now, "_fit_",names[[ind]],"_day.rds"))
  rm(fit)
},
r = dat_prep$r_age_list,
names = names(dat_prep$r_age_list),
now=now,
mc.cores = 2)

# Nowcast per federal state

res_fed = mclapply(seq_along(dat_prep$r_fed_list), function(ind, r, names, now) {
  print(names[[ind]])
  # Select model per state (here: reduced complexity for all)
  if (names[[ind]] %in% c("DE-ST", "DE-SN","DE-SL","DE-HH", "DE-MV")) {
    mod_use = mod_red  
  } else {
    mod_use = mod_red
  }
  fit = mod_use$sample(
    data = list(T=per_len,
                D=max_delay,
                r=r[[ind]],
                k_wd_haz = dim(dat_prep$x_haz_wd)[3],
                W_wd = dat_prep$x_haz_wd,
                n_week=per_len/7,
                week_ind = rep(1:(per_len/7), each=7),
                alpha=rep(1, max_delay+1)),
    seed = 129,
    chains = 4,
    parallel_chains = 4,
    # Reduced number of (warmup) samples for illustrational purpose, increase!
    iter_warmup = 100,
    iter_sampling = 400,
    refresh = 250,
    max_treedepth = 12,
    adapt_delta=0.95,
    show_messages = FALSE)
  fit$save_object(paste0("../results/daily_results/", now, "/", now, "_fit_",names[[ind]],"_day.rds"))
  rm(fit)
},
r = dat_prep$r_fed_list,
names = names(dat_prep$r_fed_list),
now=now,
mc.cores = 2)


q(save = "no")